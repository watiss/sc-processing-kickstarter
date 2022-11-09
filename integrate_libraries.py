import scanpy as sc
import anndata
from anndata import AnnData
import pandas as pd
import numpy as np
import scvi
from utils import *
import os

base_path, library_ids = get_top_level_objects()

working_dir = f"{base_path}/integrated"
os.system(f"mkdir -p {working_dir}")
logger = start_file_logger(f"{working_dir}/integrate_libraries.log")

logger.info("Reading all processed library objects")
adata_dict = {}
for library_id in library_ids:
    adata_dict[library_id] = anndata.read_h5ad(f"{base_path}/processed/{library_id}_processed.h5ad")

logger.info("Determining the set of genes for solo")
solo_genes = set()
# TODO Customize solo_genes_min_cells_th for your scenario if needed
solo_genes_min_cells_th = 5
for library_id in library_ids:
    solo_genes_idx = sc.pp.filter_genes(adata_dict[library_id], min_cells=solo_genes_min_cells_th, inplace=False)[0]
    solo_genes.update(adata_dict[library_id].var.index[solo_genes_idx].tolist())

logger.info("Concatenating all processed libraries")
adata = adata_dict[library_ids[0]].concatenate(
    [adata_dict[library_ids[j]] for j in range(1, len(library_ids))],
    join="outer",
    index_unique=None,
)
# there should be no duplicate barcodes
assert len(np.unique(adata.obs.index)) == adata.n_obs

logger.info(f"Concatenated adata object has {adata.n_obs} cells and {adata.n_vars} gene.")

# Note: If you want to estimate contamination levels from ambient RNA in your data, you
# may use decontX (DOI 10.1186/s13059-020-1950-6). This pipeline does not include that.

logger.info("Detecting highly variable genes...")
# TODO Consider if you should use other batch keys for your dataset
batch_key = "library_id"
# TODO Customize n_top_genes for your scenario if needed
n_top_genes = 4000
sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, subset=True, flavor="seurat_v3", batch_key=batch_key, span = 1.0)

def run_model(
        adata: AnnData,
        batch_key: str,
        model_save_path: str,
    ) -> str: 
    logger.info(f"Setting up adata...")
    scvi.model.SCVI.setup_anndata(adata, batch_key=batch_key)
    logger.info(f"Training SCVI model with batch key {batch_key}...")
    model = scvi.model.SCVI(adata)
    train_params = {}
    # this is justified by an experiment described here https://yoseflab.github.io/scvi-tools-reproducibility/runtime_analysis_reproducibility/
    train_params["train_size"] = 0.9 if 0.1 * adata.n_obs < 20000 else 1-(20000/adata.n_obs)
    model.train(**train_params)
    logger.info("Number of training epochs: {}...".format(len(model.history_["elbo_train"])))
    logger.info("Saving SCVI latent representation...")
    latent = model.get_latent_representation()
    adata.obsm["X_scVI"] = latent
    logger.info("Saving the model...")
    model_file_path = os.path.join(model_save_path, f"scvi_model_batch_key_{batch_key}")
    model.save(model_file_path, overwrite=True, save_anndata=True)
    return model, model_file_path

logger.info("Running scvi...")
scvi_model, _ = run_model(adata, batch_key, f"{working_dir}")

logger.info("Running solo for detecting doublets...")
batches = pd.unique(adata.obs[batch_key])
logger.info("Running solo on the following batches separately: {}".format(batches))
is_solo_singlet = np.ones((adata.n_obs,), dtype=bool)
solo_max_epochs = np.min([round((20000 / adata.n_obs) * 400), 400]) # SOLO currently has a bug where it doesn't run this heuristic
logger.info(f"Solo max_epochs: {solo_max_epochs}.")
for batch in batches:
    logger.info("Running solo on batch {}...".format(batch))
    solo_batch = scvi.external.SOLO.from_scvi_model(scvi_model, restrict_to_batch=batch)
    solo_batch.train(max_epochs=solo_max_epochs)
    is_solo_singlet[(adata.obs["library_id"] == batch).values] = solo_batch.predict(soft=False) == "singlet"
    adata.obs["is_solo_singlet"] = is_solo_singlet

logger.info("Removing doublets...")
n_obs_before = adata.n_obs
adata = adata[is_solo_singlet,].copy()
percent_removed = 100*(n_obs_before-adata.n_obs)/n_obs_before
logger.info("Removed {} estimated doublets (percent removed: {:.2f}%); {} droplets remained.".format(
    n_obs_before-adata.n_obs,
    percent_removed,
    adata.n_obs)
)

if adata.n_obs == 0:
    raise ValueError("No cells left after doublet detection")

logger.info("Calculating neighborhood graph and UMAP based on SCVI components...")
key = "scvi_neighbors"
sc.pp.neighbors(adata, use_rep="X_scVI", key_added=key) 
adata.obsm["X_umap_scvi"] = sc.tl.umap(adata, neighbors_key=key, copy=True).obsm["X_umap"]
sc.tl.leiden(adata, neighbors_key=key, key_added="scvi_leiden")

logger.info("Normalizing data...")
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=10000)
sc.pp.log1p(adata)
sc.pp.scale(adata)
adata.raw = adata

logger.info("Calculating PCA...")
sc.tl.pca(adata)
logger.info("Calculating neighborhood graph and UMAP based on PCA...")
key = "pca_neighbors"
sc.pp.neighbors(adata, use_rep="X_pca", key_added=key) 
adata.obsm["X_umap_pca"] = sc.tl.umap(adata, neighbors_key=key, copy=True).obsm["X_umap"]
sc.tl.leiden(adata, neighbors_key=key, key_added="pca_leiden")

logger.info(f"✅ Final adata:\n{adata}")

logger.info("Saving final integrated data object")
adata.write_h5ad(f"{working_dir}/integrated_adata.h5ad", compression="lzf")
