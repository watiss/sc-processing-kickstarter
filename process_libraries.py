import scanpy as sc
import shutil
from utils import *
import os

base_path, library_ids = get_top_level_objects()

working_dir = f"{base_path}/processed"
os.system(f"mkdir -p {working_dir}")
logger = start_file_logger(f"{working_dir}/process_libraries.log")

def process_lib(library_id: str):
    adata = sc.read_10x_mtx(f"{base_path}/aligned/{library_id}/filtered_feature_bc_matrix", gex_only = True)

    adata.obs["library_id"] = library_id
    # also add the lib id to the cell barcode names so that each barcode is unique to a cell+lib
    adata.obs_names = adata.obs_names + "_" + library_id
    logger.info(adata)

    logger.info("Filter out poor quality cells and genes...")
    n_cells_before = adata.n_obs
    # TODO Customize min_genes_th for your scenario if needed
    min_genes_th = 300
    sc.pp.filter_cells(adata, min_genes=min_genes_th)
    logger.info(f"Filtered out {n_cells_before-adata.n_obs} cells out of {n_cells_before} that have less than {min_genes_th} genes expressed.")
    logger.info(f"Percent cells filtered: {100*(n_cells_before-adata.n_obs)/n_cells_before:.2f}%.")

    n_cells_before = adata.n_obs
    # TODO Customize filter_cells_min_umi_th for your scenario if needed
    filter_cells_min_umi_th = 1000
    adata = adata[adata.X.sum(axis=1) > filter_cells_min_umi_th].copy()
    logger.info(f"Filtered out {n_cells_before-adata.n_obs} cells out of {n_cells_before} that have less than {filter_cells_min_umi_th} total umi's.")
    logger.info(f"Percent cells filtered: {100*(n_cells_before-adata.n_obs)/n_cells_before:.2f}%.")

    n_genes_before = adata.n_vars
    # TODO Customize min_cells_th for your scenario if needed
    min_cells_th = 0
    sc.pp.filter_genes(adata, min_cells=min_cells_th)
    logger.info(f"Filtered out {n_genes_before-adata.n_vars} genes out of {n_genes_before} that are detected in less than {min_cells_th} cells.")
    logger.info(f"Percent genes filtered: {100*(n_genes_before-adata.n_vars)/n_genes_before:.2f}%.")

    if adata.n_obs == 0:
        raise ValueError("No cells left after basic filtering steps")

    logger.info("Calculate QC metrics...")
    adata.var['mt'] = adata.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'
    adata.var['ribo'] = adata.var_names.str.startswith(("RPS","RPL")) # annotate the group of ribosomal genes genes as 'ribo'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','ribo'], percent_top=None, log1p=False, inplace=True)

    logger.info("Filter out cells based on mito gene concentrations...")
    n_cells_before = adata.n_obs
    # TODO Customize max_pct_counts_mt_th for your scenario if needed
    max_pct_counts_mt_th = 20
    adata = adata[adata.obs['pct_counts_mt'] <= max_pct_counts_mt_th, :].copy()
    logger.info(f"Filtered out {n_cells_before-adata.n_obs} cells with more than {max_pct_counts_mt_th}% counts coming from mitochondrial genes.")
    logger.info(f"Percent cells filtered: {100*(n_cells_before-adata.n_obs)/n_cells_before:.2f}%.")

    if adata.n_obs == 0:
        raise ValueError("No cells left after basic filtering steps")

    # n_genes_by_counts is the number of genes detected in a cell
    # pct_counts_mt is percent of total counts that came from mt- genes
    file_ext = ".svg" 
    sc.pl.violin(
        adata,
        ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo'],
        jitter=0.4,
        multi_panel=True,
        log=True,
        save=file_ext
    )
    figures_dir = f"{working_dir}/figures"
    os.system(f"mkdir -p {figures_dir}")
    file_name = f"{library_id}_statistics{file_ext}" 
    shutil.move(f"figures/violin{file_ext}", f"{figures_dir}/{file_name}")
    sc.pl.violin(adata, ["pct_counts_mt"], save=file_ext)
    file_name = f"{library_id}_pct_counts_mt{file_ext}" 
    shutil.move(f"figures/violin{file_ext}", f"{figures_dir}/{file_name}")

    logger.info("Saving the processed library object")
    adata.write_h5ad(f"{working_dir}/{library_id}_processed.h5ad", compression="lzf")

for library_id in library_ids:
    logger.info(f"Processing library {library_id}...")
    process_lib(library_id)
