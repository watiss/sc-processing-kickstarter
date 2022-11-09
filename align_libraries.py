import os
import argparse
from utils import *

cmd_desc = "Align the given library using cellranger.\n" \
        "After having downloaded fastq files for all of your libraries, place all fastq files that belong\n" \
        "to a library in a directory named after that library's ID. Give the full path of the latter directory\n" \
        "when calling this script."
parser = argparse.ArgumentParser(
    description=cmd_desc,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
parser.add_argument('lib_fastqs_path', type=str,
                    help='absolute path to the directory containing fastq files for the library to align')
args = parser.parse_args()
lib_fastqs_path = args.lib_fastqs_path
assert os.path.isdir(lib_fastqs_path)

base_path, _ = get_top_level_objects()
working_dir = f"{base_path}/aligned"
os.system(f"mkdir -p {working_dir}")
logger = start_file_logger(f"{working_dir}/align_libraries.log")

# TODO Provide the absolute (full) path to the cellranger executable,
# For instance it might be under cellranger-6.0.0/bin, depending on your cellranger version.
# For how to install cellranger: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation
cellranger = "your_absolute_path_to_cellranger_executable"

# TODO 1. Identify the appropriate genome reference depending on your organism (human, mouse, etc.)
# 2. Download it from cellranger: https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest
# 3. Provide the full path to the reference directory here (this is the dir that contains reference.json)
genome_reference_file = "your_absolute_path_to_cellranger_reference_directory"

# TODO This is "auto" by default, which is the recommended mode.
# If you wish to change it, see the cellranger instructions (link below) for which
# value you should select.
chem = "auto"

# TODO If you know the expected number of recovered cells in your experiment
# (i.e. average number of cells sequenced) you can let cellranger know about it
# via this argument. If you choose to use this, add "--expect-cells={expect_cells}"
# in the "cellranger count" command below.
# expect_cells = "6000"

# get the file name
lib_name = os.path.basename('lib_fastqs_path')

# TODO Feel free to adjust any of the arguments to the "cellranger count" command as you see fit
# Documentation here: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count#args
alignment_cmd = f"{cellranger} count --id={lib_name} \
                                     --transcriptome={genome_reference_file} \
                                     --fastqs={lib_fastqs_path} \
                                     --chemistry={chem}"

logger.info("Running the following alignment command...\n{}".format(alignment_cmd))
alignment_output = os.popen(alignment_cmd).read()
logger.info(alignment_output)
