import logging
import os

logging.getLogger('numba').setLevel(logging.WARNING)
logging.getLogger('matplotlib').setLevel(logging.WARNING)

def get_top_level_objects():
    # TODO Provide the full path to a local working directory where all processing results (including logs)
    # will be saved.
    base_path = "your_absolute_path_to_processing_results_dir"

    # TODO Provide your library ID's in this list 
    library_ids = ["your_lib_id_1", "your_lib_id_2", "your_lib_id_3", "your_lib_id_4"]

    return base_path, library_ids

app_icon = "⚙️"
def start_file_logger(filename):
    for handler in logging.root.handlers:
        logging.root.removeHandler(handler)
    if os.path.isfile(filename):
        os.remove(filename)
    logging.basicConfig(
        level = logging.NOTSET,
        format = f"%(asctime)s %(levelname)-8s {app_icon} %(message)s",
        datefmt = "[%Y-%m-%d %H:%M:%S]",
        filename = filename,
    )
    return logging.getLogger("sc-mini-pipeline")