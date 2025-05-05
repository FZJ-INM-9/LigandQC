import logging
from rich.logging import RichHandler

logger = logging.getLogger('ligand_inspector')
if logger.hasHandlers():
    logger.handlers.clear()
logger.propagate = False
logger.setLevel(logging.DEBUG)
file_handler = RichHandler()
file_handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(message)s')
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)