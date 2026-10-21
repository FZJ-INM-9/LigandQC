"""
log
===

This module provides logging functionality for the ligandqc package. It sets up a logger that is used throughout the package to log messages at various levels (DEBUG, INFO, WARNING, ERROR, CRITICAL). The logger is configured to output messages to both the console and a log file, with customizable logging levels.

"""

import logging

import click
from rich.logging import RichHandler

logger = logging.getLogger("ligandqc")
if logger.hasHandlers():
    logger.handlers.clear()
logger.propagate = False
logger.setLevel(logging.DEBUG)
file_handler = RichHandler()
file_handler.setLevel(logging.INFO)
formatter = logging.Formatter("%(message)s")
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)


def setup_logger(log_file, logging_level):
    # set logging level
    logging_level = logging_level.upper()
    if logging_level not in ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]:
        click.echo("Invalid logging level. Please use DEBUG, INFO, WARNING, ERROR, or CRITICAL.")
        return
    # set logging level
    logger.setLevel(logging_level)
    # set log file
    handler = logging.FileHandler(log_file)
    handler.setLevel(logging_level)
    formatter = logging.Formatter("[%(asctime)s] %(levelname)s %(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)
