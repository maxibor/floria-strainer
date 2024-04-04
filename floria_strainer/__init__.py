from importlib import metadata
import logging

__version__ = metadata.metadata("floria-strainer")["Version"]
__author__ = metadata.metadata("floria-strainer")["Author"]

logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(message)s")
logger = logging.getLogger(__name__)
