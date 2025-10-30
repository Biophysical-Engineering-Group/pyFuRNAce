from importlib.metadata import version, PackageNotFoundError
from .design import *  # noqa: F403

__all__ = ["motifs", "strand", "origami", "symbols"]  # noqa: F405

try:
    __version__ = version("pyfurnace")
except PackageNotFoundError:
    # package is not installed (e.g., running from source without build)
    __version__ = "0+unknown"
