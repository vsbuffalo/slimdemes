"""SLiM and demes integration package."""
from .main import convert as convert_demes_to_json
from .main import run_msprime as run_demes_msprime

from ._version import __version__

__all__ = ["convert_demes_to_json", "run_demes_msprime", "__version__"]
