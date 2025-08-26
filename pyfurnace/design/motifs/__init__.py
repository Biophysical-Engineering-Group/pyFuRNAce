from pathlib import Path
import inspect

# The configurations are needed to build the motifs
CONFS_PATH = Path(__file__).parent / "conf_files"

# Pyfurnace module import
from ..core.motif import Motif
from .stem import *
from .dovetail import Dovetail
from .kissing_loops import *
from .loops import TetraLoop
from .aptamers import *
from .structural import *

TL = TetraLoop
KD = KissingDimer
KD120 = KissingDimer120
KL180 = KissingLoop180
KL120 = KissingLoop120
BKL = BranchedKissingLoop
BD = BranchedDimer
DT = Dovetail

aptamers_list = [
    func_name
    for func_name, member in inspect.getmembers(aptamers, inspect.isfunction)
    if member.__module__ == aptamers.__name__
]
aptamers_list.remove("create_aptamer")
