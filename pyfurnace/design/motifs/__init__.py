from pathlib import Path
CONFS_PATH = Path(__file__).parent / 'conf_files'

from ..core.motif import Motif
from .stem import *
from .dovetail import Dovetail
from .kissing_loops import *
from .loops import TetraLoop
from .aptamers import *
from .structural import *

TL = TetraLoop
KD = KissingDimer
KL180 = KissingLoop180
KL120 = KissingLoop120
BKL = BranchedKissingLoop
BD = BranchedDimer
DT = Dovetail

