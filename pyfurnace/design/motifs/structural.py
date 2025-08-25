from . import CONFS_PATH
from ..core.coordinates_3d import Coords
from ..core.strand import Strand
from ..core.motif import Motif


def LambdaTurn(**kwargs):
    """
    Returns a Lambda turn RNA motif.

    A sharp RNA bend motif from PDB 2AW7.
    DOI: https://doi.org/10.1007/978-3-319-22186-1_45
    """
    strand1 = Strand("─CUNGAUGG─")
    strand1._coords = Coords.load_from_file(
        CONFS_PATH / "LambdaTurn_1.dat",
        dummy_ends=(True, True),
    )
    strand2 = Strand("─CCAU──GG─", start=(9, 2), direction=(-1, 0))
    strand2._coords = Coords.load_from_file(
        CONFS_PATH / "LambdaTurn_2.dat",
        dummy_ends=(True, True),
    )
    kwargs["join"] = False
    return Motif([strand1, strand2], **kwargs)


def ThreeWayJunction(**kwargs):
    """
    Returns a 3-way junction RNA motif.

    Based on PDB 2AW4, DOI: https://doi.org/10.1021/nl900261h
    """
    strand1 = Strand("─NC────UAAN─")
    strand1._coords = Coords.load_from_file(
        CONFS_PATH / "ThreeWayJunction_1.dat",
        dummy_ends=(True, True),
    )
    strand2 = Strand("─NG─AC╭A╯╭", start=(11, 2), direction=(-1, 0))
    strand2._coords = Coords.load_from_file(
        CONFS_PATH / "ThreeWayJunction_2.dat",
        dummy_ends=(True, True),
    )
    strand3 = Strand("│U╮GN─", start=(3, 4), direction=(0, -1))
    strand3._coords = Coords.load_from_file(
        CONFS_PATH / "ThreeWayJunction_3.dat",
        dummy_ends=(True, True),
    )
    kwargs["join"] = False
    return Motif([strand1, strand2, strand3], **kwargs)


# """
# 5-NC----UAAN-
#   ::    :  :
#  -NG\ /CA-GN-5
#     U=A
#     |//
#     5
# """


def Bend90(**kwargs):
    """
    Returns a 90-degree RNA bend motif.

    Based on PDB 3P59, DOI: https://doi.org/10.1073/pnas.1101130108
    """
    strand1 = Strand("─GAACUAC─")
    strand1._coords = Coords.load_from_file(
        CONFS_PATH / "Bend90.dat",
        dummy_ends=(True, True),
    )
    strand2 = Strand("─G─────C─", start=(8, 2), direction=(-1, 0))
    strand2._coords = Coords.load_from_file(
        CONFS_PATH / "Bend90_2.dat",
        dummy_ends=(True, True),
    )
    kwargs["join"] = False
    return Motif([strand1, strand2], **kwargs)
