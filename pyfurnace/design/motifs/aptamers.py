from . import CONFS_PATH
from ..core.symbols import *
from ..core.coordinates_3d import Coords
from ..core.strand import Strand
from ..core.basepair import BasePair
from ..core.motif import Motif
from .loops import Loop

###
# FLAPS
###

def Ispinach(**kwargs):
    # PDB: 5OB3
    strand1 = Strand("─CUG─UU─GA─GUAGAGUGUGGGCUC─")
    strand1._coords = Coords.load_from_file(CONFS_PATH / "Ispinach_1.dat")

    strand2 = Strand("─GUGAG──GGUCGGG──UC────CAG─", start=(26, 2), direction=(-1, 0))
    strand2._coords = Coords.load_from_file(CONFS_PATH / "Ispinach_2.dat")

    base_pairing = BasePair({(1, 0): (1, 2), (2, 0): (2, 2), (3, 0): (3, 2), (8, 0): (8, 2), (9, 0): (9, 2), (23, 0): (23, 2), (25, 0): (25, 2)})
    kwargs.setdefault('join', False)
    return Motif([strand1, strand2], basepair=base_pairing, **kwargs) 

def Mango(open_left = False, **kwargs):
    # PDB: 5V3F
    strand = Strand("─GUGC─GAA─GG─GAC─GG─UGC╰│╭────GG─AGA─GG─AGA─GCAC─", start=(23, 2), direction=(-1, 0))
    strand._coords = Coords.load_from_file(CONFS_PATH / "Mango.dat")

    return Loop(strands=strand, open_left=open_left, **kwargs)

def MalachiteGreen(open_left = False, **kwargs):
    # PDB: 1Q8N
    strand = Strand("─GGAUCC───CG─A──CUGGCGA╰│╭GAGCCAGGUAACGAAUGGAUCC─", start=(23, 2), direction=(-1, 0))
    strand._coords = Coords.load_from_file(CONFS_PATH / "MalachiteGreen.dat")

    return Loop(strands=strand, open_left=open_left, **kwargs)

def MalachiteGreenShort(**kwargs):
    # PDB: 1Q8N
    strand1 = Strand("CC───CG─A──CUG")
    strand1._coords = Coords.load_from_file(CONFS_PATH / "MalachiteGreenShort_1.dat")

    strand2 = Strand("CAGGUAACGAAUGG", start=(13, 2), direction=(-1, 0))
    strand2._coords = Coords.load_from_file(CONFS_PATH / "MalachiteGreenShort_2.dat")
    kwargs.setdefault('join', False)
    return Motif(strands=[strand1, strand2], **kwargs)
            
def Broccoli(**kwargs):
    # PDB: 7ZJ5
    strand1 = Strand("GGAGAC────GGUCGGG─UC────CAG")
    strand1._coords = Coords.load_from_file(CONFS_PATH / "Broccoli_1.dat")

    strand2 = Strand("CUG─UC─GA─GUAGAGUGUG─GGCUCC", start=(26, 2), direction=(-1, 0))
    strand2._coords = Coords.load_from_file(CONFS_PATH / "Broccoli_2.dat")
    kwargs.setdefault('join', False)
    return Motif([strand1, strand2], **kwargs) 


def Pepper(**kwargs):
    # PDB: 7ZJ5
    strand1 = Strand("UCCC─CAAUCGU─GGCGU─GUCG─GCCUGC")
    strand1._coords = Coords.load_from_file(CONFS_PATH / "Pepper_1.dat")

    strand2 = Strand("GCAGGC─ACUG─GCGCC─────────GGGA", start=(29, 2), direction=(-1, 0))
    strand2._coords = Coords.load_from_file(CONFS_PATH / "Pepper_2.dat")
    kwargs.setdefault('join', False)
    return Motif([strand1, strand2], **kwargs) 
      

###
# Molecule binding aptamers
###

def Biotin(**kwargs):
    # PDB 1F27
    strand = Strand("────GGACCGU─CA───╮│╯GAGGACACGGU─────╭U╰AAAAA─────GUCCUCU")
    strand._coords = Coords.load_from_file(CONFS_PATH / "Biotin.dat")
    kwargs.setdefault('join', False)
    return Motif(strands=strand, **kwargs)

###
# Protein binding aptamers
###

def MS2(open_left = False, **kwargs):
    # PDB: 1ZDH, extended 3 nucleotides on each side
    strand = Strand("─ACAUG─A─GG─AUCA╰│╭─────CC───CAUGU─", start=(16, 2), direction=(-1, 0))
    strand._coords =  Coords.load_from_file(CONFS_PATH / "MS2.dat", 
                                            topology_file=CONFS_PATH / "MS2.top",
                                            protein=True)
    
    return Loop(strands=strand, open_left=open_left, **kwargs)

def PP7(open_left = False, **kwargs):
    # PDB: 2QUX
    strand = Strand("─GGCAC─A─GAAG─AUAUGG╰│╭───────CUUC───GUGCC─", start=(20, 2), direction=(-1, 0))
    strand._coords = Coords.load_from_file(CONFS_PATH / "PP7.dat", 
                                            topology_file=CONFS_PATH / "PP7.top",
                                            protein=True)

    return Loop(strands=strand, open_left=open_left, **kwargs)

def TAR_TAT(open_left = False, **kwargs):
    # PDB: 6MCE
    strand = Strand("─GGGCAGA─UU─GAGC─C─UG╰│╭G──GAGCUC────UCUGCCC─", start=(21, 2), direction=(-1, 0))
    strand._coords = Coords.load_from_file(CONFS_PATH / "TAR_TAT.dat", 
                                            topology_file=CONFS_PATH / "TAR_TAT.top",
                                            protein=True)
    return Loop(strands=strand, open_left=open_left, **kwargs)

def L7Ae(**kwargs):
    # PDB: 1RLG
    strand1 = Strand("─UCU─GA────CC─")
    strand1._coords = Coords.load_from_file(CONFS_PATH / "L7Ae_1.dat")
    strand2 = Strand("─GG─CGUGA─UGA─", start=(13, 2), direction=(-1, 0))
    strand2._coords = Coords.load_from_file(CONFS_PATH / "L7Ae_2.dat",
                                            topology_file=CONFS_PATH / "L7Ae_2.top",
                                            protein=True)
    kwargs.setdefault('join', False)
    return Motif([strand1, strand2], **kwargs)

def Pip3(open_left = False, **kwargs):
    # PDB: none; generate with RNA Composer (https://doi.org/10.1093/nar/gks339, https://doi.org/10.1002/prot.26578)
    # Pubblication: https://doi.org/10.1038/ncb3473
    strand = Strand('GGGUAGACUC╰│╭──────GCUC', start=(10, 2), direction=(-1, 0))
    strand._coords = Coords.load_from_file(CONFS_PATH / "Pip3.dat")
    return Loop(strands=strand, open_left=open_left, **kwargs)
    
def Pip3_mut1(open_left = False, **kwargs):
    # PDB: none; generate with RNA Composer (https://doi.org/10.1093/nar/gks339, https://doi.org/10.1002/prot.26578)
    # Pubblication: https://doi.org/10.1038/ncb3473
    strand = Strand('GGGUCGACUC╰│╭──────GCUC', start=(10, 2), direction=(-1, 0))
    strand._coords = Coords.load_from_file(CONFS_PATH / "Pip3_mut1.dat")
    return Loop(strands=strand, open_left=open_left, **kwargs)

def Pip3_mut3(open_left = False, **kwargs):
    # PDB: none; generate with RNA Composer (https://doi.org/10.1093/nar/gks339, https://doi.org/10.1002/prot.26578)
    # Pubblication: https://doi.org/10.1038/ncb3473
    strand = Strand('GGGUAGCCUC╰│╭──────GCUC', start=(10, 2), direction=(-1, 0))
    strand._coords = Coords.load_from_file(CONFS_PATH / "Pip3_mut3.dat")
    return Loop(strands=strand, open_left=open_left, **kwargs)

def Pip3_mut5(open_left = False, **kwargs):
    # PDB: none; generate with RNA Composer (https://doi.org/10.1093/nar/gks339, https://doi.org/10.1002/prot.26578)
    # Pubblication: https://doi.org/10.1038/ncb3473
    strand = Strand('GGGUAGACGC╰│╭──────GCUC', start=(10, 2), direction=(-1, 0))
    strand._coords = Coords.load_from_file(CONFS_PATH / "Pip3_mut5.dat")
    return Loop(strands=strand, open_left=open_left, **kwargs)

def Streptavidin(open_left = False, **kwargs):
    # PDB: none; generate with RNA Composer (https://doi.org/10.1093/nar/gks339, https://doi.org/10.1002/prot.26578)
    # Pubblication: https://doi.org/10.1093/nar/gkt956
    strand = Strand('AUGCGGCCGCCGACCAGAAUCAUGCAAGUGCGUAAGAUAGU╰│╭─────────CGCG─────────────GGUCGGCGGCCGCAU',
                    start=(41, 2),
                    direction=(-1, 0))
    strand._coords = Coords.load_from_file(CONFS_PATH / "Streptavidin.dat")
    return Loop(strands=strand, open_left=open_left, **kwargs)