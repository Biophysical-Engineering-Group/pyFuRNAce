from . import CONFS_PATH
from ..core.symbols import *
from ..core.coordinates_3d import Coords
from ..core.strand import Strand
from ..core.motif import Motif
from .loops import Loop

# Shared functionality for all Aptamers
class Aptamer(Motif):
    """
    Base class for all aptamers, inheriting from the Motif class.
    This class is used to create aptamers with specific sequences and coordinates.
    At the moment, it does not add any additional functionality to the Motif class, 
    but it can be used to create a common interface for all aptamers or for screening
    motifs in an origami.
    """
    pass

# Create the Aptamer class also inheriting from the selected base class
def create_aptamer(*args, inherit_from: classmethod = None, **kwargs):
    """
    Create an Aptamer object, optionally inheriting from a specified base class.

    Parameters
    ----------
    args : list
        Positional arguments to be passed to the Aptamer constructor.
    inherit_from : classmethod, optional
        A class from which the Aptamer class should inherit. If not provided, 
        the Aptamer class will be used as is.
    kwargs : dict
        Keyword arguments to be passed to the Aptamer constructor.
    
    Returns
    -------
    Aptamer
        An instance of the Aptamer class, optionally inheriting from the specified
          inherit_from class.
    """
    aptamer = Aptamer(*args, **kwargs)
    if inherit_from:
        aptamer = type("Aptamer", (Aptamer, inherit_from), {})(*args, **kwargs)
    else:
        aptamer = Aptamer(*args, **kwargs)
    return aptamer

###
# FLAPS
###

def Ispinach(**kwargs):
    # PDB: 5OB3
    strand1 = Strand("─CUG─UU─GA─GUAGAGUGUGGGCUC─")
    strand1._coords = Coords.load_from_file(CONFS_PATH / "Ispinach_1.dat")

    strand2 = Strand("─GUGAGG─GU─CGG─G─UC────CAG─", start=(26, 2), direction=(-1, 0))
    strand2._coords = Coords.load_from_file(CONFS_PATH / "Ispinach_2.dat")

    kwargs.setdefault('join', False)
    return create_aptamer([strand1, strand2], 
                          **kwargs)

def Mango(open_left = False, **kwargs):
    # PDB: 5V3F
    strand = Strand("─GUGC─GAA─GG─GAC─GG─UGC╰│╭────GG─AGA─GG─AGA─GCAC─", start=(23, 2), direction=(-1, 0))
    strand._coords = Coords.load_from_file(CONFS_PATH / "Mango.dat")
    kwargs['inherit_from'] = Loop
    return create_aptamer(strands=strand, 
                          open_left=open_left, 
                          **kwargs)

def MalachiteGreen(open_left = False, **kwargs):
    # PDB: 1Q8N
    strand = Strand("─GGAUCC───CG─A──CUGGCGA╰│╭GAGCCAGGUAACGAAUGGAUCC─", start=(23, 2), direction=(-1, 0))
    strand._coords = Coords.load_from_file(CONFS_PATH / "MalachiteGreen.dat")
    kwargs['inherit_from'] = Loop
    return create_aptamer(strands=strand,
                          open_left=open_left,
                          **kwargs)

def MalachiteGreenShort(**kwargs):
    # PDB: 1Q8N
    strand1 = Strand("CC───CG─A──CUG")
    strand1._coords = Coords.load_from_file(CONFS_PATH / "MalachiteGreenShort_1.dat")

    strand2 = Strand("CAGGUAACGAAUGG", start=(13, 2), direction=(-1, 0))
    strand2._coords = Coords.load_from_file(CONFS_PATH / "MalachiteGreenShort_2.dat")
    kwargs.setdefault('join', False)
    return create_aptamer(strands=[strand1, strand2], 
                          **kwargs)
            
def Broccoli(**kwargs):
    # PDB: 7ZJ5
    strand1 = Strand("GGAGAC────GGUCGGG─UC────CAG")
    strand1._coords = Coords.load_from_file(CONFS_PATH / "Broccoli_1.dat")

    strand2 = Strand("CUG─UC─GA─GUAGAGUGUG─GGCUCC", start=(26, 2), direction=(-1, 0))
    strand2._coords = Coords.load_from_file(CONFS_PATH / "Broccoli_2.dat")
    kwargs.setdefault('join', False)
    return create_aptamer([strand1, strand2], 
                          **kwargs) 


def Pepper(**kwargs):
    # PDB: 7ZJ5
    strand1 = Strand("UCCC─CAAUCGU─GGCGU─GUCG─GCCUGC")
    strand1._coords = Coords.load_from_file(CONFS_PATH / "Pepper_1.dat")

    strand2 = Strand("GCAGGC─ACUG─GCGCC─────────GGGA", start=(29, 2), direction=(-1, 0))
    strand2._coords = Coords.load_from_file(CONFS_PATH / "Pepper_2.dat")
    kwargs.setdefault('join', False)
    return create_aptamer([strand1, strand2], 
                          **kwargs) 
      

###
# Substrate binding aptamers
###

def Biotin(**kwargs):
    # PDB 1F27
    strand = Strand("────GGACCGU─CA───╮│╯GAGGACACGGU─────╭U╰AAAAA─────GUCCUCU")
    strand._coords = Coords.load_from_file(CONFS_PATH / "Biotin.dat")
    kwargs.setdefault('join', False)
    return create_aptamer(strands=strand, 
                          **kwargs)

def PIP3(open_left = False, **kwargs):
    # PDB: none; generate with RNA Composer (https://doi.org/10.1093/nar/gks339, https://doi.org/10.1002/prot.26578)
    # Pubblication: https://doi.org/10.1038/ncb3473
    strand = Strand('GGGUAGACUC╰│╭──────GCUC', start=(10, 2), direction=(-1, 0))
    strand._coords = Coords.load_from_file(CONFS_PATH / "PIP3.dat")
    kwargs['inherit_from'] = Loop
    return create_aptamer(strands=strand,
                          open_left=open_left,
                          **kwargs)
    
def PIP3_mut1(open_left = False, **kwargs):
    # PDB: none; generate with RNA Composer (https://doi.org/10.1093/nar/gks339, https://doi.org/10.1002/prot.26578)
    # Pubblication: https://doi.org/10.1038/ncb3473
    strand = Strand('GGGUCGACUC╰│╭──────GCUC', start=(10, 2), direction=(-1, 0))
    strand._coords = Coords.load_from_file(CONFS_PATH / "PIP3_mut1.dat")
    kwargs['inherit_from'] = Loop
    return create_aptamer(strands=strand,
                          open_left=open_left,
                          **kwargs)

def PIP3_mut3(open_left = False, **kwargs):
    # PDB: none; generate with RNA Composer (https://doi.org/10.1093/nar/gks339, https://doi.org/10.1002/prot.26578)
    # Pubblication: https://doi.org/10.1038/ncb3473
    strand = Strand('GGGUAGCCUC╰│╭──────GCUC', start=(10, 2), direction=(-1, 0))
    strand._coords = Coords.load_from_file(CONFS_PATH / "PIP3_mut3.dat")
    kwargs['inherit_from'] = Loop
    return create_aptamer(strands=strand,
                          open_left=open_left,
                          **kwargs)

def PIP3_mut5(open_left = False, **kwargs):
    # PDB: none; generate with RNA Composer (https://doi.org/10.1093/nar/gks339, https://doi.org/10.1002/prot.26578)
    # Pubblication: https://doi.org/10.1038/ncb3473
    strand = Strand('GGGUAGACGC╰│╭──────GCUC', start=(10, 2), direction=(-1, 0))
    strand._coords = Coords.load_from_file(CONFS_PATH / "PIP3_mut5.dat")
    kwargs['inherit_from'] = Loop
    return create_aptamer(strands=strand,
                          open_left=open_left,
                          **kwargs)

###
# Protein binding aptamers
###

def MS2(open_left = False, **kwargs):
    # PDB: 1ZDH, extended 3 nucleotides on each side
    strand = Strand("─ACAUG─A─GG─AUCA╰│╭─────CC───CAUGU─", start=(16, 2), direction=(-1, 0))
    strand._coords =  Coords.load_from_file(CONFS_PATH / "MS2.dat", 
                                            topology_file=CONFS_PATH / "MS2.top",
                                            protein=True)
    kwargs['inherit_from'] = Loop
    return create_aptamer(strands=strand, 
                          open_left=open_left, 
                          **kwargs)

def PP7(open_left = False, **kwargs):
    # PDB: 2QUX
    strand = Strand("─GGCAC─A─GAAG─AUAUGG╰│╭───────CUUC───GUGCC─", start=(20, 2), direction=(-1, 0))
    strand._coords = Coords.load_from_file(CONFS_PATH / "PP7.dat", 
                                            topology_file=CONFS_PATH / "PP7.top",
                                            protein=True)
    kwargs['inherit_from'] = Loop
    return create_aptamer(strands=strand,
                          open_left=open_left,
                          **kwargs)

def TAR_TAT(open_left = False, **kwargs):
    # PDB: 6MCE
    strand = Strand("─GGGCAGA─UU─GAGC─C─UG╰│╭G──GAGCUC────UCUGCCC─", start=(21, 2), direction=(-1, 0))
    strand._coords = Coords.load_from_file(CONFS_PATH / "TAR_TAT.dat", 
                                            topology_file=CONFS_PATH / "TAR_TAT.top",
                                            protein=True)
    kwargs['inherit_from'] = Loop
    return create_aptamer(strands=strand,
                          open_left=open_left,
                          **kwargs)

def L7Ae(**kwargs):
    # PDB: 1RLG
    strand1 = Strand("─UCU─GA────CC─")
    strand1._coords = Coords.load_from_file(CONFS_PATH / "L7Ae_1.dat")
    strand2 = Strand("─GG─CGUGA─UGA─", start=(13, 2), direction=(-1, 0))
    strand2._coords = Coords.load_from_file(CONFS_PATH / "L7Ae_2.dat",
                                            topology_file=CONFS_PATH / "L7Ae_2.top",
                                            protein=True)
    kwargs.setdefault('join', False)
    return create_aptamer([strand1, strand2], 
                          **kwargs)

def Streptavidin(open_left = False, **kwargs):
    # PDB: none; generate with RNA Composer (https://doi.org/10.1093/nar/gks339, https://doi.org/10.1002/prot.26578)
    # Pubblication: https://doi.org/10.1093/nar/gkt956
    strand = Strand('AUGCGGCCGCCGACCAGAAUCAUGCAAGUGCGUAAGAUAGU╰│╭─────────CGCG─────────────GGUCGGCGGCCGCAU',
                    start=(41, 2),
                    direction=(-1, 0))
    strand._coords = Coords.load_from_file(CONFS_PATH / "Streptavidin.dat")
    kwargs['inherit_from'] = Loop
    return create_aptamer(strands=strand,
                          open_left=open_left,
                          **kwargs)

def Thrombin_exosite1(open_left = False, **kwargs):
    # PDB: none; generate with RNA Composer (https://doi.org/10.1093/nar/gks339, https://doi.org/10.1002/prot.26578)
    # Pubblication: https://doi.org/10.1111/j.1538-7836.2012.04679.x
    strand = Strand("─GGCG─GUCGAUC─ACACA╭╰GUUC╮A╯A─AC─╰││╮GUAA╰U╭AA──╯││╭─GCCAAUG╮U╯ACGAGGC╭╰────A─GA─CGACUCGCC─",
                    start=(26, 5),
                    direction=(-1, 0))
    strand._coords = Coords.load_from_file(CONFS_PATH / "Thrombin_exosite1.dat")
    
    kwargs['inherit_from'] = Loop
    return create_aptamer(strands=strand,
                          open_left=open_left,
                          **kwargs)

R9D_14T = Thrombin_exosite1

def Thrombin_exosite2(open_left = False, **kwargs):
    # PDB: 5DO4 (replace Fluorinated ribose with regular ribose)
    strand = Strand("-GGG-AACA-AAG-CUGA\\|/AGUA-CUU----A-CCC-",
                    start=(18, 2),
                    direction=(-1, 0))
    strand._coords = Coords.load_from_file(CONFS_PATH / "Thrombin_exosite2.dat",
                                            topology_file=CONFS_PATH / "Thrombin_exosite2.top",
                                            protein=True)
    kwargs['inherit_from'] = Loop
    return create_aptamer(strands=strand,
                          open_left=open_left,
                          **kwargs)