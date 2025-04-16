from . import CONFS_PATH
from ..core.symbols import *
from ..core.coordinates_3d import Coords
from ..core.strand import Strand
from ..core.motif import Motif


class Loop(Motif):
    def __init__(self, open_left = False, sequence: str = "", **kwargs):
        # create motif of Cap without basepairs by turing autobasepairing of
        if sequence:
            seq_len = len(sequence)
            ### create the strand
            strand = Strand('─' * seq_len + '╰│╭' + sequence, start=(seq_len, 2), direction=(-1, 0))
            # Add the strand to the list of strands
            kwargs['strands'] = kwargs.get('strands', []) + [strand]
        
        kwargs['join'] = False
        super().__init__(**kwargs)
        if open_left:
            self.flip(horizontally=True, vertically=True)


class TetraLoop(Loop):
    
    def __init__(self, open_left = False, sequence = "UUCG", **kwargs):
        """
        Attributes of the class Cap_UUCG, which is a daugther class of the class Motif.
        -----------------------------------------------------------------------------------
        UUCG_bool: bool (default= False)
            indicates if a UUCG sequence should be added into the cap
        """
        #create strands deascribing tetraloop
        if len(sequence) != 4:
            raise ValueError("The sequence length doesn't match the length required for a tetraloop, which is 4.")
        
        # Create new strands if the strand is not provided
        if 'strands' in kwargs:
            strands = kwargs.pop('strands')
        else:
            strand = Strand(sequence[:2] + "╰│╭" + sequence[2:4], start=(2,2), direction=(-1,0))
            ### PDB: 2KOC
            strand._coords = Coords.load_from_file(CONFS_PATH / 'TetraLoop.dat', 
                                                   dummy_ends=(True, True))
            strands = [strand]

        # create motif of Cap without basepairs by turing autobasepairing of
        kwargs.setdefault('autopairing', False)

        super().__init__(strands=strands, open_left=open_left, **kwargs)

    def set_sequence(self, new_sequence):
        """
        Set the sequence of the tetraloop
        """
        if len(new_sequence) != 4:
            raise ValueError("The sequence length doesn't match the length required for a tetraloop, which is 4.")
        self[0].sequence = new_sequence
                    

