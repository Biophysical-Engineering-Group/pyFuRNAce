import numpy as np
from ..core.symbols import *
from ..core.coordinates_3d import Coords
from ..core.sequence import Sequence
from ..core.strand import Strand
from .stem import Stem

### Experimental DAE
# DAE PARALLEL FROM PDB MEAN 

DAE_T_53 = np.array([[-0.79306138, -0.05738916, -0.60643229, -0.60152068],
                    [-0.05738916, -0.98408458,  0.16817856, -0.26217635],
                    [-0.60643229,  0.16817856,  0.77714596, -0.83494375],
                    [ 0.        ,  0.        ,  0.        ,  1.        ],])

DAE_T_35 = np.array([[-0.79306138, -0.05738916, -0.60643229, -0.99842575],
                    [-0.05738916, -0.98408458,  0.16817856, -0.15210483],
                    [-0.60643229,  0.16817856,  0.77714596,  0.32818404],
                    [ 0.        ,  0.        ,  0.        ,  1.        ],])

class Dovetail(Stem):
    
    def __init__(self, length: int = 0 , sequence: str = '', top_cross = True, bot_cross = True, sign: int = 1, wobble_interval: int = 5, wobble_tolerance: int = 2, wobble_insert : str = "middle", **kwargs):
        """
        Attributes of class DoveTail.
        The class DoveTail inherts all attributes from the parentclass Motif.
        -----------------------------------------------------------------------------------------
        length: int
            integer with a sign, 
            the integer stands for the length of nuclotides making up the Stem,
            the sign gives the direction of the dove tail (from left to right or right to left)
        sequence: str
            nucleotide sequence in core of dovetail
        top_dt: bool
            indicates if the top connection of the dovetail should be added
        bot_dt: bool
            indicates if the bottom connection of the dovetail should be added
        sign: int
            indicates the direction of the dovetail, +1 for positive and -1 for negative
        
        """
        if not isinstance(length, int):
            raise ValueError("The length parameter must be an integer.")
        
        # initialize the attributes
        self._top_cross = bool(top_cross)
        self._bot_cross = bool(bot_cross)
        if sequence:
            if sign < 0:
                self._sign = -1
            else:
                self._sign = +1
        else:
            self._sign = +1 if length >= 0 else -1

        super().__init__(length=length, sequence=sequence, wobble_interval=wobble_interval, wobble_tolerance=wobble_tolerance, wobble_insert=wobble_insert, **kwargs)

    ### 
    ### PROPERTIES
    ###

    @property
    def top_cross(self):
        """Returns boolian describing wether the dovetail has a top crossing"""
        return self._top_cross
    
    @top_cross.setter
    def top_cross(self, new_bool):
        """ Set boolian describing wether the dovetail has a top crossing"""
        self._top_cross = bool(new_bool)
        self.length = self._length

    @property
    def bot_cross(self):
        """Returns boolian describing wether the dovetail has a bottom crossing"""
        return self._bot_cross
    
    @bot_cross.setter
    def bot_cross(self, new_bool):
        """Set boolian describing wether the dovetail has a bottom crossing"""
        self._bot_cross = bool(new_bool)
        self.length = self._length

    def set_top_sequence(self, new_seq, sign = 0):
        """ Set the sequence of the top strand """
        if not isinstance(new_seq, (str, Sequence)):
            raise TypeError(f'The sequence of a stem must be a string, got {new_seq}.')
        if sign not in [-1, 0, 1]:
            raise ValueError(f'The sign of the dovetail must be -1, 0 or 1, got {sign}.')
        self._sign = sign
        if not sign:
            if self._length >= 0:
                self._sign = +1
            else:
                self._sign = -1
            
        self._length = len(new_seq) * self._sign
        self._create_strands(sequence=new_seq)

    def set_bottom_sequence(self, new_seq, sign = None):
        """ Set the sequence of the bottom strand """
        self.set_top_sequence(sequence=new_seq.translate(nucl_to_pair)[::-1], sign=sign)

    ###
    ### Protected METHODS
    ###

    def _create_strands(self, sequence: str=None, length: int=0, return_strands=False, **kwargs):
        """Create the strands of the dovetail"""
        # select the direction of the dovetail
        if sequence:
            pos = True if self._sign >= 0 else False
            seq_len = len(sequence)
        else:
            pos = True if length >= 0 else False
            seq_len = abs(length)
            self._sign = +1 if pos else -1

        top_cross = self._top_cross
        bot_cross = self._bot_cross

        kwargs.setdefault('strong_bases', top_cross and bot_cross)
        ### Create the top and bottom strands (according to wobble_insert and wobble_interval)
        top_strand, bot_strand = super()._create_strands(sequence=sequence, length=length, compute_coords=False, return_strands=True, **kwargs)

        ### Positive dovetail
        if pos:
            ### Top strands
            top_strand.strand = '──' + top_strand.strand + '╯' * top_cross + '─' * (not top_cross)
            top_strand1 = top_strand
            top_strand2 = Strand('╰' * top_cross + '─' * (not top_cross), start=(top_strand1.end[0] + 1, 0), direction=(int(not top_cross), int(top_cross)))

            ### Bottom strands
            bot_strand1 = Strand('╮' * bot_cross + '─' * (not bot_cross), start=(0, 2), direction=(-int(not bot_cross), - int(bot_cross)))
            # adjust the stem start position, strand and direction
            bot_strand.strand = '──' + bot_strand.strand + '╭' * bot_cross + '─' * (not bot_cross)
            bot_strand.start = (bot_strand.start[0] + 4, 2)
            bot_strand.direction = (-1, 0)
            bot_strand2 = bot_strand

        ### Negative dovetail
        else:
            ### Top strands
            top_strand1 = Strand('╯' * top_cross + '─' * (not top_cross), start=(0, 0), direction=(1, 0))
            # adjust the stem start position, strand and direction
            top_strand.strand = '╰' * top_cross + '─' * (not top_cross) + top_strand.strand + '──'
            top_strand.start = (1, 0)
            top_strand.direction = (int(not top_cross), int(top_cross))
            top_strand2 = top_strand

            ### Bottom strands
            bot_strand.strand = '╮' * bot_cross + '─' * (not bot_cross) + bot_strand.strand + '──'
            bot_strand.start = (bot_strand.start[0] + 3, 2)
            bot_strand.direction = (- int(not bot_cross), - int(bot_cross))
            bot_strand1 = bot_strand
            bot_strand2 = Strand('─' * (not bot_cross) + '╭' * bot_cross, start=(bot_strand1.start[0] + 1, 2), direction=(-1, 0))

        ### set up the coordinates
        coords = Coords.compute_helix_from_nucl((0,0,0), # start position
                                            (1,0,0), # base vector
                                            (0,1,0), # normal vector
                                            length= seq_len + 2, # length of the helix + dummy ends
                                            double = True)
        # leave out the first and last nucleotide to add the dummy ends
        # Here a schematic of the coordinates indexes:
        #   top_strand1; top_strand2
        #        |           |
        #        0; seq_len;  seq_len + 1;
        #        |        |  |
        #       -N--NNNNNNN--N->
        #        :  :::::::  :
        #      <-N--NNNNNNN--N-
        #        |           |
        #  seq_len * 2 + 3;  seq_len + 2
        #        |           |
        #   bot_strand1; bot_strand2


        if pos: ### the dovetail is positive
            # top strand 1
            top_coord1 = Coords(coords[1: seq_len + 1])
            if top_cross:
                top_coord1.dummy_ends = (coords[0],  # necessary dummy for 0 DT
                                         np.array(Coords.apply_transformation(DAE_T_53, coords[seq_len][0], coords[seq_len][1], coords[seq_len][2], local=True)))
            # top strand 2
            top_coord2 = Coords(np.array(()))
            if top_cross:
                top_coord2.dummy_ends = (np.array(Coords.apply_transformation(DAE_T_35, coords[seq_len + 1][0], coords[seq_len + 1][1], coords[seq_len + 1][2], local=True)),
                                         coords[seq_len + 1])
            # bot strand 1
            bot_coord1 = Coords(np.array(()))
            if bot_cross:
                bot_coord1.dummy_ends = (np.array(Coords.apply_transformation(DAE_T_35, coords[-1][0], coords[-1][1], coords[-1][2], local=True)),
                                         coords[-1])
            # bot strand 2
            bot_coord2 = Coords(coords[seq_len + 3: seq_len * 2 + 3])
            if bot_cross:
                bot_coord2.dummy_ends = (coords[seq_len + 2], # necessary dummy for 0 DT
                                         np.array(Coords.apply_transformation(DAE_T_53, coords[-2][0], coords[-2][1], coords[-2][2], local=True)))
        else: ### the dovetail is negative
            # top strand 1
            top_coord1 = Coords(np.array(()))
            if top_cross:
                top_coord1.dummy_ends = (coords[0],
                                         np.array(Coords.apply_transformation(DAE_T_53, coords[0][0], coords[0][1], coords[0][2], local=True)))
            # top strand 2
            top_coord2 = Coords(coords[1: seq_len + 1])
            if top_cross:
                top_coord2.dummy_ends = (np.array(Coords.apply_transformation(DAE_T_35, coords[1][0], coords[1][1], coords[1][2], local=True)),
                                         coords[seq_len + 1]) # coords[seq_len + 1], useful for Origami ss_assembly
            # bot strand 1
            bot_coord1 = Coords(coords[seq_len + 3: -1])
            if bot_cross:
                bot_coord1.dummy_ends = (np.array(Coords.apply_transformation(DAE_T_35, coords[seq_len + 3][0], coords[seq_len + 3][1], coords[seq_len + 3][2], local=True)),
                                         coords[-1]) #coords[-1], useful for Origami ss_assembly
                
            # bot strand 2
            bot_coord2 = Coords(np.array(()))
            if bot_cross:
                bot_coord2.dummy_ends = (np.array(coords[seq_len + 2]),
                                         np.array(Coords.apply_transformation(DAE_T_53, coords[seq_len + 2][0], coords[seq_len + 2][1], coords[seq_len + 2][2], local=True)))

        top_strand1._coords = top_coord1
        top_strand2._coords = top_coord2
        bot_strand2._coords = bot_coord2
        bot_strand1._coords = bot_coord1

        if return_strands:
            return top_strand1, top_strand2, bot_strand1, bot_strand2

        self.replace_all_strands([top_strand1, top_strand2, bot_strand1, bot_strand2], copy=False, join=True)


    def copy(self, **kwargs):
        return Dovetail(  # Stem parameters
                        sequence=self[0].sequence if self[0].sequence else self[1].sequence, wobble_interval = self.wobble_interval, wobble_insert = self._wobble_insert, 
                            # DT parameters
                        top_cross = self.top_cross, bot_cross = self.bot_cross, sign = int(self.length >= 0) - int(self.length <= 0), 
                            # Copied DT parameters
                        strands = self.copy_strands_preserve_blocks(self._strands),
                        #  Motif parameters
                        basepair = self._basepair, autopairing = self.autopairing,
                        **kwargs)