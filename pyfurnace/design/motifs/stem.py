import random
from ..core.symbols import *
from ..core.sequence import Sequence
from ..core.coordinates_3d import Coords
from ..core.strand import Strand
from ..core.motif import Motif


class Stem(Motif):
    
    def __init__(self, length: int = 0, sequence: str = "", wobble_interval: int = 4, wobble_tolerance: int = 2, wobble_insert : str = "middle", strong_bases=True, **kwargs):
        """
        Attributes of class Stem.
        The class Stem inherts all attributes from the parentclass Motif.
        -------------------------------------------------------------------------
        length: int
            the number of nucleotides in the Stem
        seq: str
            sequence of top strand of the Stem
        wobble_interval: int
            Number of nucleotides between wobble base pairing (wobbles cannot be consecutive, or at the end/start of a sequence), default is 7
        wobble_tolerance: int
            tolerance to randomize the wobble base pairing frequency
        wobble_insert: str
            position where the wobble base pairing is inserted (middle, start, end)
        """
        ### set default values
        if wobble_insert not in ["middle", "start", "end"]:
            raise ValueError(f"Invalid value for wobble_insert: {wobble_insert}. The value must be 'middle', 'start' or 'end'.")
        if not isinstance(wobble_interval, int) or wobble_interval < 0:
            raise TypeError(f'The wobble frequency must be a positive integer, got {wobble_interval}.')
        if not isinstance(wobble_tolerance, int) or wobble_tolerance < 0:
            raise TypeError(f'The wobble tolerance must be a positive integer, got {wobble_tolerance}.')
        if not isinstance(length, int):
            raise TypeError(f'The length of a stem must be an integer, got {length}.')
        if not isinstance(sequence, (str, Sequence)):
            raise TypeError(f'The sequence of a stem must be a string or a Sequence object, got {type(sequence)}.')

        self._wobble_interval = wobble_interval if not sequence else 0
        self._wobble_tolerance = wobble_tolerance if not sequence else 0
        self._wobble_insert = wobble_insert
        self._length = len(sequence) * getattr(self, '_sign', 1) if sequence else length
        
        ### If the user doesn't provide strands, update them directly
        if "strands" in kwargs:
            strands = kwargs.pop("strands")
        else:
            ### create the strands
            strands = self._create_strands(sequence=sequence, length=length, return_strands=True, strong_bases=strong_bases)

        kwargs["join"] = False
        # Initialize the motif
        super().__init__(strands=strands, **kwargs)

    ### 
    ### PROPERTIES
    ###

    @property
    def length(self):
        """ Number of nucleotides in a stem """
        return self._length
    
    @length.setter
    def length(self, new_length):
        if not isinstance(new_length, int):
            raise TypeError(f'The length of a stem must be an integer, got {new_length}.')
        self._create_strands(length=new_length)

    @property
    def wobble_interval(self):
        return self._wobble_interval
    
    @wobble_interval.setter
    def wobble_interval(self, new_freq):
        if not isinstance(new_freq, int) or new_freq < 0:
            raise TypeError(f'The wobble frequency must be a positive integer, got {new_freq}.')
        self._wobble_interval = new_freq
        # update the sequence of the top strand and the bottom strand
        self.length = self._length

    @property
    def wobble_tolerance(self):
        return self._wobble_tolerance

    @wobble_tolerance.setter
    def wobble_tolerance(self, new_tolerance):
        if not isinstance(new_tolerance, int) or new_tolerance < 0:
            raise TypeError(f'The wobble tolerance must be a positive integer, got {new_tolerance}.')
        self._wobble_tolerance = new_tolerance
        # update the sequence of the top strand and the bottom strand
        self.length = self._length

    @property
    def wobble_insert(self):
        return self._wobble_insert
    
    @wobble_insert.setter
    def wobble_insert(self, new_insert):
        if new_insert not in ["middle", "start", "end"]:
            raise ValueError(f"Invalid value for wobble_insert: {new_insert}. The value must be 'middle', 'start' or 'end'.")
        self._wobble_insert = new_insert
        # update the sequence of the top strand and the bottom strand
        self.length = self._length
    
    ###
    ### METHOD
    ###

    def set_top_sequence(self, new_seq):
        """ Set the sequence of the top strand """
        if not isinstance(new_seq, str):
            raise TypeError(f'The sequence of a stem must be a string, got {new_seq}.')
        self._create_strands(sequence=new_seq)

    def set_bottom_sequence(self, new_seq):
        """ Set the sequence of the bottom strand """
        self.set_top_sequence(sequence=new_seq.translate(nucl_to_pair)[::-1])

    def set_strong_bases(self, strong_bases):
        """ Set the strong bases of the stem """
        self._create_strands(length=self._length, strong_bases=strong_bases)

    def _create_strands(self, sequence: str=None, length: int=0, return_strands=False, compute_coords=True, strong_bases=True):
        ### Create the top and bottom 3D coordinates of the stem
        seq_len = len(sequence) if sequence else abs(length)

        if compute_coords:
            coords = Coords.compute_helix_from_nucl((0,0,0), # start position
                                                (1,0,0), # base vector
                                                (0,1,0), # normal vector
                                                length= seq_len,
                                                double = True)
            top_coord = Coords(coords[:seq_len])
            bot_coord = Coords(coords[seq_len:])
        else:
            top_coord = None
            bot_coord = None

        ### Create the top and bottom strands
        if sequence: # if a sequence is provided, it has the priority
            self._length = seq_len  * getattr(self, '_sign', 1)
            strands = [Strand(sequence, coords=top_coord), Strand(sequence.translate(nucl_to_pair)[::-1], directionality='53', start=(seq_len - 1, 2), direction=(-1, 0), coords=bot_coord)]
        else:
            self._length = length
            if seq_len <= 3 and strong_bases:
                seq = 'S' * seq_len
            elif self._wobble_interval:

                def get_wobble_interval():
                    if self._wobble_tolerance == 0:
                        return self._wobble_interval
                    min_wobble = max(1, 
                                     self._wobble_interval - self._wobble_tolerance)
                    return random.randint(min_wobble, 
                                          self._wobble_interval + self._wobble_tolerance)
                
                seq = ['N'] * seq_len
                # calculate the maximum index to calculate the wobble bases
                random_wobble = get_wobble_interval()
                i = 1 # the first and last nucleotides are always a normal nucleotide
                while i < seq_len - 1:
                    
                    if self._wobble_insert == "start":
                        seq[i] = 'K'
                    elif self._wobble_insert == "end" and i + random_wobble < seq_len - 1:
                        seq[i + random_wobble] = 'K'
                    elif self._wobble_insert == "middle" and i + random_wobble // 2 < seq_len - 1:
                        seq[i + random_wobble // 2] = 'K'
                    
                    # calculate the next index to insert a wobble base
                    i += random_wobble + 1
                    # calculate a new random wobble frequency
                    random_wobble = get_wobble_interval() 

                seq = ''.join(seq)

            else:
                seq = 'N' * seq_len
            strands = [Strand(seq, coords=top_coord), Strand(seq.translate(nucl_to_pair)[::-1], directionality='53', start=(seq_len - 1, 2), direction=(-1, 0), coords=bot_coord)]

        if return_strands:
            return strands
        
        self.replace_all_strands(strands, copy=False, join=False)
