import os
import RNA
from . import CONFS_PATH
from ..core.symbols import *
from ..core.coordinates_3d import Coords
from ..core.strand import Strand
from .loops import Loop

### File Location for the kissing loop energy dictionaries
__location__ = os.path.realpath(
    os.path.join(os.getcwd(), os.path.dirname(__file__)))

# https://doi.org/10.2210/pdb2D1B/pdb
class KissingLoop(Loop):
    
    def __init__(self, open_left = False, sequence: str = "", seq_len=0, pk_index: str|int = "0", energy: float = -9, energy_tolerace: float = 1.0, **kwargs):

        self._pk_index = self._check_pk_index(pk_index)
        self._seq_len = seq_len
        self._energy = energy
        self._energy_tolerance = energy_tolerace
        if 'strands' in kwargs:
            strands = kwargs.pop('strands')
        else:
            strands = self._create_strands(sequence=sequence, return_strand = True, pk_index=pk_index)
        # create motif with the strands making up the external kissing loop structure
        super().__init__(strands=strands, open_left=open_left, **kwargs)
        # insert the pk_index in the strand

    ### 
    ### Properties
    ###
    
    @property
    def pk_index(self):
        """ Returns the pseudoknot symbol of the kissing loop """
        return self._pk_index
    
    @pk_index.setter
    def pk_index(self, new_index):
        self._create_strands(sequence=self.get_kissing_sequence(), pk_index = new_index)

    @property
    def energy_tolerance(self):
        """ Returns the energy tolerance of the internal kissing loop """
        return self._energy_tolerance
    
    @energy_tolerance.setter
    def energy_tolerance(self, new_energy_tolerance):
        """ Set the energy tolerance of the internal kissing loop """
        if not isinstance(new_energy_tolerance, (float, int)) or new_energy_tolerance < 0:
            raise ValueError(f"The energy tolerance should be a positive number.")
        self._energy_tolerance = new_energy_tolerance
        for strand in self:
            if hasattr(strand, 'pk_info'):
                strand.pk_info['dE'] = [new_energy_tolerance]
        self._trigger_callbacks()

    @property
    def energy(self):
        """ Returns the energy of the internal kissing loop """
        return self._energy
    
    @energy.setter
    def energy(self, new_energy):
        """ Set the energy of the internal kissing loop """
        new_energy = round(float(new_energy), 2)
        self._energy = new_energy
        for strand in self:
            if hasattr(strand, 'pk_info'):
                strand.pk_info['E'] = [new_energy]
        self._trigger_callbacks()

    ###
    ### METHODS
    ###

    def get_kissing_sequence(self):
        """ Returns the kissing sequence of the kissing loop """
        return self[0].sequence

    def set_sequence(self, new_seq):
        """ Set the sequence of the strand """
        new_seq = str(new_seq)
        self._create_strands(sequence = new_seq, pk_index = self._pk_index)

    def _check_pk_index(self, pk_index):
        if pk_index is None:
            pk_index = '0'
        elif not isinstance(pk_index, (int, str)):
            raise ValueError(f"The pk_index should be an integer or a string.")
        elif isinstance(pk_index, int):
            pk_index = str(pk_index) + "'" * (pk_index < 0)
        return pk_index

    def _create_strands(self, sequence = "", return_strand = False, pk_index = None):
        # check the pk_index
        self._pk_index = self._check_pk_index(pk_index)
    
        seq_len = self._seq_len
        if sequence:
            if seq_len and len(sequence) != seq_len:
                raise ValueError(f"The sequence length doesn't match the length required for this kissing loop, which is {seq_len}.")
            elif not seq_len: # if the sequence length is not defined, calculate the length
                seq_len = len(sequence)
            if all([s in 'ACGU' for s in sequence]):
                self._energy = round(RNA.fold(sequence + '&' + sequence.translate(nucl_to_pair)[::-1])[1], 2)
                self._energy_tolerance = 0
        else:
            sequence = 'N' * seq_len
        
        ### if the strands are already created, just update the sequence
        if hasattr(self, '_strands'):
            self._strands[0].sequence = sequence
            return self._strands
        
        ### create the strand
        strand = Strand('─' * seq_len + '╰│╭' + sequence, start=(seq_len, 2), direction=(-1, 0))
        pk_info = {"id": [self._pk_index], 'ind_fwd': [(0, seq_len - 1)], 'E': [self._energy], 'dE': [self._energy_tolerance]}
        setattr(strand, 'pk_info', pk_info)

        ### if we don't want to replace the strands, just return the strand, otherwise replace the strands
        if return_strand:
            return [strand]
        # replace the strands
        self.replace_all_strands([strand], copy=False, join=False)

    def copy(self, **kwargs):
        # return a copy of the KissingLoop or subclass
        return  type(self)(strands = self.copy_strands_preserve_blocks(self._strands), 
                           pk_index= self._pk_index, energy=self._energy, energy_tolerace=self._energy_tolerance,
                            #  Motif parameters
                            basepair = self._basepair, autopairing = self.autopairing,
                            **kwargs)

class KissingLoop120(KissingLoop):

    def __init__(self, open_left = False, sequence: str = "", pk_index: str|int = '0', energy: float = -9.0, energy_tolerace: float = 1.0, **kwargs):
        kwargs['seq_len'] = 7
        super().__init__(open_left = open_left, sequence = sequence, pk_index = pk_index, energy = energy, energy_tolerace = energy_tolerace, **kwargs)
        self[0]._coords = Coords.load_from_file(CONFS_PATH / 'KissingLoop120.dat', 
                                                dummy_ends=(True, True))

class KissingLoop180(KissingLoop):

    def __init__(self, open_left = False, sequence: str = "", pk_index: str|int = '0', energy: float = -9.0, energy_tolerace: float = 1.0, **kwargs):
        kwargs['seq_len'] = 6
        super().__init__(open_left = open_left, sequence = sequence, pk_index = pk_index, energy = energy, energy_tolerace = energy_tolerace, **kwargs)

    def get_kissing_sequence(self):
        return super().get_kissing_sequence()[2:-1]

    def _create_strands(self, sequence = "", return_strand = False, pk_index = 0):
        self._pk_index = self._check_pk_index(pk_index)
        if sequence:
            if self._seq_len and len(sequence) != self._seq_len:
                raise ValueError(f"The sequence length doesn't match the length required for this kissing loop, which is {self._seq_len}.")
            if all([s in 'ACGU' for s in sequence]):
                self._energy = round(RNA.fold('A' + sequence + 'A&A' + sequence.translate(nucl_to_pair)[::-1] + 'A')[1], 2)
                self._energy_tolerance = 0
            sequence = sequence
        else:
            sequence = 'N' * self._seq_len
        
        # if the strands are already created, just update the sequence
        if hasattr(self, '_strands'):
            self._strands[0].sequence = 'AA' + sequence + 'A'
            return self._strands
        # create the strand
        strand = Strand(f"A╭╯────A───╰╭{sequence}─╯│╭─A", start=(10, 2), direction=(-1, 0))
        ### COORDINATES FROM OXVIEW HELIX 
        strand._coords = Coords.load_from_file(CONFS_PATH / 'KissingLoop180.dat')
        pk_info = {"id": [self._pk_index], 'ind_fwd': [(2, 7)], 'E': [self._energy], 'dE': [self._energy_tolerance]}
        setattr(strand, 'pk_info', pk_info)

        # if we don't want to replace the strands, just return the strand
        if return_strand:
            return [strand]
        # replace the strands
        self.replace_all_strands([strand], copy=False, join=False)


class BranchedKissingLoop(KissingLoop):

    def __init__(self, open_left = False, sequence: str = "", pk_index: str|int = '0', energy: float = -9.0, energy_tolerace: float = 1.0, **kwargs):
        kwargs['seq_len'] = 6
        super().__init__(open_left = open_left, sequence = sequence, pk_index = pk_index, energy = energy, energy_tolerace = energy_tolerace, **kwargs)

    def get_kissing_sequence(self):
        """ Returns the kissing sequence of the kissing loop """
        strand_ind = [i for i, s in enumerate(self) if len(s.sequence) == 7][0]
        return self[strand_ind].sequence[:-1]

    def _create_strands(self, sequence = "", return_strand = False, pk_index = 0):
        self._pk_index = self._check_pk_index(pk_index)
        if sequence:
            if self._seq_len and len(sequence) != self._seq_len:
                raise ValueError(f"The sequence length doesn't match the length required for this kissing loop, which is {self._seq_len}.")
            if all([s in 'ACGU' for s in sequence]):
                self._energy = round(RNA.fold(sequence + '&' + sequence.translate(nucl_to_pair)[::-1])[1], 2)
                self._energy_tolerance = 0
            sequence = sequence
        else:
            sequence = 'N' * self._seq_len
        
        # if the strands are already created, just update the sequence
        if hasattr(self, '_strands'):
            strand_ind = [i for i, s in enumerate(self) if len(s.sequence) == 7][0]
            self._strands[strand_ind].sequence = sequence + 'A'
            return self._strands
        # create the strand
        kissing_strand = Strand(f'╮──╰╭{sequence}─╯│╭A', start=(3, 3), direction=(0, -1))
        kissing_strand._coords = Coords.load_from_file(CONFS_PATH / 'BranchedKissingLoop_1.dat',
                                                       dummy_ends=(True, False))
        pk_info = {"id": [self._pk_index], 'ind_fwd': [(0, 5)], 'E': [self._energy], 'dE': [self._energy_tolerance]}
        setattr(kissing_strand, 'pk_info', pk_info)

        connect_strand = Strand('╭╯────╭', start=(9, 2), direction=(-1, 0))
        connect_strand._coords = Coords.load_from_file(CONFS_PATH / 'BranchedKissingLoop_2.dat',
                                                       dummy_ends=(True, True))
        strands = [kissing_strand, connect_strand]
        # if we don't want to replace the strands, just return the strand
        if return_strand:
            return strands
        # replace the strands
        self.replace_all_strands(strands, copy=False, join=False)

class KissingDimer(KissingLoop180):
    
    def __init__(self, sequence: str = "", pk_index: str|int = '0', energy: float = -9.0, energy_tolerace: float = 1.0, **kwargs):
        """
        Attributes of the class KissingDimer, which is a daugther class of the class Motif.
        -----------------------------------------------------------------------------------
        sequence: str
            nucelotide sequnce in the internal KL
        """
        super().__init__(sequence = sequence, pk_index=pk_index, energy = energy, energy_tolerace = energy_tolerace, **kwargs)

    ### 
    ### METHODS
    ###

    def _create_strands(self, sequence="", return_strand=False, pk_index = 0):
        new_pk_index = self._check_pk_index(pk_index)
        # the bottom pk_index is the inverse of the top one
        if "'" == new_pk_index[-1]:
            bottom_pk_index = new_pk_index[:-1]
        else:
            bottom_pk_index = new_pk_index + "'"

        bottom_strand = super()._create_strands(sequence, return_strand=True, pk_index=bottom_pk_index)[0]
        seq = bottom_strand.sequence[2:-1]
        rev_comp = seq.translate(nucl_to_pair)[::-1]
    
        self._pk_index = new_pk_index # add the pk_index to override the pk_index of the bottom strand

        ### if the strands are already created, just update the sequence and return the strands
        if hasattr(self, '_strands'):
            self._strands[1].sequence = 'AA' + rev_comp + 'A'
            return self._strands
        
        ### shift the second strand to make space for the second one
        bottom_strand.start = (13, 3)
        bottom_strand.sequence = 'AA' + rev_comp + 'A'  # the second strand is the reverse complement of the first

        ### create the second
        top_strand = Strand(f"A╯╭────A───╮╯{seq}─╭│╯─A", directionality='53', start=(0, 1), direction=(1, 0))
        ## COORDINATES FROM OXVIEW HELIX
        top_strand._coords = Coords.load_from_file(CONFS_PATH / 'KissingLoop180_2.dat')
        pk_info = {"id": [self._pk_index], 'ind_fwd': [(2, 7)], 'E': [self._energy], 'dE': [self._energy_tolerance]}
        setattr(top_strand, 'pk_info', pk_info)

        strands = [top_strand, bottom_strand]

        ### if we don't want to replace the strands, just return the strand, otherwise replace the strands
        if return_strand:
            return strands
        # replace the strands
        self.replace_all_strands(strands, copy=False, join=False)

    def set_top_sequence(self, new_seq):
        """ Set the sequence of the top strand """
        self.set_sequence(new_seq)

    def set_bot_sequence(self, new_seq):
        """ Set the sequence of the top strand"""
        self.set_sequence(new_seq.translate(nucl_to_pair)[::-1])

class BranchedDimer(BranchedKissingLoop):
    # strand 0: branched KL
    # strand 1: bkl connection
    # strand 2: Dimer strand
    
    def __init__(self, sequence: str = "", pk_index: str|int = '0', energy: float = -9.0, energy_tolerace: float = 1.0, **kwargs):
        """
        Attributes of the class KissingDimer, which is a daugther class of the class Motif.
        -----------------------------------------------------------------------------------
        sequence: str
            nucelotide sequnce in the internal KL
        """
        super().__init__(sequence = sequence, pk_index=pk_index, energy = energy, energy_tolerace = energy_tolerace, **kwargs)

    ### 
    ### METHODS
    ###

    def _create_strands(self, sequence="", return_strand=False, pk_index=0):
        new_pk_index = self._check_pk_index(pk_index)
        # the bottom pk_index is the inverse of the top one
        if "'" == new_pk_index[-1]:
            bottom_pk_index = new_pk_index[:-1]
        else:
            bottom_pk_index = new_pk_index + "'"

        if not sequence:
            sequence = sequence = 'N' * self._seq_len
        rev_comp = sequence.translate(nucl_to_pair)[::-1]
        strands = super()._create_strands(rev_comp, return_strand=True, pk_index=bottom_pk_index)
        
        self._pk_index = new_pk_index # add the pk_index to override the pk_index of the bottom strand

        ### if the strands are already created, just update the sequence and return the strands
        if hasattr(self, '_strands'):
            kl_index = [i for i, s in enumerate(self) if len(s.sequence) == 9][0]
            self._strands[kl_index].sequence = 'AA' + sequence + 'A'
            return self._strands
        
        ### shift the bottom branched KL
        strands[0].start = (strands[0].start[0] + 3, strands[0].start[1] + 1)
        strands[1].start = (strands[1].start[0] + 3, strands[1].start[1] + 1)

        ### create the top strand
        top_strand = Strand(f"A╯╭────A───╮╯{sequence}─╭│╯─A", directionality='53', start=(0, 1), direction=(1, 0))
        top_strand._coords = Coords.load_from_file(CONFS_PATH / 'BranchedKissingLoop_3.dat')
        pk_info = {"id": [self._pk_index], 'ind_fwd': [(2, 7)], 'E': [self._energy], 'dE': [self._energy_tolerance]}
        setattr(top_strand, 'pk_info', pk_info)

        strands.insert(0, top_strand)
        ### if we don't want to replace the strands, just return the strand, otherwise replace the strands
        if return_strand:
            return strands
        # replace the strands
        self.replace_all_strands(strands, copy=False, join=False)

    def set_top_sequence(self, new_seq):
        """ Set the sequence of the top strand """
        self.set_sequence(new_seq)

    def set_bot_sequence(self, new_seq):
        """ Set the sequence of the top strand"""
        self.set_sequence(new_seq.translate(nucl_to_pair)[::-1])