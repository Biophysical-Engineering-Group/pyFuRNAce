import warnings
import random
from .symbols import *
from .callback import Callback


class Sequence(Callback):
    ### 
    ### MAGIC METHODS
    ###

    def __init__(self, sequence: str = '', directionality: str = '53', **kwargs):
        super().__init__(**kwargs)
        if directionality not in ('53', '35'):
            raise ValueError(f"Sequence directionality not allowed. It must be either '53' or '35', got {directionality} instead.")
        self._directionality = directionality
        self._check_line(sequence)
        self._sequence = str(sequence).upper().translate(only_nucl)

    def __repr__(self):
        return self._directionality[0] + ' ' + self._sequence + ' ' + self._directionality[1]
    
    def __str__(self):
        return self._sequence
        
    def __len__(self):
        return len(self._sequence)
    
    def __getitem__(self, idx):
        dir_slice = 1
        if isinstance(idx, slice):
            if idx.step is not None and idx.step < 0:
                dir_slice = -1

        return Sequence(self._sequence[idx], 
                        directionality=self.directionality[::dir_slice],
                        callbacks=self._callbacks)
    
    def __bool__(self):
        return bool(self._sequence)

    def __setitem__(self, idx, val):
        self._check_line(val)
        if isinstance(idx, slice):
            seq_line = list(self._sequence)
            seq_line[idx] = val
            self._sequence = "".join(seq_line)
        else:
            if idx < 0:
                idx = len(self) + idx
            self._sequence = str(self._sequence[: idx] + val + self._sequence[idx+1:])
        self._trigger_callbacks()

    def __add__(self, other):
        self._check_addition(other)
        return Sequence(str(self)+ str(other), self.directionality)
    
    def __mul__(self, other):
        if isinstance(other, int):
            return Sequence(self._sequence * other, self.directionality)
        raise ValueError(f'Can only multiply sequence by an integer, got {type(other)} instead')
    
    def __rmul__(self, other):
        return self.__mul__(other)
    
    def __iadd__(self, other):
        self._check_addition(other)
        self._sequence = self._sequence + str(other)
        self._trigger_callbacks()
        return self

    def __radd__(self, other):
        if other == 0:
            return self
        elif isinstance(other, str):
            return Sequence(other + str(self), self.directionality)
        else:
            return self.__add__(other)
        
    def __eq__(self, other):
        """ Check that the two sequence have same string and directionality"""
        if isinstance(other, str):
            return str(self) == other
        elif isinstance(other, Sequence):
            if self.directionality != other.directionality:
                return str(self) == str(other)[::-1]
            else:
                return str(self) == str(other)
        return False
    
    def __contains__(self, other):
        """Check if a subsequence is included in the sequence"""
        if isinstance(other, (str, Sequence)):
            str_seq = str(other)
            if isinstance(other, Sequence) and other.directionality != self.directionality:
                str_seq = str_seq[::-1]
            return str_seq in str(self)
        return False
    
    def __iter__(self):
        return iter(self._sequence)
    
    def __hash__(self):
        return hash(str(self))

    ### 
    ### PROPERTIES
    ###
 
    @property
    def directionality(self):
        """ directionality of the sequence (either '53' or '35'): the directionality is set in the initialization of the sequence. """
        return self._directionality
    
    @directionality.setter
    def directionality(self, new_directionality):
        #check that the new directionality is either '53' or '35' which are the only allowed strings
        if new_directionality not in ('53', '35'):
                raise ValueError(f"Sequence directionality not allowed. It must be either '53' or '35', got {new_directionality} instead.")
        self._directionality = new_directionality
        self._trigger_callbacks()

    ### 
    ### CHECK METHODS
    ###

    def _check_line(self, line):
        if not isinstance(line, (str, Sequence)):
            raise ValueError(f"The sequence must be a string or a sequence object. Got {type(line)} instead.")
        if isinstance(line, str) and line.translate(nucl_to_none):
            warnings.warn(f"Warning: The string '{line}' contains nucleotides not allowed in ROAD that will be removed. The allowed nucleotides are: {nucl_to_none}.", AmbiguosStructure, stacklevel=3)
        if self.directionality[0] in line[1:] or self.directionality[1] in line[:-1]:
            raise ValueError(f"The start/end symbols '{self.directionality}' are not at the end of the sequence: '{line}'") 
        return True

    def _check_addition(self, other):
        if isinstance(other, str):
            self._check_line(other)
        elif not isinstance(other, Sequence):
            raise ValueError(f'{other} is not a valid type for addition')
        elif self.directionality != other.directionality:
            raise ValueError(f'Cannot add two sequences with different directionalitys')
        return True

    ### 
    ### METHODS
    ###

    def translate(self, dictionary, inplace=False):
        """ Translate the sequence using a dictionary
        --------------------------------------------------------------------------------------
        dictionary: dict
            a dictionary that contains the translation of the sequence
        """
        if not inplace:
            return Sequence(self._sequence.translate(dictionary), directionality=self.directionality)
        new_sequence = str(self._sequence.translate(dictionary))
        self._check_line(new_sequence)
        self._sequence = new_sequence
        self._trigger_callbacks()
        return self

    def reverse(self, inplace=True):
        """ Reverse the sequence: it changes the directionality of the sequence."""
        if not inplace:
            return Sequence(self._sequence, directionality=self._directionality[::-1], callback=self.callbacks)
        self._directionality = self._directionality[::-1]
        self._trigger_callbacks()
        return self

    def complement(self):
        """ Return the complement of the sequence: the sequence is not changed."""
        return Sequence(self._sequence.translate(nucl_to_pair), self.directionality)
    
    def reverse_complement(self):
        """ Return the reverse complement of the sequence: the sequence is not changed."""
        return Sequence(self._sequence.translate(nucl_to_pair)[::-1], self.directionality)

    def gc_content(self, extended_alphabet=False):
        """Calculate the percentage of G, C, S and half K in the sequence.
        --------------------------------------------------------------------------------------
        Return: float
            The percentage of G, C, or S in the sequence.
        """
        total_count = len(self)
        if not total_count:
            return 0
        seq = self._sequence
        gc_count = seq.count('G') + seq.count('C') 
        if extended_alphabet:
            gc_count += seq.count('S') + sum(map(seq.count, ['M', 'R', 'Y', 'K'])) / 2 + sum(map(seq.count, ['D', 'H'])) / 3 + sum(map(seq.count, ['V', 'B'])) * 2 / 3 + sum(map(seq.count, ['N', 'X']))  / 4
        percentage = (gc_count / total_count) * 100
        return percentage
    
    def molecular_weight(self):
        """Calculate the molecular weight of the sequence.
        --------------------------------------------------------------------------------------
        Return: float
            The molecular weight of the sequence.
        """
        molecular_weight_table = {
                                'A': 347.2,
                                'G': 363.2,
                                'C': 323.2,
                                'U': 324.2}
        total_weight = 0
        for nucleotide in self._sequence:
            total_weight += molecular_weight_table.get(nucleotide, 0)
        return total_weight

    def copy(self, **kwargs):
        """ Return a copy of the sequence
        --------------------------------------------------------------------------------------
        Return: sequence
            a copy of the current sequence
        """        
        return Sequence(str(self), self.directionality, **kwargs)

    def pop(self, idx):
        """ Pop the element at index
        --------------------------------------------------------------------------------------
        idx: int
            the index at which adding the character
        val: str
            the characters to add
        """        
        seq_line = list(self._sequence)
        popped_val = seq_line.pop(idx)
        self._sequence = "".join(seq_line)
        self._trigger_callbacks()
        return popped_val
    
    def replace(self, old, new):
        """ Replace the old character with the new one
        --------------------------------------------------------------------------------------
        old: str
            the character to replace
        new: str
            the character to add
        """        
        self._check_line(new)
        self._sequence = str(self._sequence.replace(old, new))
        self._trigger_callbacks()
        return self
    
    def upper(self):
        """ Return the sequence in uppercase
        --------------------------------------------------------------------------------------
        Return: str
            the sequence in uppercase
        """        
        return self._sequence.upper()
    
    def split(self, sep=None):
        """ Split the sequence
        --------------------------------------------------------------------------------------
        sep: str
            the separator to split the sequence
        Return: list
            the list of the splitted sequence
        """        
        return self._sequence.split(sep)
    
    def get_random_sequence(self, target_structure=None, pair_map=None):
        if not target_structure and not pair_map:
            return ''.join([random.choice(list(iupac_code[nucleotide])) for nucleotide in self._sequence])
        elif target_structure and len(target_structure) != len(self):
            raise ValueError(f"The target dot-bracket must have the same length as the sequence. Got {len(target_structure)}, expected {len(self)}.")
        if not pair_map:
            pair_map = dot_bracket_to_pair_map(target_structure)
        # make a first random sequence
        seq = [random.choice(list(iupac_code[nucleotide])) for nucleotide in self._sequence]
        # paired the nucleotied that are paired in the target structure
        for k, v in pair_map.items():
            if v is None:
                continue
            # the paired nucleotide has the symbol for wobble pairings
            if self._sequence[v] == "K":
                if seq[k] == "G": seq[v] = "U"
                elif seq[k] == "U": seq[v] = "G"
            else: # normal pairing
                pair_sym = seq[k].translate(nucl_to_pair) 
                # take the nucleotides that are allowed by IUPAC code of the paired position and the pair symbols
                possible_paired_nucleotides = (iupac_code[self._sequence[v]] & iupac_code[pair_sym]) 
                if possible_paired_nucleotides:
                    seq[v] = random.choice(list(possible_paired_nucleotides))
        return "".join(seq)
            
    def find_repeated_subsequence(self, min_length=8):
        """ Find the repeated subsequence in the sequence
        --------------------------------------------------------------------------------------
        min_length: int
            the minimum length of the repeated subsequence
        Return: list
            the list of the repeated subsequences
        """        
        repeated_subsequences = set()
        length = len(self)
        for i in range(length):
            for j in range(i+min_length, length):
                subsequence = self[i:j]
                if 'N' in subsequence: # skip the subsequence that contains N
                    continue
                if subsequence in self[j:]:
                    repeated_subsequences.add(str(subsequence))
        return list(repeated_subsequences)
    
    def distance(self, other):
        """Calculate the distance between two sequences.
        --------------------------------------------------------------------------------------
        other: Sequence
            The other sequence to calculate the distance with.
        Return: int
            The distance between the two sequences.
        """
        if not isinstance(other, (Sequence, str)):
            raise ValueError("Invalid type for 'other'. Expected Sequence or String object.")
        if isinstance(other, str):
            other = Sequence(other, self.directionality)
        if len(self) != len(other):
            raise ValueError("Sequences must have the same length.")
        distance = 0 # Initialize the distance
        for ind, (nt1, nt2) in enumerate(zip(self, other)):
            # the symbols are not compatible
            if nt2 not in iupac_code[nt1] and nt1 not in iupac_code[nt2]:
                distance += 1
                print(ind, nt1, nt2)
        return distance
    