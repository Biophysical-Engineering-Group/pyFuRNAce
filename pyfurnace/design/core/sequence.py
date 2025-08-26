import warnings
import random
from typing import Literal, Union, List, Optional, Set
from .symbols import (
    iupac_code,
    only_nucl,
    nucl_to_none,
    nucleotides,
    nucl_to_pair,
    dot_bracket_to_pair_map,
    AmbiguosStructure,
)
from .callback import Callback


class Sequence(Callback):
    """
    Represents a nucleotide sequence with 5' to 3' or 3' to 5' directionality.

    This class provides functionality for manipulating RNA sequences
    with support for slicing, translation, complementarity, GC content,
    and integration with structure-aware tools.

    Parameters
    ----------
    sequence : str, optional
        The nucleotide sequence.
    directionality : {'53', '35'}, optional
        Directionality of the sequence (default is '53').
    **kwargs : dict
        Additional arguments passed to the Callback base class.

    Attributes
    ----------
    directionality : str
        Direction of the sequence ('53' or '35').
    """

    def __init__(
        self, sequence: str = "", directionality: Literal["53", "35"] = "53", **kwargs
    ) -> None:
        """
        Initializes a Sequence object.

        Parameters
        ----------
        sequence : str, optional
            The nucleotide sequence.
        directionality : str, optional
            The directionality of the sequence ('53' or '35') (default is '53').
        **kwargs : Any
            Arbitrary keyword arguments to pass to the Callback class.

        Raises
        ------
        ValueError
            If the directionality is not '53' or '35'.
        """
        super().__init__(**kwargs)

        if directionality not in ("53", "35"):
            raise ValueError(
                f"Sequence directionality not allowed. "
                f"It must be either '53' or '35', "
                f"got {directionality} instead."
            )

        self._directionality = directionality
        sequence = self._check_line(sequence)
        self._sequence = str(sequence).translate(only_nucl)

    def __str__(self) -> str:
        """Return the string representation of the sequence."""
        return self._sequence

    def __repr__(self) -> str:
        """Return a string representation of the sequence object."""
        return f"{self._directionality[0]} {self._sequence} {self._directionality[1]}"

    def __getitem__(self, idx: Union[int, slice]) -> "Sequence":
        """Return a subsequence of the sequence."""
        dir_slice = 1
        if isinstance(idx, slice):
            if idx.step is not None and idx.step < 0:
                dir_slice = -1

        return Sequence(
            self._sequence[idx],
            directionality=self.directionality[::dir_slice],
            callbacks=self._callbacks,
        )

    def __setitem__(self, idx: Union[int, slice], val: Union[str, "Sequence"]) -> None:
        """Set a subsequence of the sequence."""
        val = self._check_line(val)
        if isinstance(idx, slice):
            seq_line = list(self._sequence)
            seq_line[idx] = val.upper()
            self._sequence = "".join(seq_line)
        else:
            if idx < 0:
                idx = len(self) + idx
            self._sequence = str(self._sequence[:idx] + val + self._sequence[idx + 1 :])
        self._trigger_callbacks(new_sequence=self._sequence)

    def __len__(self) -> int:
        """Return the length of the sequence."""
        return len(self._sequence)

    def __add__(self, other: Union[str, "Sequence"]) -> "Sequence":
        """Concatenate two sequences or a sequence and a string"""
        self._check_addition(other)
        return Sequence(str(self) + str(other), self.directionality)

    def __iadd__(self, other: Union[str, "Sequence"]) -> "Sequence":
        """Concatenate a sequence or a string to the current sequence"""
        self._check_addition(other)
        self._sequence = self._sequence + str(other)
        self._trigger_callbacks(new_sequence=self._sequence)
        return self

    def __radd__(self, other: Union[str, "Sequence"]) -> "Sequence":
        """Concatenate two sequences or a sequence and a string"""
        if other == 0:
            return self
        return Sequence(str(other) + str(self), self.directionality)

    def __mul__(self, other: int) -> "Sequence":
        """Multiply the sequence by an integer."""
        if isinstance(other, int):
            return Sequence(self._sequence * other, self.directionality)
        raise ValueError(
            f"Can only multiply sequence by an integer, " f"got {type(other)} instead"
        )

    def __rmul__(self, other: int) -> "Sequence":
        """Multiply the sequence by an integer."""
        return self.__mul__(other)

    def __bool__(self) -> bool:
        """Check if the sequence is not empty."""
        return bool(self._sequence)

    def __contains__(self, other: Union[str, "Sequence"]) -> bool:
        """Check if a subsequence is included in the sequence"""
        if isinstance(other, (str, Sequence)):
            str_seq = str(other)
            if (
                isinstance(other, Sequence)
                and other.directionality != self.directionality
            ):
                str_seq = str_seq[::-1]
            return str_seq in str(self)
        return False

    def __iter__(self):
        """Return an iterator over the sequence."""
        return iter(self._sequence)

    def __eq__(self, other: Union[str, "Sequence"]) -> bool:
        """Check that the two sequence have same string and directionality"""
        if isinstance(other, str):
            return str(self) == other
        elif isinstance(other, Sequence):
            if self.directionality != other.directionality:
                return str(self) == str(other)[::-1]
            else:
                return str(self) == str(other)
        return False

    def __hash__(self):
        """Return a hash of the sequence."""
        return hash(self.__repr__())

    ###
    ### PROPERTIES
    ###

    @property
    def directionality(self):
        """directionality of the sequence (either '53' or '35')"""
        return self._directionality

    @directionality.setter
    def directionality(self, new_directionality: Literal["53", "35"]) -> None:
        """Set the directionality of the sequence."""
        if new_directionality not in ("53", "35"):
            raise ValueError(
                f"Sequence directionality not allowed. "
                f"It must be either '53' or '35', got "
                f"{new_directionality} instead."
            )
        self._directionality = new_directionality
        self._trigger_callbacks(new_sequence=self._sequence)

    ###
    ### PRETECTED METHODS
    ###

    def _check_addition(self, other: Union[str, "Sequence"]) -> bool:
        """Check that the sequence is valid and have same directionality"""
        if isinstance(other, str):
            other = self._check_line(other)

        elif not isinstance(other, Sequence):
            raise ValueError(f"{other} is not a valid type for addition")

        elif self.directionality != other.directionality:
            raise ValueError("Cannot add two sequences with different directionality")
        return True

    def _check_line(self, line: Union[str, "Sequence"]) -> bool:
        """Check that the sequence is valid and contains only allowed nucleotides."""
        if not isinstance(line, (str, Sequence)):
            raise ValueError(
                f"The sequence must be a string or a sequence object. "
                f"Got {type(line)} instead."
            )

        if isinstance(line, str):
            line = line.upper()
            if line.translate(nucl_to_none).replace("&", ""):
                warnings.warn(
                    f"Warning: The string '{line}' contains nucleotides not"
                    f" allowed in ROAD that will be removed. The allowed "
                    f"nucleotides are: {nucleotides.union('&')}.",
                    AmbiguosStructure,
                    stacklevel=3,
                )

        if self.directionality[0] in line[1:] or self.directionality[1] in line[:-1]:
            raise ValueError(
                f"The start/end symbols '{self.directionality}' "
                f"are not at the end of the sequence: '{line}'"
            )

        return line

    ###
    ### METHODS
    ###

    def complement(self) -> "Sequence":
        """
        Return the complement of the sequence.

        Returns
        -------
        Sequence
            Complementary sequence.
        """
        return Sequence(self._sequence.translate(nucl_to_pair), self.directionality)

    def copy(self, **kwargs) -> "Sequence":
        """
        Create a copy of the sequence.

        Returns
        -------
        Sequence
            A new Sequence instance.
        """
        return Sequence(str(self), self.directionality, **kwargs)

    def distance(self, other: Union[str, "Sequence"]) -> int:
        """
        Compute the number of mismatched bases between this sequence and another.

        Parameters
        ----------
        other : str or Sequence
            Sequence to compare against.

        Returns
        -------
        int
            Number of incompatible positions.

        Raises
        ------
        ValueError
            If the input is invalid or lengths do not match.
        """
        other = self._check_line(other)

        if len(self) != len(other):
            raise ValueError("Sequences must have the same length.")

        distance = 0  # Initialize the distance
        for ind, (nt1, nt2) in enumerate(zip(self, other)):
            # the symbols are not compatible
            if nt2 not in iupac_code[nt1] and nt1 not in iupac_code[nt2]:
                distance += 1

        return distance

    def find(self, sub: str) -> int:
        """
        Find the first occurrence of a subsequence in the sequence.

        Parameters
        ----------
        sub : str
            Subsequence to find.

        Returns
        -------
        int
            Index of the first occurrence of the subsequence, or -1 if not found.
        """
        return self._sequence.find(sub)

    def find_repeated_subsequence(self, min_length: int = 8) -> Set[str]:
        """
        Find all subsequences of minimum length that appear more than once.

        Parameters
        ----------
        min_length : int, optional
            Minimum length of subsequence (default is 8).

        Returns
        -------
        set of str
            Repeated subsequences.
        """
        repeated_subsequences = set()
        length = len(self)
        for i in range(length):
            for j in range(i + min_length, length):

                subsequence = self[i:j]

                if subsequence in self[j:]:
                    repeated_subsequences.add(str(subsequence))

        return repeated_subsequences

    def gc_content(self, extended_alphabet: bool = False) -> float:
        """
        Calculate the GC content of the sequence.

        Parameters
        ----------
        extended_alphabet : bool, optional
            Whether to include ambiguous base codes (e.g., S, M, R) in the calculation.

        Returns
        -------
        float
            GC content percentage.
        """
        total_count = len(self)
        if not total_count:
            return 0
        seq = self._sequence
        gc_count = seq.count("G") + seq.count("C")
        if extended_alphabet:
            gc_count += (
                seq.count("S")
                + sum(map(seq.count, ["M", "R", "Y", "K"])) / 2
                + sum(map(seq.count, ["D", "H"])) / 3
                + sum(map(seq.count, ["V", "B"])) * 2 / 3
                + sum(map(seq.count, ["N", "X"])) / 4
            )
        percentage = (gc_count / total_count) * 100
        return percentage

    def get_random_sequence(self, structure: str = "") -> "Sequence":
        """
        Generate a randomized sequence compatible with the IUPAC symbols and
        optionally a dot-bracket structure to be respected.

        Parameters
        ----------
        structure : str, optional
            Dot-bracket structure.

        Returns
        -------
        str
            A new randomized sequence.
        """
        # make a first random sequence
        seq = [
            random.choice(list(iupac_code[nucleotide])) for nucleotide in self._sequence
        ]

        # fully random sequence
        if not structure:
            return Sequence("".join(seq), self.directionality)

        # small sanity check of the structure
        elif structure and len(structure) != len(self):
            raise ValueError(
                f"The target dot-bracket must have the same length as "
                f"the sequence. Got {len(structure)}, "
                f"expected {len(self)}."
            )

        # build the pair map
        pair_map = dot_bracket_to_pair_map(structure)

        # paired the nucleotied that are paired in the target structure
        for k, v in pair_map.items():

            # unpaired nucleotide
            if v is None:
                continue

            # the paired nucleotide has the symbol for wobble pairings
            if self._sequence[v] == "K":
                if seq[k] == "G":
                    seq[v] = "U"
                elif seq[k] == "U":
                    seq[v] = "G"

            # normal pairing
            else:
                # the iupac code that pairs the nucleotides
                pair_sym = seq[k].translate(nucl_to_pair)
                # the iupac code of the nucleotide at position v
                sym_at_pos = iupac_code[self._sequence[v]]
                # take the nucleotides that are allowed by both
                possible_paired_nucleotides = sym_at_pos & iupac_code[pair_sym]
                if possible_paired_nucleotides:
                    seq[v] = random.choice(list(possible_paired_nucleotides))

        return Sequence("".join(seq), self.directionality)

    def lower(self) -> str:
        """
        Convert sequence to lowercase.

        Returns
        -------
        str
            Lowercase sequence string.
        """
        return self._sequence.lower()

    def molecular_weight(self) -> float:
        """
        Calculate the molecular weight of the sequence.

        Returns
        -------
        float
            Total molecular weight in Daltons.
        """
        molecular_weight_table = {"A": 347.2, "G": 363.2, "C": 323.2, "U": 324.2}
        total_weight = 0
        for nucleotide in self._sequence:
            total_weight += molecular_weight_table.get(nucleotide, 0)
        return total_weight

    def pop(self, idx: int) -> str:
        """
        Remove and return a nucleotide at the given index.

        Parameters
        ----------
        idx : int
            Index of nucleotide to remove.

        Returns
        -------
        str
            The removed character.
        """
        seq_line = list(self._sequence)
        popped_val = seq_line.pop(idx)
        self._sequence = "".join(seq_line)
        self._trigger_callbacks(new_sequence=self._sequence)
        return popped_val

    def replace(self, old: str, new: str) -> "Sequence":
        """
        Replace all occurrences of `old` with `new` in the sequence.

        Parameters
        ----------
        old : str
            Character to replace.
        new : str
            Replacement character.

        Returns
        -------
        Sequence
            Updated sequence.
        """
        new = self._check_line(new)
        self._sequence = str(self._sequence.replace(old, new))
        self._trigger_callbacks(new_sequence=self._sequence)
        return self

    def reverse(self, inplace: bool = False) -> "Sequence":
        """
        Reverse the directionality of the sequence.

        Parameters
        ----------
        inplace : bool, optional
            If True, reverse in-place. Otherwise, return a new Sequence.

        Returns
        -------
        Sequence
            Reversed sequence.
        """
        if not inplace:
            return Sequence(
                self._sequence,
                directionality=self._directionality[::-1],
                callback=self.callbacks,
            )

        self._directionality = self._directionality[::-1]
        self._trigger_callbacks(new_sequence=self._sequence)
        return self

    def reverse_complement(self) -> "Sequence":
        """
        Return the reverse complement of the sequence.

        Returns
        -------
        Sequence
            Reverse-complemented sequence.
        """
        return Sequence(
            self._sequence.translate(nucl_to_pair)[::-1], self.directionality
        )

    def split(self, sep: Optional[str] = None) -> List["Sequence"]:
        """
        Split the sequence by a separator.

        Parameters
        ----------
        sep : str, optional
            Separator to use (default is None, meaning split on whitespace).

        Returns
        -------
        list of Sequence
            List of subsequences.
        """
        return [Sequence(s, self.directionality) for s in self._sequence.split(sep)]

    def translate(self, dictionary: dict, inplace: bool = False) -> "Sequence":
        """
        Translate the sequence using a mapping dictionary.

        Parameters
        ----------
        dictionary : dict
            Dictionary to translate each character.
        inplace : bool, optional
            If True, modify the current sequence. Otherwise, return a new instance.

        Returns
        -------
        Sequence
            Translated sequence.
        """
        dictionary = str.maketrans(dictionary)

        if not inplace:
            return Sequence(
                self._sequence.translate(dictionary), directionality=self.directionality
            )

        new_sequence = str(self._sequence.translate(dictionary))
        new_sequence = self._check_line(new_sequence)
        self._sequence = new_sequence
        self._trigger_callbacks(new_sequence=self._sequence)
        return self

    def upper(self) -> str:
        """
        Convert sequence to uppercase.

        Returns
        -------
        str
            Uppercase sequence string.
        """
        return self._sequence.upper()
