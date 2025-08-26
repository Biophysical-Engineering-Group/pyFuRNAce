from typing import Optional, List, Union, Literal
import random
from ..core.symbols import nucl_to_pair
from ..core.sequence import Sequence
from ..core.coordinates_3d import Coords
from ..core.strand import Strand
from ..core.motif import Motif


class Stem(Motif):
    """
    Represents a helical RNA stem motif consisting of complementary base-paired strands.

    The `Stem` class models a double-stranded RNA region with optional support for
    wobble base pairs and short-sequences set to strong bases.
    It can be created from a user-defined sequence or automatically
    generated based on a target length.

    Parameters
    ----------
    length : int, optional
        Number of base pairs in the stem. Ignored if `sequence` is provided.
    sequence : str or Sequence, optional
        Sequence for the top strand. If provided, the bottom strand is inferred
        as its reverse complement.
    wobble_interval : int, optional
        Spacing interval for inserting wobble base pairs, e.g., every N bases.
        Default is 5. Ignored if `sequence` is given.
    wobble_tolerance : int, optional
        Maximum random deviation from the defined wobble interval. Default is 2.
    wobble_insert : {"middle", "start", "end"}, optional
        Positioning strategy for wobble base insertions. Default is "middle".
    strong_bases : bool, optional
        If True, use strong (GC) base pairing for short stems (≤ 3 bp).
        Default is True.
    **kwargs : dict, optional
        Additional keyword arguments passed to the `Motif` superclass.

    Attributes
    ----------
    length : int
        Signed length of the stem (positive for forward orientation).
    wobble_interval : int
        Interval used for wobble base insertion. Zero if sequence is specified.
    wobble_tolerance : int
        Allowed range of randomness around the wobble interval.
    wobble_insert : str
        Strategy for where to place wobble base pairs.
    strands : list of Strand
        Two strands forming the helical stem, including sequence and 3D coordinates.
    """

    def __init__(
        self,
        length: int = 0,
        sequence: Union[str, Sequence] = "",
        wobble_interval: int = 5,
        wobble_tolerance: int = 2,
        wobble_insert: Literal["middle", "start", "end"] = "middle",
        strong_bases: bool = True,
        **kwargs,
    ) -> None:
        """
        Initialize a Stem motif, representing a double-stranded helical region.

        Parameters
        ----------
        length : int, default 0
            Number of base pairs in the stem. Ignored if `sequence` is provided.
        sequence : str or Sequence, default ""
            Nucleotide sequence for the top strand. If provided, wobble settings
            are ignored (default is "").
        wobble_interval : int, default 5
            Number of bases between wobble base pair insertions.
        wobble_tolerance : int, default 2
            Random variation range for wobble base pair placement
            (0 to `wobble_tolerance`).
        wobble_insert : ["middle", "start", "end"], default "middle"
            Strategy for wobble insertion: "middle", "start", or "end",
            default is "middle".
        strong_bases : bool, default True
            If True, use strong bases (G or C) for short stems shorter than 4 bases.
        **kwargs : dict, optional
            Additional arguments passed to the Motif superclass.

        Raises
        ------
        TypeError
            If parameter types are invalid.
        ValueError
            If `wobble_insert` is not one of {"middle", "start", "end"}.

        Returns
        -------
        None
        """
        ### set default values
        if wobble_insert not in ["middle", "start", "end"]:
            raise ValueError(
                f"Invalid value for wobble_insert: {wobble_insert}. "
                "The value must be 'middle', 'start' or 'end'."
            )
        if not isinstance(wobble_interval, int) or wobble_interval < 0:
            raise TypeError(
                f"The wobble frequency must be a positive integer,"
                f" got {wobble_interval}."
            )
        if not isinstance(wobble_tolerance, int) or wobble_tolerance < 0:
            raise TypeError(
                f"The wobble tolerance must be a positive integer, "
                f"got {wobble_tolerance}."
            )
        if not isinstance(length, int):
            raise TypeError(f"The length of a stem must be an integer, got {length}.")
        if not isinstance(sequence, (str, Sequence)):
            raise TypeError(
                f"The sequence of a stem must be a string or a Sequence "
                f"object, got {type(sequence)}."
            )

        self._wobble_interval = wobble_interval if not sequence else 0
        self._wobble_tolerance = wobble_tolerance if not sequence else 0
        self._wobble_insert = wobble_insert
        self._length = length
        if sequence:
            self._length = len(sequence) * getattr(self, "_sign", 1)

        ### If the user doesn't provide strands, update them directly
        if "strands" in kwargs:
            strands = kwargs.pop("strands")
        else:
            ### create the strands
            strands = self._create_strands(
                sequence=sequence,
                length=length,
                return_strands=True,
                strong_bases=strong_bases,
            )

        kwargs["join"] = False
        # Initialize the motif
        super().__init__(strands=strands, **kwargs)

    ###
    ### PROPERTIES
    ###

    @property
    def length(self):
        """Number of nucleotides in a stem"""
        return self._length

    @length.setter
    def length(self, new_length):
        if not isinstance(new_length, int):
            raise TypeError(
                f"The length of a stem must be an integer, " f"got {new_length}."
            )
        self._create_strands(length=new_length)

    @property
    def wobble_interval(self):
        return self._wobble_interval

    @wobble_interval.setter
    def wobble_interval(self, new_freq):
        if not isinstance(new_freq, int) or new_freq < 0:
            raise TypeError(
                f"The wobble frequency must be a positive integer, " f"got {new_freq}."
            )
        self._wobble_interval = new_freq
        # update the sequence of the top strand and the bottom strand
        self.length = self._length

    @property
    def wobble_tolerance(self):
        return self._wobble_tolerance

    @wobble_tolerance.setter
    def wobble_tolerance(self, new_tolerance):
        if not isinstance(new_tolerance, int) or new_tolerance < 0:
            raise TypeError(
                f"The wobble tolerance must be a positive integer, "
                f"got {new_tolerance}."
            )
        self._wobble_tolerance = new_tolerance
        # update the sequence of the top strand and the bottom strand
        self.length = self._length

    @property
    def wobble_insert(self):
        return self._wobble_insert

    @wobble_insert.setter
    def wobble_insert(self, new_insert):
        if new_insert not in ["middle", "start", "end"]:
            raise ValueError(
                f"Invalid value for wobble_insert: {new_insert}. "
                f"The value must be 'middle', 'start' or 'end'."
            )
        self._wobble_insert = new_insert
        # update the sequence of the top strand and the bottom strand
        self.length = self._length

    ###
    ### METHOD
    ###

    def set_up_sequence(self, new_seq):
        """Set the sequence of the top strand"""
        if not isinstance(new_seq, str):
            raise TypeError(f"The sequence of a stem must be a string, got {new_seq}.")
        self._create_strands(sequence=new_seq)

    def set_down_sequence(self, new_seq):
        """Set the sequence of the bottom strand"""
        self.set_up_sequence(sequence=new_seq.translate(nucl_to_pair)[::-1])

    def set_strong_bases(self, strong_bases):
        """Set wether to use strong bases for short stems"""
        self._create_strands(length=self._length, strong_bases=strong_bases)

    def _create_strands(
        self,
        sequence: Optional[str] = None,
        length: int = 0,
        return_strands: bool = False,
        compute_coords: bool = True,
        strong_bases: bool = True,
    ) -> Optional[List[Strand]]:
        """
        Internal method to create the top and bottom strands for the stem motif.

        Parameters
        ----------
        sequence : str, optional
            Nucleotide sequence for the top strand. If provided, it takes priority
            over `length`.
        length : int, default 0
            Length of the stem in nucleotides, used if `sequence` is not provided.
        return_strands : bool, default False
            If True, return the generated strands instead of assigning them.
        compute_coords : bool, default True
            Whether to compute 3D coordinates for the strands.
        strong_bases : bool, default True
            Whether to enforce strong base pairs for very short sequences.

        Returns
        -------
        list of Strand or None
            The generated strands if `return_strands` is True, otherwise None.
        """
        ### Create the top and bottom 3D coordinates of the stem
        seq_len = len(sequence) if sequence else abs(length)

        if compute_coords:
            coords = Coords.compute_helix_from_nucl(
                (0, 0, 0),  # start position
                (1, 0, 0),  # base vector
                (0, 1, 0),  # normal vector
                length=seq_len,
                double=True,
            )
            top_coord = Coords(coords[:seq_len])
            bot_coord = Coords(coords[seq_len:])
        else:
            top_coord = None
            bot_coord = None

        ### Create the top and bottom strands
        if sequence:  # if a sequence is provided, it has the priority

            if not isinstance(sequence, Sequence):
                sequence = Sequence(sequence, directionality="53")

            self._length = seq_len * getattr(self, "_sign", 1)
            strands = [
                Strand(sequence, coords=top_coord),
                Strand(
                    sequence.translate(nucl_to_pair)[::-1],
                    directionality="53",
                    start=(seq_len - 1, 2),
                    direction=(-1, 0),
                    coords=bot_coord,
                ),
            ]
        else:
            self._length = length
            if seq_len <= 3 and strong_bases:
                seq = "S" * seq_len
            elif self._wobble_interval:

                def get_wobble_interval():
                    if self._wobble_tolerance == 0:
                        return self._wobble_interval
                    min_wobble = max(1, self._wobble_interval - self._wobble_tolerance)
                    return random.randint(
                        min_wobble, self._wobble_interval + self._wobble_tolerance
                    )

                seq = ["N"] * seq_len
                # calculate the maximum index to calculate the wobble bases
                random_wobble = get_wobble_interval()
                i = 1  # the first and last nucleotides are always a normal nucleotide
                while i < seq_len - 1:

                    if self._wobble_insert == "start":
                        seq[i] = "K"

                    elif (
                        self._wobble_insert == "end" and i + random_wobble < seq_len - 1
                    ):
                        seq[i + random_wobble] = "K"

                    elif (
                        self._wobble_insert == "middle"
                        and i + random_wobble // 2 < seq_len - 1
                    ):
                        seq[i + random_wobble // 2] = "K"

                    # calculate the next index to insert a wobble base
                    i += random_wobble + 1
                    # calculate a new random wobble frequency
                    random_wobble = get_wobble_interval()

                seq = "".join(seq)

            else:
                seq = "N" * seq_len
            strands = [
                Strand(seq, coords=top_coord),
                Strand(
                    seq.translate(nucl_to_pair)[::-1],
                    directionality="53",
                    start=(seq_len - 1, 2),
                    direction=(-1, 0),
                    coords=bot_coord,
                ),
            ]

        if return_strands:
            return strands

        self.replace_all_strands(strands, copy=False, join=False)
