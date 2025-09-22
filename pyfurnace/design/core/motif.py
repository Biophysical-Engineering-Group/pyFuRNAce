import warnings
import copy
from functools import wraps
from inspect import signature
from collections.abc import Iterable
from typing import (
    Any,
    List,
    Dict,
    Tuple,
    Union,
    Optional,
    Literal,
    Callable,
    TYPE_CHECKING,
)

# OAT IMPORTS
try:
    from oxDNA_analysis_tools.external_force_utils.force_reader import write_force_file
    from oxDNA_analysis_tools.external_force_utils.forces import mutual_trap
    from oxDNA_analysis_tools.UTILS.RyeReader import (
        get_confs,
        describe,
        strand_describe,
        inbox,
    )
    from oxDNA_analysis_tools.oxDNA_PDB import oxDNA_PDB

    oat_installed = True
except ImportError:
    oat_installed = False

# pyFuRNAce IMPORTS
from .symbols import (
    pair_map_to_dot_bracket,
    dot_bracket_to_pair_map,
    Node,
    tree_to_dot_bracket,
    dot_bracket_to_tree,
    accept_symbol,
    base_pairing,
    bp_symbols,
    rotate_90,
    MotifStructureError,
    folding_barriers as fold_bar,
)
from .callback import Callback
from .sequence import Sequence
from .strand import Strand, StrandsBlock
from .basepair import BasePair
from .position import Position, Direction

if TYPE_CHECKING:  # only for type checkers / linters, not at runtime
    from .origami import Origami


class Motif(Callback):
    """
    Represents a structural motif in an RNA Origami design.

    A `Motif` is composed of multiple `Strand` objects, defining an RNA secondary
    structure with base pair interactions. It supports coordinate locking,
    strand joining, sequence assignment, and 2D motif visualization.

    Parameters
    ----------
    *strands : list of Strand, optional
        A list of `Strand` objects to initialize the motif.
    basepair : BasePair, dict, optional
        A dictionary defining base pair connections, with `Position` objects as keys
        and values representing paired positions.
    structure : str, optional
        A dot-bracket notation representation of the motif's structure.
        If this is provided, the `basepair` parameter is ignored.
    autopairing : bool, default True
        If True, automatically pairs bases based on complementary interactions.
        If a basepair or structure is given, autopairing is turned off.
    lock_coords : bool, default True
        If True, locks the coordinate system for all strands in the motif.
        When a strand coordinate transformed, all other strands are transformed.
    join : bool, default True
        If True, attempts to merge consecutive strands where possible.
    copy : bool, default False
        If True, creates copies of the input strands instead of using references.
    **kwargs : dict
        Additional arguments to pass to the parent `Callback` class.

    Attributes
    ----------
    autopairing : bool
        Indicates whether the motif automatically pairs bases.
    basepair : BasePair
        The base-pairing map of the motif.
        A double dictionary with `Position` objects as keys and values
        representing paired positions in the 2D representation.
    junctions : dict
        A dictionary with 2D Direction as key and Positions of the junctions as values.
    lock_coords : bool
        Indicates whether the strands coordinates are part of the same coordinate
        system.
    max_pos : tuple of int
        The maximum x, y coordinates occupied by the motif.
    min_pos : tuple of int
        The minimum x, y coordinates occupied by the motif.
    num_char : int
        The maximum length of the motif lines.
    num_lines : int
        The number of lines in the motif structure.
    pair_map : dict
        A mapping of paired bases indices (alternative to dot-bracket notation).
    positions : tuple of Position
        A tuple of the positions of each character of each strand in 2D space.
    seq_positions : tuple of Position
        The positions of each nucleotide in the motif sequence (x,y coordinates).
        Note: the sequence is always 5' to 3' in the motif.
    sequence : Sequence
        The nucleotide sequence of the motif. It's the concatenation of the sequence.
        of each strand, with 5' to 3' directionality.
    strands : list of Strand
        A list of strands composing the motif.
    structure : str
        The motif structure represented in dot-bracket notation.

    Notes
    -----
    - The motif is built from strands and supports strand merging.
    - Base pairs can be set manually or inferred automatically.
    - Supports transformations, file export, and visual representation.
    """

    def __init__(
        self,
        *strands: List[Strand],
        basepair: Optional[Dict[Tuple[int, int], Tuple[int, int]]] = None,
        structure: Optional[str] = None,
        autopairing: bool = True,
        lock_coords: bool = True,
        join: bool = True,
        copy: bool = False,
        **kwargs,
    ) -> None:
        """
        Initialize a Motif instance.

        Parameters
        ----------
        *strands : list of Strand
            List of strands to add to the motif.
        basepair : dict, optional
            The base-pairing map of the motif.
            A double dictionary with `Position` objects as keys and values
            representing paired positions in the 2D representation.
        structure : str, optional
            The structure of the motif in dot-bracket notation.
            If provided, the `basepair` parameter is ignored.
        autopairing : bool, default True
            Whether to enable automatic base pairing.
            If a basepair or structure is given, autopairing is turned off.
        lock_coords : bool, default True
            Whether to lock coordinates in the same coordinate system.
        join : bool, default True
            Whether to join consecutive strands.
        copy : bool, default False
            Whether to copy strands.
        """
        super().__init__(**kwargs)

        ### Initialize the attributes
        self._positions = None
        self._seq_positions = None
        self._max_pos = None
        self._min_pos = None
        self._junctions = {}
        self._sequence = Sequence()
        self._pair_map = BasePair()
        self._structure = None
        self._strands = []
        self.lock_coords = lock_coords
        self._strands_block = StrandsBlock()

        ### STRANDS INITIALIZATION
        if (
            len(strands) == 1
            and isinstance(strands[0], Iterable)
            and not isinstance(strands[0], Strand)
        ):
            # if a single list of strands is passed, unpack it
            strands = tuple(strands[0])

        extra = kwargs.pop("strands", tuple())

        if isinstance(extra, Strand):
            strands += (extra,)
        elif isinstance(extra, Iterable):
            strands += tuple(extra)
        else:
            raise TypeError(
                "`strands` keyword argument must be a Strand or an "
                "iterable of Strands."
            )

        if not all(isinstance(s, Strand) for s in strands):
            raise ValueError(
                "All the elements in the strands input must be of type"
                f" Strand. Got types {[type(s) for s in strands]}."
            )

        ### COPY/REGISTER STRANDS
        if copy:
            strands = self.copy_strands_preserve_blocks(strands, self)
        else:
            for s in strands:
                s.register_callback(self._updated_strands)

        ### JOIN STRANDS
        if join:
            strands = self.join_strands(strands=strands)
        self._strands = list(strands)  # set the strands

        ### LOCK COORDINATES
        if self.lock_coords:
            # IMPORTANT: look at the _coords to avoid locking strands
            #            that don't have coordinates
            # (therefore should not be locked)
            self._strands_block = StrandsBlock(*[s for s in self._strands])
            #  if not s._coords.is_empty()])

        ### BASEPAIR/DOT-BRACKET INITIALIZATION
        # prioritize the structure
        if structure:
            self.structure = structure
            autopairing = False

        else:
            # if basepair is given, set off autopairing
            if basepair:
                autopairing = False
            else:
                autopairing = autopairing
                basepair = BasePair()
            # initialize basepair property
            self.basepair = BasePair(basepair, callback=self._updated_basepair)

        # set autopairing
        self._autopairing = autopairing

    def __str__(self) -> str:
        """Return the blueprint of the motif."""
        if not self:
            return ""
        # create the canvas
        canvas_repr = [" " * self.num_char] * self.num_lines

        # draw each strand to the canvas
        for strand in self:
            strand.draw(canvas_repr, return_string=False)

        # add base pairing at position (if the base pairing position is free)
        for pos1, pos2 in self.basepair.items():

            # select the horizontal basepair symbol
            if pos1[1] - pos2[1] == 0 and abs(pos1[0] - pos2[0]) == 2:
                pos = ((pos1[0] + pos2[0]) // 2, pos1[1])
                sym = "="
            # select the vertical basepair symbol
            elif pos1[0] - pos2[0] == 0 and abs(pos1[1] - pos2[1]) == 2:
                pos = (pos1[0], (pos1[1] + pos2[1]) // 2)
                sym = "┊"
            # skip the basepair if it is not horizontal or vertical
            else:
                continue

            # if the position is free, add the basepair symbol
            if (
                pos[1] < len(canvas_repr)
                and pos[0] < len(canvas_repr[0])
                and canvas_repr[pos[1]][pos[0]] == " "
            ):
                canvas_repr[pos[1]] = (
                    canvas_repr[pos[1]][: pos[0]]
                    + sym
                    + canvas_repr[pos[1]][pos[0] + 1 :]
                )

            # send a warning if the position is already occupied
            else:
                warnings.warn(
                    f"Hidden basepair at position {pos}."
                    f" Trying to pair {pos1} with {pos2}.",
                    stacklevel=2,
                )

        return "\n".join(canvas_repr)

    def __repr__(self) -> str:
        """Return the list of strands."""
        return str(list(self))

    def __getitem__(self, idx: int) -> "Strand":
        """Get the strand at index."""
        return self._strands[idx]

    def __setitem__(self, idx: int, strand: "Strand") -> None:
        """Set the strand at index and try to join."""
        if not isinstance(strand, Strand):
            raise ValueError(f"{strand} is not a Strand object.")
        # register the callback to update the motif
        strand.register_callback(self._updated_strands)
        self._strands[idx] = strand
        # indicate that the strands have been updated
        self._strands = self.join_strands(self._strands)
        if self.lock_coords:
            self._strands_block = StrandsBlock(*[s for s in self._strands])
        self._updated_strands()

    def __len__(self) -> int:
        """Get the number of strands."""
        return len(self._strands)

    def __add__(self, other: "Motif") -> "Motif":
        """
        Add two motifs together by stacking them horizontally.

        Parameters
        ----------
        other : Motif
            The motif to be added.

        Returns
        -------
        Motif
            The new motif resulting from the addition.

        Raises
        ------
        ValueError
            If `other` is not a `Motif` instance.
        """
        self._check_addition(other)

        # edge cases
        if not self and not other:
            return Motif()
        elif not self:
            return other.copy()
        elif not other:
            return self.copy()

        # create a copy of the motif to add.
        # All operations will be performed on the copies
        self_copy = self.copy()
        other_copy = other.copy()

        # align the motifs horizontally
        self.align(self_copy, other_copy)
        # get the x shifts to move the motifs
        x_shifts = self.get_sequential_shift(
            [self_copy, other_copy], position_based=False
        )

        # move motif to the right of the first motif
        other_copy.shift((x_shifts[1], 0))

        # if one of the two motifs doesn't use autopairing,
        # set it false and update the basepair dictionary
        if self_copy.autopairing and other_copy.autopairing:
            new_basepair = BasePair()
        else:
            new_basepair = self_copy.basepair
            new_basepair.update(other_copy.basepair)

        # return the new motif collecting the strands
        return Motif(strands=list(self_copy) + list(other_copy), basepair=new_basepair)

    def __iadd__(self, other: "Motif") -> "Motif":
        """
        Add another motif to the current motif in place and join the strands.

        Parameters
        ----------
        other : Motif
            The motif to be added.

        Returns
        -------
        Motif
            The updated motif.
        """
        self._check_addition(other)

        # edge cases
        if not other:
            return self

        # copy the other motif to prevent modifying it
        other_copy = other.copy()

        if self:
            # align the motifs horizontally
            self.align(self, other_copy)
            x_shifts = self.get_sequential_shift(
                [self, other_copy], position_based=False
            )
            # move motif to the right of the first motif
            other_copy.shift((x_shifts[1], 0))

        # if other has autopairing off, set it off
        if not other_copy.autopairing:
            self._autopairing = False
            self._basepair.update(other_copy.basepair)

        # add the other strands to self
        for strand in other_copy:
            self._strands.append(strand)

        # join the strands
        self._strands = self.join_strands(self._strands)
        if self.lock_coords:
            self._strands_block = StrandsBlock(*[s for s in self._strands])
        for s in self:
            s.register_callback(self._updated_strands)

        self._updated_strands()
        return self

    def __radd__(self, other: "Motif") -> "Motif":
        """
        Enable right addition for motifs.
        """
        if other == 0:
            return self
        else:
            return self.__add__(other)

    def __bool__(self) -> bool:
        """True if the motif contains at least one valid strand, False otherwise."""
        if not self._strands:
            return False
        for s in self._strands:
            if s:
                return True
        return False

    def __contains__(self, other: Union["Strand", "Sequence", str, "BasePair"]) -> bool:
        """Check if the motif contains a strand, sequence, or base pair."""
        if isinstance(other, Strand):
            return other in self._strands
        elif isinstance(other, Sequence):
            return other in self.sequence
        elif isinstance(other, str):
            return any(other in s for s in self._strands)
        elif isinstance(other, BasePair):
            return all(self.basepair.get(k) == v for k, v in other.items())
        return False

    def __eq__(self, other: Any) -> bool:
        """True if two motifs have equal strands and equal basepair."""
        if not isinstance(other, Motif):
            return False

        elif len(self) != len(other):
            return False

        for s1, s2 in zip(self, other):
            if s1 != s2:
                return False

        if self.basepair != other.basepair:
            return False

        return True

    ###
    ### PROPERTIES
    ###

    @property
    def autopairing(self) -> bool:
        """
        Indicates whether the motif automatically pairs bases.
        """
        return self._autopairing

    @autopairing.setter
    def autopairing(self, autopairing: bool) -> None:
        """
        Set the autopairing property.
        """
        self._autopairing = bool(autopairing)
        if autopairing:
            self._basepair = BasePair()
            self._calculate_basepair()

    @property
    def basepair(self) -> BasePair:
        """
        A dictionary with positions as key and the paired position as values.
        """
        if self._autopairing and (
            not self._positions or not self._sequence or not self._basepair
        ):
            # calculate the basepair dictionary
            self._calculate_basepair()
            # reset the dot bracket
            self._structure = None

        return self._basepair

    @basepair.setter
    def basepair(
        self, basepair_dict: Union[Dict[Position, Position], BasePair]
    ) -> None:
        """
        Set the basepair dictionary and turn autopairing off.

        Parameters
        ----------
        basepair_dict : dict or BasePair
            A dictionary with positions as key and the paired position as values.
        """

        if not isinstance(basepair_dict, (dict, BasePair)):
            raise ValueError(
                f"{basepair_dict} must be a dictionary or a BasePair"
                f" object. Got {type(basepair_dict)} instead."
            )
        try:
            basepair_dict = {Position(k): Position(v) for k, v in basepair_dict.items()}
        except Exception as e:
            raise ValueError(
                "Error converting the basepair dictionary to"
                f" a dictionary of Position objects: {e}"
            )

        self._autopairing = False
        self._basepair = BasePair(basepair_dict, callback=self._updated_basepair)
        self._trigger_callbacks()

    @property
    def junctions(self) -> Dict[Tuple[int, int], Tuple[Tuple[int, int]]]:
        """
        A dictionary with Direction as key (right: (1,0), bottom: (0,1),
        left: (-1,0), top: (0,1)) and an ordered list of junction positon as value.
        """
        ### if the junctions are already calculated return them
        if self._positions is None:
            self._calculate_positions()
        return self._junctions

    @property
    def lock_coords(self) -> bool:
        """
        Boolean indicating if the coordinates are locked in the same block.
        """
        return self._lock_coords

    @lock_coords.setter
    def lock_coords(self, lock_coords):
        """
        Set wether the strands coordinates must be locked in the same block.
        """
        self._lock_coords = bool(lock_coords)

    @property
    def max_pos(self) -> Tuple[int, int]:
        """
        The maximum x, y coordinates occupied of the motif.
        """
        if not self._positions:
            self._calculate_positions()
        return self._max_pos

    @property
    def min_pos(self) -> Tuple[int, int]:
        """
        The minimum x, y coordinates occupied of the motif.
        """
        if not self._positions:
            self._calculate_positions()
        return self._min_pos

    @property
    def num_char(self) -> int:
        """
        The maximum length of the motif lines.
        """
        if not self._positions:
            self._calculate_positions()
        return self._max_pos[0] + 1

    @property
    def num_lines(self) -> int:
        """
        The number of lines in the motif structure.
        """
        if not self._positions:
            self._calculate_positions()
        return self._max_pos[1] + 1

    @property
    def pair_map(self) -> Dict[int, Union[int, None]]:
        """
        The dictionary of the paired indexes (alternative to the dot bracket
        notation).
        """
        if self._positions and self._pair_map:
            return self._pair_map

        pos_to_ind = {pos: ind for ind, pos in enumerate(self.seq_positions)}
        self._pair_map = BasePair(
            {
                ind: pos_to_ind.get(self.basepair.get(pos))
                for pos, ind in pos_to_ind.items()
            }
        )

        return self._pair_map

    @property
    def positions(self) -> Tuple[Position]:
        """
        A tuple of the positions of each character of each strand in 2D space
        (x,y coordinates).
        """
        if not self._positions:
            self._calculate_positions()
        return self._positions

    @property
    def sequence(self) -> str:
        """
        Return the sequence of the motif.
        """
        if self._sequence:
            self._sequence

        ### calculate the sequence
        tot_seq = ""
        for s in self:
            if not s.sequence:  # skip the strands without sequence
                continue
            if s.directionality == "35":
                tot_seq += str(s.sequence)[::-1]
            else:
                tot_seq += str(s.sequence)
            tot_seq += "&"

        self._sequence = Sequence(tot_seq.strip("&"))  # remove separator at the end
        return self._sequence

    @sequence.setter
    def sequence(self, seq_list: Union[str, List[str], Sequence] = None) -> None:
        """
        Set the sequence of each strand the motif

        Parameters
        ----------
        seq_list : str or Sequence or list of str or Sequence
            The list of sequences to set to the strands.
            If a single string is passed, it is split by '&'.

        Raises
        ------
        ValueError
            If the number of sequences is different from the number of strands.
        """
        if not isinstance(seq_list, Iterable) or (
            not isinstance(seq_list, (str, Sequence))
            and not isinstance(seq_list, (tuple, list))
            and any(not isinstance(s, (str, Sequence)) for s in seq_list)
        ):
            raise ValueError(
                f"{seq_list} must be a string, a Sequence object or "
                "a list of strings or Sequence objects. "
                f"Got {type(seq_list)} instead."
            )

        if isinstance(seq_list, (str, Sequence)):
            seq_list = seq_list.split("&")

        if len(seq_list) != len(self):
            raise ValueError(
                "The number of sequences must be equal to the number of"
                f" strands. Got {len(seq_list)} sequences "
                f" for {len(self)} strands."
            )

        seq_list = seq_list[:]  # copy the list
        for s in self:
            seq = seq_list.pop(0)
            direct = 1 if s.directionality == "53" else -1
            s.sequence = seq[::direct]

    @property
    def seq_positions(self) -> Tuple[Position]:
        """
        The positions of each nucleotide in the motif sequence (x,y coordinates).
        The sequence has always the directionality 5' to 3'
        """
        if not self._positions:
            self._calculate_positions()
        return self._seq_positions

    @property
    def strands(self) -> List["Strand"]:
        """
        The list of strands in the motif.
        """
        return self._strands

    @property
    def structure(self) -> str:
        """
        Return the dot bracket representation of the motif.
        """
        ### if the dot bracket is already calculated return it
        if self._structure:
            return self._structure

        # calculate the break points of the strands
        break_points = [i for i, sym in enumerate(self.sequence) if sym == "&"]

        ### CREATE THE DOT BRACKET NOTATION ###
        dotbracket = pair_map_to_dot_bracket(
            self.pair_map, len(self.sequence.replace("&", ""))
        )

        ### ADD THE BREAK POINTS ###
        for i, bp in enumerate(break_points):
            dotbracket = dotbracket[:bp] + "&" + dotbracket[bp:]
        self._structure = dotbracket

        return self._structure

    @structure.setter
    def structure(self, structure: str) -> None:
        """Set the dot bracket notation of the motif, it updates the basepair."""
        if not isinstance(structure, str):
            raise ValueError(
                f"{structure} must be a string. " f"Got {type(structure)} instead."
            )
        if len(structure) != len(self.sequence):
            raise ValueError(
                f"The length of the dot bracket must be equal to the"
                f" length of the sequence."
                f"Got {len(structure)} for {len(self.sequence)}"
            )

        basepair_dict = BasePair()
        # get the dictionary of paired indexes without the separation symbol
        pair_map = dot_bracket_to_pair_map(structure.replace("&", ""))
        # get the list of base positions in simple order

        # iterate over the base positions
        for index, pos in enumerate(self.seq_positions):
            paired = pair_map[index]
            # the index of the base position is paired to something
            if paired is not None:
                basepair_dict[Position(pos)] = Position(self.seq_positions[paired])

        # set the basepair dictionary
        self.basepair = basepair_dict
        self._pair_map = pair_map
        self._structure = structure

    ###
    ### CLASS METHODS
    ###

    @classmethod
    def concat(
        cls,
        *motifs: List["Motif"],
        axis: Literal[0, 1] = 1,
        extend: bool = False,
        copy: bool = True,
        align: bool = True,
        position_based: bool = False,
        align_junctions: Optional[List[Tuple[int, int]]] = None,
        unlock_strands: bool = False,
        return_shifts: bool = False,
        **kwargs,
    ) -> "Motif":
        """
        Concatenate multiple motifs along a specified axis.

        Parameters
        ----------
        *motifs : list of Motif
            List of motifs to concatenate.
        axis : int, default 1
            The numpy axis along which motifs are aligned.
            1 for horizontal alignment, 0 for vertical alignment.
        extend : bool, default False
            Whether to extend junctions before concatenation.
        copy : bool, default True
            Whether to copy motifs before concatenation.
        align : bool, default True
            Whether to align motifs before concatenation.
        position_based : bool, default False
            Whether to align motifs based on positions instead of junctions.
        align_junctions : list of tuple, optional
            List of junction indices for alignment.
        unlock_strands : bool, default False
            Whether to unlock all strands, so they are not part
            of the same Strand block (different coordinate systems).
        return_shifts : bool, default False
            Whether to return the shift values for each motif to be aligned
            and concatenated. If True, a list of tuples with the shifts is returned,
            along with the concatenated motif (useful for Origami alignment).

        **kwargs : dict
            Additional arguments to pass to the Motif constructor.

        Returns
        -------
        Motif
            The concatenated Motif object.
        List[Tuple[int, int]], optional
            The list of shifts applied to each motif during alignment.
            Only returned if return_shifts is True.
        """
        if (
            len(motifs) == 1
            and isinstance(motifs[0], Iterable)
            and not isinstance(motifs[0], Motif)
        ):
            # if a single list of motifs is passed, unpack it
            motifs = tuple(motifs[0])

        extra = tuple(kwargs.pop("motifs", tuple()))
        if isinstance(extra, Motif):
            motifs += (extra,)
        elif isinstance(extra, Iterable):
            motifs += tuple(extra)

        # check if all the elements in the motifs input are of type Motif
        if not all(isinstance(m, Motif) for m in motifs):
            raise ValueError(
                "All the elements in the motifs input must be of type"
                f" Motif. Got types {[type(m) for m in motifs]}."
            )

        # make a motif list and remove the empty motifs
        motifs = [m for m in motifs if m]

        if copy:
            motifs = [m.copy() for m in motifs]

        if return_shifts:
            # create a list of shifts for each motif
            shifts = [Position.zero()] * len(motifs)

        if align:
            aligned = Motif.align(
                *motifs,
                axis=axis,
                align_junctions=align_junctions,
                return_shifts=return_shifts,
            )
            if return_shifts:
                aligned, shifts = aligned

        else:
            aligned = motifs

        # trick to handle all the axis at once
        x_pos = axis
        y_pos = int(not axis)

        if extend:
            # Extend the junctions
            extend_until = [list(m.max_pos) for m in aligned]
            max_extend = max(extend_until, key=lambda x: x[axis], default=[0, 0, 0])[
                axis
            ]
            for extend_to in extend_until:
                extend_to[axis] = max_extend
            aligned = [
                m.extend_junctions(until=extend_to)
                for m, extend_to in zip(aligned, extend_until)
            ]

        # prepare the motif shifting them two by two
        for ind, m1 in enumerate(aligned[:-1]):
            m2 = aligned[ind + 1]

            # calculate the shift based on the positions
            if position_based:
                max_pos_m1 = m1.max_pos[y_pos]
                min_pos_m2 = m2.min_pos[y_pos]

            # calculate the shift based on the junctions
            else:
                m1_junct = m1.junctions[Position((x_pos, y_pos))]
                m2_junct = m2.junctions[Position((-x_pos, -y_pos))]

                if m1_junct and m2_junct:
                    ind1, ind2 = 0, 0
                    if align_junctions and align_junctions[ind]:
                        ind1, ind2 = align_junctions[ind]
                    max_pos_m1 = m1_junct[ind1][y_pos]
                    min_pos_m2 = m2_junct[ind2][y_pos]

                # No junctions, just place the motifs one after the other
                else:
                    max_pos_m1 = m1.max_pos[y_pos]
                    min_pos_m2 = 0
                    if m2.min_pos[y_pos] > max_pos_m1:
                        # if the second motif is already shifted, don't shift
                        max_pos_m1 = -1

            # Calculate the shift for the axis
            shift_pos = Position(
                (
                    x_pos * (max_pos_m1 - min_pos_m2 + 1),
                    y_pos * (max_pos_m1 - min_pos_m2 + 1),
                )
            )

            # Apply the shift
            m2.shift(shift_pos)

            # update the shifts
            if return_shifts:
                shifts[ind + 1] = shifts[ind + 1] + shift_pos

        ### BASEPAIRS HANDLING
        basepair = BasePair()
        # if all the motifs have autopairing on, the new motif has autopairing on
        autopairing = all([m.autopairing for m in aligned])

        if not autopairing:
            for m in aligned:
                basepair.update(m.basepair)

        ### STRANDS LOCKING
        if unlock_strands:
            for m in aligned:
                for s in m:
                    s.strands_block = None

        # create a object of the class
        new_motif = cls(**kwargs)
        new_motif._autopairing = autopairing
        new_motif._basepair = basepair
        new_motif.replace_all_strands(
            [strand for m in aligned for strand in m], copy=False, join=True
        )

        # used in Origami assembly
        if return_shifts:
            return new_motif, shifts

        return new_motif

    @classmethod
    def from_file(cls, file_path: str, **kwargs) -> "Motif":
        """
        Create a Motif object from a text file containing a motif sketch.
        Each Strand is read starting from the `5` symbol. If you want to add the 5'
        terminal symbol, start the strand with `55`. Only one symbol should be placed
        next to the `5` start of the strand in order to guess the right direction.
        Alternatively, you can use symbols `^v><` to start a strand and indicate
        the start direction (up, down, right, left); in this case you can place
        multiple symbols next to the start symbol.

        Parameters
        ----------
        file_path : str
            Path to the file containing the motif structure.
        **kwargs : dict
            Additional arguments to pass to the Motif constructor.

        Returns
        -------
        Motif
            The constructed Motif object.
        """
        with open(file_path, "r", encoding="utf-8") as f:
            motif_text = f.read()
        return cls.from_text(motif_text, **kwargs)

    @classmethod
    def from_list(cls, motif_list: List[str], **kwargs) -> "Motif":
        """
        Create a Motif object from a list of strings representing a motif sketch.
        Each Strand is read starting from the `5` symbol. If you want to add the 5'
        terminal symbol, start the strand with `55`. Only one symbol should be placed
        next to the `5` start of the strand in order to guess the right direction.
        Alternatively, you can use symbols `^v><` to start a strand and indicate
        the start direction (up, down, right, left); in this case you can place
        multiple symbols next to the start symbol.

        Parameters
        ----------
        motif_list : list of str
            The list of strings representing the motif.
        **kwargs : dict
            Additional arguments to pass to the Motif constructor.

        Returns
        -------
        Motif
            The constructed Motif object.
        """
        ### initialize the variables
        strand_list = []  # list of strands
        mapped_pos = set()  # set of visited positions

        ### Get the maximum x and y positions
        max_x = max([len(line) for line in motif_list], default=1)
        max_y = len(motif_list) - 1

        ### Uniform the line length
        for index, line in enumerate(motif_list):
            diff = max_x - (len(line) - 1)
            if diff:
                motif_list[index] = line + " " * diff

        ### TRACE THE STRANDS STARTING WITH ^v><
        for y, line in enumerate(motif_list):
            for x, char in enumerate(line):
                # check a start of a strand
                if char in "^v><":  # Special direction characters
                    if char == "^":
                        direction = (0, -1)
                    elif char == "v":
                        direction = (0, 1)
                    elif char == ">":
                        direction = (1, 0)
                    elif char == "<":
                        direction = (-1, 0)

                    # get the strand direction
                    start_pos = (x + direction[0], y + direction[1])
                    # get the strand characters
                    strand_chars = Motif.trace(
                        motif_list, start_pos, direction, limits=(max_x, max_y)
                    )
                    # create the Strand object
                    strand = Strand(
                        strand_chars,
                        directionality="53",
                        start=start_pos,
                        direction=direction,
                    )
                    # add the Strand to the list
                    strand_list.append(strand)
                    # update the visited positions
                    mapped_pos.update(set(strand.positions))

        ### TRACE THE STRANDS STARTING WITH 5
        for y, line in enumerate(motif_list):
            for x, char in enumerate(line):
                # check a start of a strand
                if char == "5" and (x, y) not in mapped_pos:
                    right_sym = None
                    bot_sym = None
                    left_sym = None
                    top_sym = None

                    # try to guess the direction based on the next symbol
                    if x + 1 <= max_x:
                        right_sym = motif_list[y][x + 1].capitalize()
                    if y + 1 <= max_y:
                        bot_sym = motif_list[y + 1][x].capitalize()
                    if x > 0:
                        left_sym = motif_list[y][x - 1].capitalize()
                    if y > 0:
                        top_sym = motif_list[y - 1][x].capitalize()

                    direction = None
                    if right_sym in accept_symbol and motif_list[y][x + 1] not in " 3":
                        direction = (1, 0)
                    elif bot_sym in accept_symbol and bot_sym not in " 3":
                        direction = (0, 1)
                    elif left_sym in accept_symbol and left_sym not in " 3":
                        direction = (-1, 0)
                    elif top_sym in accept_symbol and top_sym not in " 3":
                        direction = (0, -1)

                    if direction is None:
                        raise MotifStructureError(
                            f"Strand at position {x, y}" " has no direction."
                        )

                    # get the strand direction
                    start_pos = (x + direction[0], y + direction[1])
                    # get the strand characters
                    strand_chars = Motif.trace(
                        motif_list, start_pos, direction, limits=(max_x, max_y)
                    )
                    # create the Strand object
                    strand = Strand(
                        strand_chars,
                        directionality="53",
                        start=start_pos,
                        direction=direction,
                    )
                    # add the Strand to the list
                    strand_list.append(strand)
                    # update the visited positions
                    mapped_pos.update(set(strand.positions))

        ### TRACE THE BASEPAIRS ###
        basepair = BasePair()
        for y, line in enumerate(motif_list):
            for x, sym in enumerate(line):
                pos1 = None

                # if found a basepair symbol, add the position to the dictionary
                if sym in "┊!*:":  # vertical basepair
                    pos1 = (x, y - 1)
                    pos2 = (x, y + 1)

                elif sym in "=":  # horizontal basepair
                    pos1 = (x - 1, y)
                    pos2 = (x + 1, y)

                if pos1 is not None:
                    basepair[Position(pos1)] = Position(pos2)

        return cls(strands=strand_list, basepair=basepair, **kwargs)

    @classmethod
    def from_structure(
        cls,
        structure: Optional[Union[str, dict, BasePair, Node]] = None,
        sequence: Optional[str] = None,
        pk_energy=-8.5,
        pk_denergy=0.5,
        **kwargs,
    ) -> "Motif":
        """
        Parse a structure or sequence representation to a Motif object.
        If a structure is not provided, it is calculated from the sequence
        with RNAfold. If a sequence is not provided, it is assumed to be
        a sequence of 'N's of the same length as the structure.

        Parameters
        ----------
        structure : Union[str, dict, BasePair, Node]
            The structure representation to convert.
        sequence : str, optional
            The sequence or sequence constraints of the motif.
        pk_energy : float, optional
            The energy of the pseudoknots (if present).
        pk_denergy : float, optional
            The energy tolerance of the pseudoknots (if present).
        **kwargs : dict
            Additional arguments to pass to the Motif constructor.


        Returns
        -------
        Motif
            The Motif object created from the structure representation.
        """
        # import here to avoid circular imports
        from .origami import Origami
        from RNA import fold

        if not structure:
            # if only sequence is provided, fold it to get the structure
            structure = fold(sequence)[0]
            for i, sym in enumerate(sequence[::-1]):
                if sym == "&":
                    structure = structure[:-i] + "&" + structure[-i:]

        # input dot-bracket notation
        if isinstance(structure, str):
            node = dot_bracket_to_tree(structure, sequence=sequence)
            pair_map = dot_bracket_to_pair_map(structure)
        # input pair map
        elif isinstance(structure, (BasePair, dict)):
            pair_map = structure.copy()
            node = dot_bracket_to_tree(
                pair_map_to_dot_bracket(structure), sequence=sequence
            )
            structure = pair_map_to_dot_bracket(structure)
        # input tree
        elif isinstance(structure, Node):
            node = structure
            pair_map = dot_bracket_to_pair_map(tree_to_dot_bracket(node))
            structure = tree_to_dot_bracket(node)
        else:
            raise ValueError(f"Invalid structure representation: {structure}")

        if not sequence:
            sequence = "".join("N" if sym != "&" else "&" for sym in structure)
        else:
            sequence = str(sequence).replace("T", "U").upper()

        if isinstance(structure, str) and len(structure.strip("& ")) != len(
            sequence.strip("& ")
        ):
            raise ValueError(
                f"The sequence length must be equal to the structure "
                f"length. Got sequence len {len(sequence)} for structure"
                f" len {len(structure)}."
            )

        # initialize the origami object
        origami = Origami([[]], align="first", ss_assembly=True)
        current_index = [0, 0]

        def recursive_build_origami(node, insert_at=None, flip=False):
            """
            Recursively build the origami from the tree representation.
            """
            # nonlocal current_index

            # initialize the variables
            if insert_at is None:
                insert_at = current_index
            motif = None

            ### BASE CASES: sequence break, unpaired nucleotide, stem
            if node.label == "&":
                return

            if node.label == "(":
                motif = Motif(
                    Strand(node.seq if node.seq else "N"),
                    Strand(
                        sequence[pair_map[node.index]],
                        start=(0, 2),
                        directionality="35",
                    ),
                    basepair={(0, 0): (0, 2)},
                )

            elif node.label == ".":
                motif = Motif(
                    Strand(node.seq if node.seq else "N"), Strand("-", start=(0, 2))
                )

            # add the motif and update the current index
            if motif:
                origami.insert(insert_at, motif.flip(flip, flip))
                current_index[1] += 1  # increment the x index

            # recursive call for the children
            if node.children:
                child_inds = []

                # check each child before recursive call
                for i, child in enumerate(node.children):
                    insert_at = None
                    flip = False

                    # bulge after a stem
                    if child.label == "." and any(
                        c.label == "(" for c in node.children[:i]
                    ):
                        insert_at = child_inds.pop()
                        flip = True

                    # sequence break + only unpaired
                    elif child.label in ".&" and all(
                        c.label in ".&" for c in node.children
                    ):
                        if "&" in [c.label for c in node.children[:i]]:
                            insert_at = child_inds.pop()
                            flip = True

                    # sequence break or multiple stems
                    elif child.label == "&" or (
                        child.label == "("
                        and any(c.label == "(" for c in node.children[:i])
                    ):
                        connect_down = Motif(
                            Strand("──"),
                            Strand("╮", start=(0, 2), directionality="35"),
                            Strand("╭", start=(1, 2), direction=(0, -1)),
                        )
                        connect_up = Motif(
                            Strand("││╰─", direction=(0, 1), directionality="35"),
                            Strand("╰", start=(1, 0), direction=(0, 1)),
                        )
                        if child_inds:
                            insert_connect = child_inds.pop()
                        else:
                            insert_connect = current_index

                        # insert the top connector
                        origami.insert(insert_connect, connect_down)

                        shift_x = sum(
                            [
                                m.num_char
                                for m in origami[insert_connect[0], : insert_connect[1]]
                            ]
                        )

                        connect_up.shift((shift_x, 0))
                        origami.append([connect_up])

                        # increment the y index
                        current_index[0] += 1
                        # set the x index to the end of the line
                        current_index[1] = len(origami[-1])

                        for i in range(insert_connect[0] + 1, current_index[0]):
                            # add the vertical connector
                            origami.insert(
                                (i, 0),
                                Motif(
                                    Strand("│", direction=(0, 1), directionality="35"),
                                    Strand("│", direction=(0, 1), start=(1, 0)),
                                ).shift((shift_x, 0)),
                            )
                            # shift all the motifs until you reach the first connector
                            for m in origami[i, 1:]:
                                m.shift((2, 0))
                                if "││╰─" in m:
                                    break

                    if insert_at is None:
                        insert_at = current_index.copy()

                    child_inds.append(insert_at)
                    recursive_build_origami(child, insert_at=insert_at, flip=flip)

                # this could not work in the case a stem doesn't end with at least
                # one unpaired nucleotide, but that does never happen in natural
                # structures, so we can ignore this case
                if not any(c.children or c.label == "&" for c in node.children):
                    origami.append(Motif(Strand("╮│╯")))
                    current_index[1] -= 1  # decrement the x index

        # call the recursive function
        recursive_build_origami(node)
        # get the motif from the origami object
        motif = origami.assembled

        ### ADD THE PSEUDOKNOTS ###
        seq_offset = 0
        # dictionary with index as key and pseudoknot id as value
        full_map = dict()
        pair_map = dot_bracket_to_pair_map(structure.replace("&", ""))

        # iterate over the  subsequences
        for i, struct in enumerate(structure.split("&")):
            new_pk_info = {"id": [], "ind_fwd": [], "E": [], "dE": []}
            j = 0
            ss_len = len(struct)

            # iterate over the subsequences structure
            while j < ss_len:
                sym = struct[j]

                # found pseudoknot
                if sym not in ".()" and (seq_offset + j) not in full_map:
                    # get the length of the pseudoknot
                    length = 1
                    while struct[j + length] == sym:
                        length += 1

                    # get the pseudoknot id of get a new one
                    if pair_map[seq_offset + j] in full_map:
                        pk_id = full_map[pair_map[seq_offset + j]] + "'"
                    else:
                        inds = [k.split("_")[1].strip("'") for k in full_map.values()]
                        pk_id = "100_" + str(int(max(inds, default="-1")) + 1)

                    # add the pseudoknot to the motif
                    new_pk_info["id"].append(pk_id)
                    new_pk_info["ind_fwd"].append((j, j + length - 1))
                    indices = range(seq_offset + j, seq_offset + j + length)
                    # update the full map
                    full_map.update({k: pk_id for k in indices})

                    new_pk_info["E"].append(pk_energy)
                    new_pk_info["dE"].append(pk_denergy)
                    j += length
                j += 1

            # add the pseudoknots info to the strand
            motif[i].pk_info = new_pk_info
            seq_offset += len(struct)

        kwargs.setdefault("lock_coords", False)
        obj = cls(**kwargs)
        obj.replace_all_strands(motif, copy=False, join=False)
        return obj

    @classmethod
    def from_text(cls, motif_text: str, **kwargs) -> "Motif":
        """
        Create a Motif object from a text string representing a motif sketch.
        Each Strand is read starting from the `5` symbol. If you want to add the 5'
        terminal symbol, start the strand with `55`. Only one symbol should be placed
        next to the `5` start of the strand in order to guess the right direction.
        Alternatively, you can use symbols `^v><` to start a strand and indicate
        the start direction (up, down, right, left); in this case you can place
        multiple symbols next to the start symbol.

        Parameters
        ----------
        motif_text : str
            The motif structure in text format.
        **kwargs : dict
            Additional arguments to pass to the Motif constructor.

        Returns
        -------
        Motif
            The constructed Motif object.
        """
        motif_list = [line for line in motif_text.split("\n")]
        return cls.from_list(motif_list, **kwargs)

    ###
    ### STATIC METHODS
    ###

    @staticmethod
    def align(
        *motifs: List["Motif"],
        axis: Literal[0, 1] = 1,
        extend: bool = False,
        align_junctions: Optional[List[Tuple[int, int]]] = None,
        align_to: Literal["first", "last", "center"] = "first",
        return_shifts: bool = False,
        **kwargs,
    ) -> List["Motif"] | Tuple[List["Motif"], List[Tuple[int, int]]]:
        """
        Align motifs along a given axis by shifting them (without concatenating).

        Parameters
        ----------
        *motifs : List[Motif]
            List of motifs to align.
        axis : int, default 1
            The numpy axis along which motifs are aligned
            (1 for horizontal, 0 for vertical).
        extend : bool, default False
            If True, junctions are extended to accommodate the shift.
        align_junctions : Optional[List[Tuple[int, int]]], default None
            List of tuples specifying which junctions to align.
        align_to : str, default "first"
            Specifies how to align the motifs. Options: "first", "last", "center".
            "first" aligns the first junction of the first motif with the first junction
            of the second motif. "last" aligns the last junction of the first motif
            with the last junction of the second motif. "center" aligns the center
            of the first motif with the center of the second motif.
        return_shifts : bool, default False
            If True, returns the shifts applied to each motif during alignment,
            along with the aligned motifs (useful for Origami alignment).
            If False, only the aligned motifs are returned.
        **kwargs : dict
            Parse if calling the motifs list explicitly.

        Returns
        -------
        List[Motif]
            The list of aligned motifs.
        List[Tuple[int, int]], optional
            The list of shifts applied to each motif during alignment.
            Only returned if return_shifts is True.
        """

        ### get the alignment direction
        if axis in (1, 0):
            direction = Position((axis, int(not axis)))
        else:
            raise ValueError(
                f"{axis} is not a valid value for the axis parameter."
                " The axis parameter must be 0 or 1"
            )

        # the direction to shift is the opposite of the direction of the axis
        shift_direction = Position((int(not axis), axis))

        ### check alignment type
        if align_to not in ("first", "last", "center"):
            raise ValueError(
                f"{align_to} is not a valid value for the align_to"
                "parameter. The align_to parameter must be"
                ' "first", "last" or "center"'
            )

        ### check if the motifs are of the right type
        if (
            len(motifs) == 1
            and isinstance(motifs[0], Iterable)
            and not isinstance(motifs[0], Motif)
        ):
            # an iterable of motifs is passed
            motifs = tuple(motifs[0])

        extra = tuple(kwargs.pop("motifs", tuple()))
        if isinstance(extra, Motif):
            motifs += (extra,)
        elif isinstance(extra, Iterable):
            motifs += tuple(extra)

        # check if all the motifs are of type Motif
        if not all(isinstance(m, Motif) for m in motifs):
            raise ValueError("All the motifs must be of type Motif.")

        # if shifts are requested, save them
        if return_shifts:
            shifts = [Position((0, 0)) for _ in motifs]

        ### align all the motifs
        ind1 = 0
        ind2 = 1
        n_motifs = len(motifs)
        # start with the assumption that all the motifs are connected
        connected_motifs = [True] * n_motifs

        while ind2 < n_motifs:  # stop when the second motif is the last one
            m1 = motifs[ind1]
            m2 = motifs[ind2]

            if align_junctions and align_junctions[ind1]:
                # align the junction with the given index
                junct_ind1, junct_ind2 = align_junctions[ind1]

            else:
                # align the first junction of the first motif
                # with the first junction of the second motif
                junct_ind1, junct_ind2 = 0, 0  # by default align the first junction

                if align_to == "last":
                    junct_ind1, junct_ind2 = -1, -1

            junctions1 = m1.junctions[direction]
            junctions2 = m2.junctions[-direction]

            ### calculate the shift of the two motifs
            # and shift the motif with the lower junction
            if align_to == "center":
                # if align_to is center align to the center of the motifs
                shift = int(
                    (m1.max_pos[axis] + m1.min_pos[axis]) / 2
                    - (m2.max_pos[axis] + m2.min_pos[axis]) / 2
                )

            elif (
                junctions1
                and junctions2
                and junct_ind1 < len(junctions1)
                and junct_ind2 < len(junctions2)
            ):
                # if junctions are detected, calculate the shift based on the junctions
                shift = junctions1[junct_ind1][axis] - junctions2[junct_ind2][axis]

            else:  # nothing to align, the motifs are not connected
                connected_motifs[ind1] = False  # mark the first motif as not connected
                shift = 0

            shift_pos = shift_direction * shift

            if shift > 0:
                # shift only the newly appended motif
                m2.shift(shift_pos, extend=extend)

                # update the shifts
                if return_shifts:
                    shifts[ind2] = shifts[ind2] + shift_pos

            elif shift < 0:
                shift_pos = -shift_pos
                # shift the left motif and all the previous one if they are connected
                for j in range(ind1, -1, -1):
                    if not connected_motifs[j]:
                        break

                    motifs[j].shift(shift_pos, extend=extend)

                    # update the shifts
                    if return_shifts:
                        shifts[j] = shifts[j] + shift_pos

            # update the indices
            ind1 += 1
            ind2 += 1

        if return_shifts:
            # if the shifts are requested, return them
            return motifs, shifts

        # if the shifts are not requested, return the motifs
        return motifs

    @staticmethod
    def copy_strands_preserve_blocks(
        strands: List["Strand"], motif: Optional["Motif"] = None
    ) -> List["Strand"]:
        """
        Copy a list of strands while preserving strand blocks.
        Ensures that strands that belong to the same motif remain linked together.

        Parameters
        ----------
        strands : list of Strand
            The list of strands to be copied.
        motif : Motif, optional
            If provided, associates the strands with the given motif and registers
            callback functions.

        Returns
        -------
        list of Strand
            A list of copied strands.
        """
        # IMPORTANT: keep strand that are part of the same motif linked
        # get the set of strands block id
        motifs_id = {id(s.strands_block) for s in strands}
        # for each strands block, make a new one and link to the old id
        new_motifs_dict = {key: StrandsBlock() for key in motifs_id}

        strands_copy = []
        # copy all strands
        for strand in strands:
            # check that the strands are in the strand class
            if not isinstance(strand, Strand):
                raise ValueError(f"{strand} is not a Strand object.")
            # collect the new strands block
            new_strands_block = new_motifs_dict[id(strand.strands_block)]
            # copy the strand and add the callback
            if motif is not None:
                copied = strand.copy(callback=motif._updated_strands)
            else:
                copied = strand.copy()
            strands_copy.append(copied)
            new_strands_block.add(copied)  # add the origami to the new strands block
        return strands_copy

    @staticmethod
    def get_sequential_shift(
        motifs: List["Motif"],
        *args,
        axis: Literal[0, 1] = 1,
        position_based: bool = True,
    ) -> List[int]:
        """
        Calculate the shift needed to align motifs sequentially.

        Parameters
        ----------
        motifs : List[Motif]
            List of motifs to align.
        *args: List[Motif]
            Additional motifs to align.
        axis : int, default 1
            The numpy axis along which the shift is computed
            (1 for horizontal, 0 for vertical).
        position_based : bool, default True
            If True, the shift is calculated based on motif positions
            rather than junctions.

        Returns
        -------
        List[int]
            A list of shift integers for each motif along the axis.

        Examples
        --------
        .. code-block:: text

            m1:
                --NN--
                  ::
                --NN--
            m2:
                --SK--
                  ::
                --SK--
            Resulting shift: 2
                [0, 2]
            To have:
                --NN----SK--
                  ::    ::
                --NN----SK--
        """
        if isinstance(motifs, Motif):
            motifs = [motifs]
        # put the motifs in a list and remove the empty motifs
        motifs = [m for m in motifs + list(args) if m]

        # trick to handle all the axis at once
        x_pos = axis
        y_pos = int(not axis)

        shifts = [0]
        current_shift = 0

        for ind, m1 in enumerate(motifs[:-1]):
            m2 = motifs[ind + 1]

            # calculate the shift based on the positions
            if position_based:
                max_pos_m1 = m1.max_pos[y_pos]
                min_pos_m2 = m2.min_pos[y_pos]

            # calculate the shift based on the junctions
            else:
                m1_junct = m1.junctions
                m2_junct = m2.junctions

                # if the junctions are detected, calculate the shift
                # based on the junctions
                junc_dir = Position((x_pos, y_pos))
                if m1_junct[junc_dir] and m2_junct[-junc_dir]:
                    max_pos_m1 = max((pos[y_pos] for pos in m1_junct[junc_dir]))
                    min_pos_m2 = min((pos[y_pos] for pos in m2_junct[-junc_dir]))

                else:  # no junctions, just concatenate the motifs
                    max_pos_m1 = m1.max_pos[y_pos]
                    min_pos_m2 = 0

            # update the shift
            current_shift += max_pos_m1 - min_pos_m2 + 1
            shifts.append(current_shift)

        return shifts

    @staticmethod
    def join_strands(strands: List["Strand"]) -> List["Strand"]:
        """
        Attempt to join consecutive strands and return the list of joined strands.

        Parameters
        ----------
        strands : List[Strand]
            The list of strands to be joined.

        Returns
        -------
        List[Strand]
            The list of joined strands.
        """
        # remove the empty strands
        strands = [s for s in strands if s]

        # ordering strand is:
        # - the strands with 5' end
        # - strand with the lowest y start position
        # - strand with the lowest x start position
        strands = sorted(strands, key=lambda s: (-int("5" in s.strand), *s.start[::-1]))

        # Comments for coders: I tried different ways to join the strands,
        # especially building a graph of joinable strands then using Depth First Search
        # to join the adjacent strands; but this version of the code is faster

        joined_strands = set()
        for ind1, s1 in enumerate(strands[:-1]):
            # if the strand is empty, or is already joined skip it
            if ind1 in joined_strands:
                continue

            # join the strand with the consecutive strands
            # until no more strands are joined
            while True:
                current_len = len(joined_strands)
                for ind2, s2 in enumerate(strands):

                    # skip s2 is the same strand or already joined
                    if ind1 == ind2 or ind2 in joined_strands:
                        continue

                    # skip if the motifs are not joinable for sure
                    if (
                        len(
                            set(
                                (
                                    s1.prec_pos,
                                    s1.start,
                                    s1.end,
                                    s1.next_pos,
                                    s2.prec_pos,
                                    s2.start,
                                    s2.end,
                                    s2.next_pos,
                                )
                            )
                        )
                        == 8
                    ):
                        continue

                    # try to join the strands
                    joined = Strand.join_strands(s1, s2)

                    # the strands are joined
                    if joined is not None:
                        # add the index of strand2 to the joined strands
                        joined_strands.add(ind2)
                        # replace the first strand with the joined strand
                        strands[ind1] = joined
                        # update the first strand
                        s1 = joined

                # if nothing is joined, break the loop
                if current_len == len(joined_strands):
                    break

        # return the strands excluded the joined ones
        return [s for i, s in enumerate(strands) if i not in joined_strands]

    @staticmethod
    def trace(
        motif_list: List[str],
        pos: Tuple[int, int],
        direction: Tuple[int, int],
        limits: Tuple[int, int],
    ) -> str:
        """
        Trace strands in a motif representation in a recursive manner.

        Parameters
        ----------
        motif_list : list of str
            The list of strings of the motif sketch.
        pos : tuple of int
            The starting position.
        direction : tuple of int
            The tracing direction.
        limits : tuple of int
            The boundary limits.

        Returns
        -------
        str
            The traced strand sequence.
        """
        # pos is a tuple (x,y) of the character we are considering
        # so it is a tuple (char_ind, line_ind)
        x = pos[0]
        y = pos[1]

        # We reached an adge
        if x > limits[0] or x == -1 or y > limits[1] or y == -1:
            return ""

        # Read the new symbol
        sym = motif_list[y][x]

        # Terminal symbols
        if sym == "3":
            return sym
        elif sym == " " or sym not in accept_symbol:
            return ""

        # Direction turns
        elif sym in "╰╮\\":  # reverse the tuple
            direction = direction[::-1]
        elif sym in "╭╯/":  # reverse the tuple and change sign
            direction = (-direction[1], -direction[0])

        # Error
        elif sym in bp_symbols:
            raise MotifStructureError("The strand leads to a base pairing")

        # Go to the next position
        pos = (pos[0] + direction[0], pos[1] + direction[1])

        # Recursive call
        return sym + Motif.trace(motif_list, pos, direction, limits)

    ###
    ### PROTECTED METHODS
    ###

    def _calculate_basepair(self) -> None:
        """
        Calculate and store base pairings in the `_basepair` attribute.
        The method determines base pairing by identifying complementary bases
        that are positioned one step away (horizontally or vertically). The
        base pairing dictionary (using positions as key/value) is stored in
        the `_basepair` attribute as a Basepair object.
        """
        basepair = BasePair(callback=self._updated_basepair)
        pos_to_ind = {pos: ind for ind, pos in enumerate(self.seq_positions)}
        sequence = str(self.sequence).replace("&", "")

        for pos1, ind1 in pos_to_ind.items():
            base1 = sequence[ind1]
            # skip already paired pos
            if pos1 in basepair:
                continue

            # control if there is a complementary base in all the directions
            for d in Direction:
                # calculate the second position
                pos2 = pos1 + d * 2

                # skip the positions already paired or not in the base map
                if pos2 in basepair or pos2 not in pos_to_ind:
                    continue

                # if the basepair position is already taken, skip
                bp_pos = pos1 + d
                if bp_pos in self._positions:
                    continue

                # check the bases can pair
                base2 = sequence[pos_to_ind[pos2]]
                if base2 in base_pairing[base1]:
                    basepair[Position(pos1)] = Position(pos2)
                    break

        self._basepair = basepair

    def _calculate_positions(self) -> None:
        """
        Calculate the 2D properties of the strand.
        This function calculates the properties:
            - positions
            - seq_positions
            - max_pos
            - min_pos
            - junctions
        """
        ### initialize the variables
        positions = ()
        seq_positions = ()
        max_pos = [0, 0]
        min_pos = [float("inf"), float("inf")]
        ### initialize the junctions
        junctions_dict = {direct: [] for direct in Direction}

        for ind, s in enumerate(self):
            if not s:
                continue
            positions += s.positions

            seq_dir = 1 if s.directionality == "53" else -1
            seq_positions += s._seq_positions[::seq_dir]

            # update MAX/MIN positions
            if s._max_pos[0] > max_pos[0]:
                max_pos[0] = s._max_pos[0]
            if s._max_pos[1] > max_pos[1]:
                max_pos[1] = s._max_pos[1]
            if s._min_pos[0] < min_pos[0]:
                min_pos[0] = s._min_pos[0]
            if s._min_pos[1] < min_pos[1]:
                min_pos[1] = s._min_pos[1]

            # skip if the the start is '5' or '3' and the strand is not just the symbol
            if s[0] not in "35":
                # invert the start direction to have the direction of the junction
                junctions_dict[-s.direction].append(s.start)

            # skip if the the end is '5' or '3' and the strand is not just the symbol
            if s[-1] not in "35":
                junctions_dict[s.end_direction].append(s.end)

        ### order the junctions
        for key, val in junctions_dict.items():
            # order the bottom/top junctions
            if key in (Direction.DOWN, Direction.UP):
                val.sort()  # order the junctions in growing x and y
            # order the left/right junctions
            elif key in (Direction.LEFT, Direction.RIGHT):
                val.sort(key=lambda pos: pos.swap_xy())  # order in growing y and x

        #### convert the junctions to tuple
        self._junctions = {key: tuple(val) for key, val in junctions_dict.items()}

        # save the other positional properties
        self._positions = positions
        self._seq_positions = seq_positions
        self._max_pos = Position(max_pos)
        self._min_pos = Position((c if c != float("inf") else 0 for c in min_pos))

    def _check_addition(
        self, other: "Motif", direction: Union[Position, Tuple[int, int]] = None
    ) -> None:
        """
        Check whether two motifs can be added together in a given direction.
        The function ensures that the motifs have compatible junctions for addition.
        It raises an error if they cannot be added.

        Parameters
        ----------
        other : Motif
            The motif to be added.
        direction : tuple of int, optional
            The direction in which the motifs should be checked for addition.
            Default is Direction(1, 0) (right).

        Raises
        ------
        ValueError
            If `other` is not an instance of `Motif`.
        MotifStructureError
            If the motifs do not have compatible junctions.
        """
        # if one of the two motifs is empty, return without error
        if not self or not other:
            return

        if not isinstance(other, Motif):
            raise ValueError(f"{other} is not a valid type for addition")

        if direction is None:
            direction = Direction.RIGHT
        elif not isinstance(direction, (Direction, Position)):
            direction = Position(direction)

        # take the junctions of the left and right side of the motifs
        Strand._check_position(input_pos_dir=direction, direction=True)
        junction1 = self.junctions[direction]
        junction2 = other.junctions[direction * -1]

        ### SANITY CHECKS
        if not junction1 or not junction2:
            raise MotifStructureError(
                f"The motifs cannot be added in the direction {direction}, "
                "the junctions are missing. Junctions"
                f" motif1: {junction1}, Junctions motif2: {junction2}. If you want"
                " to concatenate the motifs, use pf.Motif.concat() method."
            )

    def _updated_basepair(self, **kwargs) -> None:
        """
        Update the base-pair dictionary when the sequence changes.

        This is a callback function triggered when modifications occur to the sequence
        or base-pair mappings, ensuring consistency in structural relationships.

        Parameters
        ----------
        **kwargs
            Additional arguments that may be passed during the update process.

        Notes
        -----
        This function resets the dot-bracket notation representation of the motif.
        """
        self._structure = None  # reset the dot bracket
        self._pair_map = BasePair()
        self._trigger_callbacks(**kwargs)

    def _updated_strands(self, **kwargs) -> None:
        """
        Update the motif when the strands are changed.

        This is a callback function that is triggered when the strands of the motif
        are modified, ensuring that cached properties (e.g., structure, base pairs)
        remain consistent.

        Parameters
        ----------
        **kwargs
            Additional arguments that may be passed during the update process.

        Notes
        -----
        This method resets cached attributes such as:
        - Sequence
        - Position mappings
        - Base-pair structures
        """
        self._sequence = ""
        self._positions = None
        if self._autopairing:
            self._basepair = BasePair()
            self._structure = None
        self._updated_basepair()
        self._trigger_callbacks(**kwargs)

    ###
    ### METHODS
    ###

    def append(
        self, strand: "Strand", join: bool = True, copy: bool = False
    ) -> "Motif":
        """
        Append a strand to the motif.

        Parameters
        ----------
        strand : Strand
            The strand to be added.
        join : bool, default True
            Whether to attempt joining the new strand with existing strands.
        copy : bool, default False
            If True, a copy of the strand is appended.

        Returns
        -------
        Motif
            The updated motif.

        Raises
        ------
        ValueError
            If the provided `strand` is not a `Strand` instance.
        """

        if not isinstance(strand, Strand):
            raise ValueError(f"{strand} is not a Strand object.")

        if copy:
            strand = strand.copy()

        strand.register_callback(self._updated_strands)
        self._strands.append(strand)

        if join:
            self._strands = self.join_strands(self._strands)
        if self.lock_coords:
            self._strands_block = StrandsBlock(*[s for s in self._strands])

        self._updated_strands()
        return self

    def copy(self, callback=None, **kwargs) -> "Motif":
        """
        Return a copy of the motif.

        This method creates a custom deep copy of the current motif, preserving
        its strand arrangements, base pairings, and additional attributes.

        Parameters
        ----------
        callback : callable, default None
            A callback function to be registered in the copied motif.
        **kwargs : dict
            If the keyword 'deepcopy' with a list of attributes names is passed,
            the attributes will deepcopyed instead referenced.
            Same for the keyword 'copy' with a list of attributes names.

        Returns
        -------
        Motif or subclass
            A new instance of the motif or its subclass, with the same properties
            and attributes as the original.
        """
        new_motif = type(self).__new__(type(self))

        # attributes that are calculated fresh or immutable that
        # can be just reassigned
        basic_properties = {
            "_autopairing",
            "_lock_coords",
            "_structure",
            "_positions",
            "_seq_positions",
            "_max_pos",
            "_min_pos",
            "_sequence",
        }
        for prop in basic_properties:
            setattr(new_motif, prop, getattr(self, prop))

        # attributes that have to be teatred differently
        avoid_to_copy = {
            "_basepair",
            "_junctions",
            "_pair_map",
            "_strands_block",
            "_strands",
            "_callbacks",
        }

        new_motif._strands = self.copy_strands_preserve_blocks(self._strands, new_motif)
        if self._lock_coords and new_motif._strands:
            new_motif._strands_block = new_motif[0].strands_block
        new_motif._junctions = {k: tuple(v) for k, v in self._junctions.items()}
        new_motif._basepair = self._basepair.copy()
        new_motif._pair_map = self._pair_map.copy()
        new_motif._callbacks = [callback] if callback else []

        ### The other attributes should be freshly calculated
        # or immutable, so just reassign them
        for attr in self.__dict__:
            # deepcopy the attribute if it is in the deepcopy list
            if attr in kwargs.get("deepcopy", ()):
                setattr(new_motif, attr, copy.deepcopy(getattr(self, attr)))
                continue

            # copy the attribute if it is in the copy list
            elif attr in kwargs.get("copy", ()):
                setattr(new_motif, attr, getattr(self, attr).copy())
                continue

            # check if the attribute is to avoid or already copied
            if attr in avoid_to_copy or attr in basic_properties:
                continue

            # reference the attribute
            setattr(new_motif, attr, getattr(self, attr))

        return new_motif

    def extend_junctions(
        self,
        skip_axis: Literal[None, 0, 1] = None,
        skip_directions: List[Direction] = None,
        until: Optional[Position] = None,
    ) -> "Motif":
        """
        Extend the junctions of the motif in the direction of the junctions.
        Unless the skip_axis is specified, the junctions are extended
        in all directions. The until parameter allows to specify
        the position until which to extend the junctions in the positive direction.

        Parameters
        ----------
        skip_axis: int, optional
            The numpy axis to skip when extending the junctions
            (1 for horizontal, 0 for vertical).
        skip_directions: List[Direction], optional
            The list of directions to skip when extending the junctions.
        until : Position, optional
            The position until which to extend the junctions
            If None, extend until the end of the motif in all directions
            (to the origin and to the maximum position).

        Returns
        -------
        Motif
            The motif with the extended junctions.

        Raises
        ------
        ValueError
            If the axis parameter is not None, 0 or 1.
        MotifStructureError
            If there is an error in extending the junction
        """
        if skip_directions is None:
            skip_directions = []

        if skip_axis in (0, 1):
            skip_directions.append((skip_axis, int(not skip_axis)))
            skip_directions.append((-skip_axis, -int(not skip_axis)))
        elif skip_axis:
            raise ValueError(
                f"{skip_axis} is not a valid value for the axis."
                " You can skip the axis 0 or 1"
            )

        # make a dictionary from positions to index
        index_map = {p: ind for ind, s in enumerate(self) for p in s.positions}
        # create a dictionry from junction direction to strand index
        j_to_ind = {
            d: {index_map[pos] for pos in pos_list}
            for d, pos_list in self.junctions.items()
            if d not in skip_directions
        }

        if until is None:
            until = self.max_pos

        try:
            # this code doens't check for eventual overlaps,
            # so it can lead to MotifStructureError if the junctions are not extendable
            for d, inds in j_to_ind.items():
                for i in inds:
                    s = self[i]
                    if -1 in d:
                        s.extend(direction=d, check=False)
                    # if it's a positive direction and until is None, extend to max_pos
                    if 1 in d:
                        s.extend(direction=d, until=until, check=False)

        # The following code works, it's more efficient, and check for overlaps,
        # but it's hard to maintain, doesn't keep directionality check and understand.
        # So better using the simpler version above.
        # Here the previous code:
        # axis = [d for d in Direction]
        # for axis, sym in ((1, '─'), (0, '│')): # consider the two axis
        #     naxis = int(not axis)
        #     neg_direction = Position((- axis, - naxis))
        #     pos_direction = Position((axis, naxis))

        #     if neg_direction not in skip_directions:
        #         # consider all the junction positions in the negative direction
        #         for pos in juncts[neg_direction]:

        #             # if the position is not at the minimum border
        #             if pos[naxis] != 0:
        #                 # add a strand from the position and going to the border
        #                 extend_strand = Strand(sym * pos[naxis],
        #                                        start=(pos[0] * naxis,
        #                                        pos[1] * axis),
        #                                        direction=pos_direction)

        #                 # check that the extension doesn't overlap
        #                 # with other strands
        #                 if not (set(extend_strand.positions)
        #                                 & set(index_map)):
        #                     # get the strand at the position
        #                     strand_at_pos = self[index_map[pos]]
        #                     strand_at_pos.join(extend_strand)

        #     if pos_direction not in skip_directions:
        #         for pos in juncts[pos_direction]:

        #             strand_at_pos = self[index_map[pos]]
        #             max_pos = self.max_pos[naxis]

        #             # set the maximum position to the until parameter
        #             if until[naxis]:
        #                 max_pos = until[naxis]

        #             if pos[naxis] < max_pos:
        #                 # add a strand from the position and going to the border,
        #                 extend_strand = Strand(sym * (max_pos - pos[naxis]),
        #                                        start=(pos[0] + 1 * axis,
        #                                               pos[1] + 1 * naxis),
        #                                        direction=pos_direction)

        #                 # check that the extension doesn't overlap
        #                 # with other strands
        #                 if not (set(extend_strand.positions)
        #                                 & set(index_map)):
        #                     strand_at_pos.join(extend_strand)

        except MotifStructureError as e:
            print(f"Error in extending junction at strand {s}. Full error:")
            print("\t", e)

        return self

    def flip(
        self,
        horizontally: bool = True,
        vertically: bool = True,
        strand_index: list = None,
        reorder=False,
    ) -> "Motif":
        """
        Flip the strands of the motif horizontally and/or vertically inplace.

        Parameters
        ----------
        horizontally : bool, default True
            If True, flip the strands horizontally.
        vertically : bool, default True
            If True, flip the strands vertically.
        strand_index : list, default None
            The list of indices of the strands to flip.
            If None, all the strands are flipped.
        reorder : bool, default False
            If True, reorder the strands after flipping,
            so that the first strand is the one with the 5' end,
            then the strand with the lowest y start position,
            and finally the strand with the lowest x start position.

        Returns
        -------
        Motif
            The motif with the flipped strands.
        """
        # save the initial index of character and lines,
        # (it changes every time you change a strand otherwise)
        idx_char = self.num_char - 1
        idx_lines = self.num_lines - 1

        # create new basepair dictionary in case autopairing is off
        new_basepair_dict = BasePair()
        for pos1, pos2 in self.basepair.items():
            if horizontally:
                pos1 = (idx_char - pos1[0], pos1[1])  # flip the horizontal position
                pos2 = (idx_char - pos2[0], pos2[1])  # flip the horizontal position
            if vertically:
                pos1 = (pos1[0], idx_lines - pos1[1])  # flip the vertical position
                pos2 = (pos2[0], idx_lines - pos2[1])  # flip the vertical position
            new_basepair_dict[Position(pos1)] = Position(pos2)

        self._basepair = new_basepair_dict  # save the new basepair
        for ind, strand in enumerate(self):
            # flip only the strands in the strand_index list, if it is not None
            if strand_index and ind not in strand_index:
                continue

            # flip the start position of the strands:
            # the new start is the border - the previous position
            new_start = list(strand.start)
            if horizontally:
                new_start[0] = idx_char - strand.start[0]
            if vertically:
                new_start[1] = idx_lines - strand.start[1]

            strand.start = new_start
            strand.flip(horizontally, vertically, flip_start=False)

        if reorder:
            # reorder the strands in order to have:
            # - the strands with 5' end
            # - strand with the lowest y start position
            # - strand with the lowest x start position
            self._strands.sort(key=lambda s: (-int("5" in s.strand), *s.start[::-1]))

        return self

    @wraps(fold_bar)  # inherit the documentation from the function
    def folding_barriers(self, kl_delay: int = 150) -> Tuple[str, int]:
        return fold_bar(self.structure, kl_delay)

    def get_strand_index_map(self) -> Dict[Tuple[int, int], int]:
        """
        Get a dictionary of positions as keys and the index of the
        strand in the motif as values.

        Returns
        -------
        Dict[Tuple[int, int], int]
            A dictionary mapping positions to strand indexes
        """
        return {pos: ind for ind, strand in enumerate(self) for pos in strand.positions}

    def get_position_map(self) -> Dict[Tuple[int, int], str]:
        """
        Get a dictionary of positions as keys and the corresponding
        characters as values.

        Returns
        -------
        Dict[Tuple[int, int], str]
            A dictionary mapping positions to characters.
        """
        return {
            pos: sym
            for strand in self
            for pos, sym in zip(strand.positions, strand.strand)
        }

    def insert(
        self, index: int, strand: "Strand", join: bool = True, copy: bool = False
    ) -> "Motif":
        """
        Insert a strand at a given index in the motif.

        Parameters
        ----------
        index : int
            The index at which to insert the strand.
        strand : Strand
            The strand to be inserted.
        join : bool, default True
            Whether to attempt joining the new strand with existing strands.
        copy : bool, default False
            If True, a copy of the strand is inserted.

        Returns
        -------
        Motif
            The updated motif.

        Raises
        ------
        ValueError
            If `strand` is not a `Strand` instance.
        """
        if not isinstance(strand, Strand):
            raise ValueError(f"{strand} is not a Strand object.")

        if copy:
            strand = strand.copy()

        strand.register_callback(self._updated_strands)
        self._strands.insert(index, strand)

        if join:
            self._strands = self.join_strands(self._strands)

        if self.lock_coords:
            self._strands_block = StrandsBlock(*[s for s in self._strands])

        self._updated_strands()
        return self

    def pop(self, index: int = -1) -> "Strand":
        """
        Remove and return the strand at the specified index.

        Parameters
        ----------
        index : int, default -1
            The index of the strand to remove.

        Returns
        -------
        Strand
            The removed strand.
        """
        strand = self._strands.pop(index)

        # REMOVE THE PAIRS IN WHICH THE STRAND IS INVOLVED
        new_basepair = BasePair()
        for k, v in self._basepair.items():
            if k in strand.positions or v in strand.positions:
                continue
            new_basepair[Position(k)] = Position(v)

        self._basepair = new_basepair
        strand._clear_callbacks()
        self._updated_strands()
        self._updated_basepair()

        if self.lock_coords:
            self._strands_block = StrandsBlock(*[s for s in self._strands])

        return strand

    def replace_all_strands(
        self, strands: List["Strand"], copy: bool = True, join: bool = False
    ) -> "Motif":
        """
        Replace all the strands in the motif with the provided list of strands.

        Parameters
        ----------
        strands : List[Strand]
            The list of strands to replace the current strands.
        copy : bool, default True
            If True, the strands are copied before replacing the current strands.
        join : bool, default False
            Whether to attempt joining the new strands.

        Returns
        -------
        Motif
            The updated motif.

        Raises
        ------
        ValueError
            If any of the provided strands is not a `Strand` instance.
        """
        # copy the strands if needed, registart callbacks
        if copy:
            new_strands = self.copy_strands_preserve_blocks(strands, self)

        else:
            new_strands = strands
            for s in strands:
                if not isinstance(s, Strand):
                    raise ValueError(f"{s} is not a Strand object.")
                s.register_callback(self._updated_strands)

        # join the strands if needed
        if join:
            new_strands = self.join_strands(new_strands)

        self._strands = new_strands
        if self.lock_coords:
            self._strands_block = StrandsBlock(*[s for s in self._strands])

        # this will also update the basepair dictionary
        self._updated_strands()

        return self

    def rotate(self, times: int = 1) -> "Motif":
        """
        Rotate the motif 90 degrees clockwise.

        Parameters
        ----------
        times : int, default 1
            The number of times to rotate the motif.

        Returns
        -------
        Motif
            The rotated motif.
        """
        for _ in range(times % 4):
            # collect the num_lines before is updated
            num_lines = self.num_lines

            for s in self:
                sign = +1
                # when rotating from vertical to horizontal: change sign
                if s._direction[1]:
                    sign = -1
                s._start = s._start.replace(x=num_lines - 1 - s.start[1], y=s.start[0])
                s._direction = s._direction.replace(
                    x=sign * s._direction[1], y=sign * s._direction[0]
                )
                # adjust the sequence symbols
                s.strand = s.strand.translate(rotate_90)

            # rotate the basepair dictionary
            if not self._autopairing and self._basepair:
                new_bp = BasePair()
                for k, v in self._basepair.items():
                    new_k = (num_lines - 1 - k[1], k[0])
                    new_v = (num_lines - 1 - v[1], v[0])
                    new_bp[Position(new_k)] = Position(new_v)
                self._basepair = new_bp

        return self

    def save_3d_model(
        self: Union["Motif", "Origami"],
        filename: str = "motif",
        config: bool = True,
        topology: bool = True,
        forces: bool = False,
        pk_forces: bool = False,
        return_text: bool = False,
        sequence: str = None,
        pdb: bool = False,
        **kwargs,
    ) -> Optional[Tuple[str, str]]:
        """
        Save the motif 3D structure in a file for structural 3D modeling.
        This function save the 3D oxDNA motif representation in conformation and
        topology xoDNA-format files, and optionally in PDB format.
        It can also save additional force constraints for oxDNA simulation
        and protein configurations if specified.

        Parameters
        ----------
        filename : str, optional, default="motif"
            The base filepath for the generated oxDNA files (without extension).
        config : bool, default=True
            If True, generates a configuration (.dat) file with nucleotide positions.
        topology : bool, default=True
            If True, generates a topology (.top) file specifying strand connectivity.
        forces : bool, default=False
            If True, saves force constraints for base-pair interactions.
        pk_forces : bool, default=False
            If True, saves force constraints specifically for pseudoknots.
        return_text : bool, default=False
            If True, returns the generated oxDNA configuration and topology as strings
            instead of writing to files (used for real-time visualization).
        sequence : str, optional
            If provided, uses the given sequence to generate the topology.
        pdb : bool, optional, default=False
            If True, exports a PDB (Protein Data Bank) file for visualization.
        **kwargs
            Additional arguments for customizing force constraints
            and PDB export settings.

        Returns
        -------
        Optional[Tuple[str, str]]
            If `return_text` is True, returns a tuple containing:
            - The configuration file content as a string.
            - The topology file content as a string.
            Otherwise, returns None.

        Notes
        -----
        This function requires the oxDNA analysis tools package to be installed
        for the generation of force constraints and pdb export.
        """

        def get_kwargs_names(func):
            sig = signature(func)
            # Extract parameters that have default values (kwargs)
            kwargs = [
                param.name
                for param in sig.parameters.values()
                if param.default != param.empty
            ]
            return kwargs

        # remove the extension from the filename
        filename = filename.split(".")[0]
        strands = [s for s in self if s.sequence]
        n_nucleotides = sum([len(s.sequence) for s in strands])
        n_strands = len(strands)

        # create the conformation and topology text
        conf_text = "t = 0\nb = 1000 1000 1000\nE = 0 0 0\n"
        topology_text = ""

        ### ADD THE STRANDS TO THE CONFORMATION AND TOPOLOGY TEXTS ###
        for s in strands:
            # check for the sequence direction
            dir = 1 - 2 * (s.directionality == "35")
            seq = str(s.sequence[::dir])
            coord_array = s.coords[::dir]

            # add the coordinates to the conformations text
            if config:
                for pos, a1, a3 in coord_array:
                    conf_text += (
                        f"{pos[0]} {pos[1]} {pos[2]} "
                        f"{a1[0]} {a1[1]} {a1[2]} "
                        f"{a3[0]} {a3[1]} {a3[2]}\n"
                    )

            # add the sequence to the topology text
            if topology:
                topology_text += seq + " type=RNA circular=false \n"

            # add the proteins to the conformation text
            try:
                for protein in s.coords.proteins:
                    if config:
                        for pos, a1, a3 in protein.coords:
                            conf_text += (
                                f"{pos[0]} {pos[1]} {pos[2]} "
                                f"{a1[0]} {a1[1]} {a1[2]} "
                                f"{a3[0]} {a3[1]} {a3[2]}\n"
                            )

                    # add the proteins to the topology text
                    if topology:
                        topology_text += (
                            f"{protein.sequence} " "type=peptide circular=false \n"
                        )

                    n_nucleotides += len(protein)
                    n_strands += 1
            except Exception as e:
                print("Problem with proteins", s, e)

        if sequence:
            topology_text = (
                "\n".join(
                    [f"{seq} type=RNA circular=false" for seq in sequence.split("&")]
                )
                + "\n"
            )

        topology_text = f"{n_nucleotides} {n_strands} 5->3\n" + topology_text

        if return_text:
            return conf_text, topology_text

        # save the files
        conf_file = f"{filename}.dat"
        if config:
            with open(conf_file, "w", encoding="utf-8") as f:
                f.write(conf_text)
        top_file = f"{filename}.top"
        if topology:
            with open(top_file, "w", encoding="utf-8") as f:
                f.write(topology_text)

        ### save the external forces
        if not forces and not pk_forces:
            pass
        elif not oat_installed:
            warnings.warn(
                "oxDNA_analysis_tools is not installed. " "Skipping force writing.",
                UserWarning,
            )
        else:
            trap_kw_names = get_kwargs_names(mutual_trap)
            trap_kwargs = {k: v for k, v in kwargs.items() if k in trap_kw_names}
            pair_map = dot_bracket_to_pair_map(self.structure.replace("&", ""))
            trap_kwargs.setdefault("stiff", 0.09)
            trap_kwargs.setdefault("r0", 1.2)
            trap_kwargs.setdefault("PBC", True)
            trap_kwargs.setdefault("rate", 0)
            trap_kwargs.setdefault("stiff_rate", 0)
            force_list = []
            pk_force_list = []
            ss = self.structure.replace("&", "")
            i = 0
            for n1, n2 in pair_map.items():
                if n2 is None:
                    continue
                i += 1
                trap1 = mutual_trap(n1, n2, **trap_kwargs)
                trap2 = mutual_trap(n2, n1, **trap_kwargs)
                if forces:
                    force_list.append(trap1)
                    force_list.append(trap2)
                if pk_forces and ss[n1] not in ".()":
                    pk_force_list.append(trap1)
                    pk_force_list.append(trap2)
            if forces:
                write_force_file(force_list, f"{filename}_forces.txt")
            if pk_forces:
                write_force_file(pk_force_list, f"{filename}_pk_forces.txt")

        if pdb:
            if not oat_installed:
                warnings.warn(
                    "oxDNA_analysis_tools is not installed. " "Skipping PDB export.",
                    UserWarning,
                )
            else:
                # Read oxDNA configuration
                system, _ = strand_describe(top_file)
                ti, di = describe(top_file, conf_file)
                conf = get_confs(ti, di, 0, 1)[0]
                conf = inbox(conf, center=True)
                # remove the proteins from the configuration if no pdb files provided
                if not kwargs.get("protein_pdb_files"):
                    strand_offset = 0
                    to_pop = []
                    conf_to_keep = []
                    for i, strand in enumerate(system.strands):
                        strand_end = strand_offset + strand.get_length()
                        if strand.type == "peptide":
                            to_pop.append(i)
                        else:
                            conf_to_keep.append((strand_offset, strand_end))
                        strand_offset = strand_end
                    for i in to_pop[::-1]:
                        system.strands.pop(i)

                oxdna_pdb_kw_names = get_kwargs_names(oxDNA_PDB)
                oxdna_pdb_kwargs = {
                    k: v for k, v in kwargs.items() if k in oxdna_pdb_kw_names
                }
                oxdna_pdb_kwargs.setdefault("uniform_residue_names", True)
                oxDNA_PDB(conf, system, filename, **oxdna_pdb_kwargs)

    def save_fasta(self, filename: str = "motif") -> None:
        """
        Save the motif sequences in a FASTA file.

        Parameters
        ----------
        filename : str, optional
            The filepath to save (without extension). Default is 'motif'.
        """
        seqs = self.sequence.split("&")
        dotb = self.structure.split("&")
        with open(f"{filename}.fasta", "w", encoding="utf-8") as f:
            for i, seq in enumerate(seqs):
                f.write(f">strand_{i}\n")
                f.write(f"{seq}\n")
                f.write(f"{dotb[i]}\n")

    def save_text(self, filename_path: str = "motif") -> None:
        """
        Save the motif representation as a text file.

        Parameters
        ----------
        filename_path : str, optional
            The filepath to save (without extension). Default is 'motif'.
        """
        path = filename_path.split(".")[0]
        name = path.split("/")[-1].split("\\")[-1]
        with open(f"{path}.txt", "w", encoding="utf-8") as f:
            f.write(f">{str(name)}\n")
            f.write(f"{self.sequence}\n")
            f.write(f"{self.structure}\n\n")
            f.write(str(self))

    def shift(self, shift_vect: Tuple[int, int], extend: bool = False) -> "Motif":
        """
        Shift the motif of the given shift vector.

        Parameters
        ----------
        shift_vect : Tuple[int, int]
            The (x, y) shift values.
        extend : bool, default False
            Whether to extend junctions when shifting
            (in the direction opposite to the shifting direction).

        Returns
        -------
        Motif
            The shifted motif.

        Raises
        ------
        ValueError
            If shifting moves strands to negative positions.
        """
        Strand._check_position(input_pos_dir=shift_vect)
        shift = Position(shift_vect)
        min_pos = self.min_pos

        # check if the shift will bring the strands to negative positions
        if min_pos[0] + shift[0] < 0 or min_pos[1] + shift[1] < 0:
            raise MotifStructureError(
                f"The motif cannot be shifted. The strands cannot be drawn"
                f" at negative positions. Attempt to draw the motif at "
                f"position ({min_pos[0] + shift_vect[0]}, "
                f"{min_pos[1] + shift_vect[1]})"
            )

        # shift the strands
        for s in self:

            # save the initial start and end positions
            s_start = s.start
            s_end = s.end

            # shift the strand
            s.shift(shift, check=False)

            if extend:
                # extend the strand in the opposite direction
                s.extend(direction=s.end_direction, until=s_end, check=False)
                s.extend(direction=-s.direction, until=s_start, check=False)

        # shift every basepair too
        if not self._autopairing:
            self._basepair = self.basepair.shift(shift)
        return self

    def sort(
        self, key: Optional[Callable[["Strand"], Any]] = None, reverse: bool = False
    ) -> "Motif":
        """
        Sort the strands in the motif.

        Parameters
        ----------
        key : function, optional
            The function to use to sort the strands. If ``None``, the strands are
            sorted by:

            - strands with a 5' end first
            - lowest starting y position
            - lowest starting x position

        reverse : bool, default False
            If ``True``, sort in descending order.

        Returns
        -------
        Motif
            The motif with the sorted strands.
        """
        if not key:
            # sort the strand according to the lowest start position
            def key(s):
                return (-int("5" in s.strand), *s.start[::-1])

        self._strands.sort(key=key, reverse=reverse)
        return self

    def strip(self, skip_axis: Literal[None, 0, 1] = None) -> "Motif":
        """
        Remove the empty lines/columns in the motif structure.

        Parameters
        ----------
        skip_axis: int, optional
            The numpy axis to skip when stripping the motif
            (1 for horizontal, 0 for vertical).

        Returns
        -------
        Motif
            The stripped motif.

        """
        min_pos = self.min_pos
        shift = (-min_pos[0] * int(skip_axis != 1), -min_pos[1] * int(skip_axis != 0))
        self.shift(shift)
        return self
