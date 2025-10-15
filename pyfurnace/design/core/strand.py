import warnings
import copy
from typing import Optional, Union, Literal, List, Dict, Tuple, Any, Set
import numpy as np

try:
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

from .symbols import (
    nucleotides,
    symb_to_road,
    nucl_to_none,
    symb_to_none,
    accept_symbol,
    horiz_flip,
    verti_flip,
    MotifStructureError,
    AmbiguosStructure,
)
from .callback import Callback
from .position import Position, Direction
from .sequence import Sequence
from .coordinates_3d import Coords


class Strand(Callback):
    """
    Represents a single RNA strand in a hybrid 1D/2D/3D representation.

    Parameters
    ----------
    strand : str, optional
        The strand text characters (default is empty).
    directionality : {'53', '35'}, optional
        The directionality of the strand (default is '53').
    start : Position, optional
        The starting position of the strand in 2D space (default is Position.zero()).
        If set to 'minimal', the start position is calculated as the minimal
        start position required to draw the strand in 2D to avoid negative coordinates.
    direction : Direction, optional
        The initial direction of the strand in 2D space (default is Direction.RIGHT).
    coords : Coords, optional
        The 3D coordinates for the strand (if not provided, they will be calculated).
    strands_block : StrandsBlock, optional
        The block of strands to which the strand belongs.
    **kwargs : dict
        Arbitrary keyword arguments passed to the Callback class.

    Attributes
    ----------
    coords : Coords
        The 3D coordinates for the strand.
    direction : Direction
        The direction of the strand in 2D space.
    directionality : str
        The directionality of the sequence ('53' or '35').
    directions : Tuple[Direction]
        The directions of each strand character in 2D space (x,y coordinates).
    end : Position
        The last position of the strand in 2D space.
    end_direction : Direction
        The last direction of the strand in 2D space.
    max_pos : Position
        The maximum x, y coordinates of the strand in 2D space.
    min_pos : Position
        The minimum x, y coordinates of the strand in 2D space.
    minimal_dimensions : Position
        The minimum start Positions required to draw the strand in 2D
        without reaching negative positions.
        Calculating the minimal dimensions doesn't take into
        acount structural errors.
    next_pos : Position
        The next position of the strand in 2D space.
    pk_info : Optional[dict]
        The pseudoknot information of the strand (if any).
    positions : Tuple[Position]
        The positions of each strand character in 2D space (x,y coordinates).
    prec_pos : Position
        The preceding position of the strand in 2D space.
    seq_positions : Tuple[Position]
        The positions of each nucleotide in the strand sequence (x,y coordinates).
    sequence : Sequence
        The nucleotide sequence of the strand.
    sequence_list : List[Sequence]
        The list of consecutive sequences in the strand.
    sequence_slice : List[slice]
        The list of slices of the sequence in the strand.
    start : Position
        The starting position of the strand.
    strand : str
        The strand string representation.
    strands_block : StrandsBlock
        The block of strands to which the strand belongs.
    """

    def __init__(
        self,
        strand: str = "",
        directionality: Literal["53", "35"] = "53",
        start: Optional[Union[Position, tuple, Literal["minimal"]]] = None,
        direction: Optional[Direction] = None,
        coords: Coords = None,
        strands_block: Optional["StrandsBlock"] = None,
        **kwargs,
    ) -> None:
        """
        Initialize the Strand object with given properties.

        Parameters
        ----------
        strand : str, optional
            The characters strand representation (default is empty).
        directionality : {'53', '35'},
            The directionality of the strand (default is '53').
        start : Position, default Position.zero()
            The starting position of the strand in 2D space.
            If start is set to 'minimal', the start position is calculated
            as the minimal start position required to draw the strand in 2D
            to avoid negative coordinates.
        direction : Direction, default is Direction.RIGHT
            The initial direction of the strand in 2D space.
        coords : Coords, optional
            The 3D coordinates for the strand.
        strands_block : StrandsBlock, optional
            The block of strands to which the strand belongs.
        **kwargs : dict
            Arbitrary keyword arguments passed to the Callback class.
        """
        # register the callback
        super().__init__(**kwargs)
        self.pk_info = None
        self._coords = Coords()
        self._positions = None
        self._directions = None
        self._seq_positions = None
        self._prec_pos = None
        self._next_pos = None
        self._max_pos = None
        self._min_pos = None
        self._seq_slice = []

        strand = self._check_line(strand)

        # Set the strand and sequence
        self._update_sequence_insertion(strand, directionality)

        ### check the 2D direction
        if direction is None:
            self._direction = Direction.RIGHT
        else:
            self._direction = self._check_position(direction, True)

        ### check the start position
        if start is None:
            self._start = Position.zero()
        elif start == "minimal":
            self._start = self.minimal_dimensions
        else:
            self._start = self._check_position(start)

        ### set 2D properties
        self._reset_positions()

        if coords is not None:
            self.coords = coords

        # set the strand block to lock the 3D coordinates
        if strands_block is not None:
            self.strands_block = strands_block
        else:
            self.strands_block = StrandsBlock(self)

    def __str__(self) -> str:
        """Return a string representation of the Strand object."""
        return self.strand

    def __repr__(self) -> str:
        """Return a string representation of the Strand object."""
        return self.strand

    def __getitem__(self, idx) -> str:
        """Return the character at the given index."""
        return self.strand[idx]

    def __setitem__(self, idx: Union[int, slice], val: str) -> None:
        """Set the character at the given index in the strand."""
        val = val.upper()
        # convert to list to accept slice indexing and negative indexing
        new_strand_line = list(self.strand)
        new_strand_line[idx] = val
        new_strand_str = "".join(new_strand_line)
        ### update the strand
        self.strand = new_strand_str

    def __add__(self, other: Union[str, Sequence, "Strand"]) -> "Strand":
        """Add two strands together, check the directionality."""
        # check and copy the other strand
        other = self._check_addition(other)
        # combine the 3D coordinates
        coords = Coords.combine_coords(self, other)
        # create a new strand combined strand
        new_strand = Strand(
            str(self) + str(other),
            self.directionality,
            self.start,
            self.direction,
            callbacks=self._callbacks,
            coords=coords,
        )

        # combine the pseudoknots information and add them
        new_strand.pk_info = self._combine_pk_info(other)

        return new_strand

    def __iadd__(self, other: Union[str, Sequence, "Strand"]) -> "Strand":
        """In place addition of two strands together."""
        # check and copy the other strand
        other = self._check_addition(other)
        # add the strand to the current strand
        self.strand = self.strand + other.strand
        # combine the 3D coordinates
        self._coords = Coords.combine_coords(self, other)

        # combine the pseudoknots information and add them
        self.pk_info = self._combine_pk_info(other)

        # notify the upper class that the strand has changed
        self._trigger_callbacks()
        return self

    def __radd__(self, other: Union[str, Sequence]) -> "Strand":
        """Right addition a Strand to a string or a Sequence."""
        if other == 0:
            return self
        self._check_addition(other, copy=False)
        # the only accepted options are str and Sequence
        return Strand(
            str(other) + str(self),
            directionality=self.directionality,
            start=self.start,
            direction=self.direction,
            callbacks=self._callbacks,
        )

    def __len__(self) -> int:
        """Return the length of the strand."""
        return len(self.strand)

    def __eq__(self, other) -> bool:
        """Check if two strands are equal, or if the strand is equal to a
        string or a Sequence."""
        if isinstance(other, Strand):
            return (
                self.positions == other.positions
                and self.directionality == other.directionality
            )
        elif isinstance(other, Sequence):
            return self.sequence == other
        elif isinstance(other, str):
            return self.strand == other
        return False

    def __bool__(self):
        """Return true if the strand has strand characters."""
        return bool(self.strand)

    def __contains__(self, other):
        """Check if a position or a substrand is included in the strand."""
        # check a position (a tuple or list containing two integers)
        if isinstance(other, (Position, tuple, list)):
            return other in self.positions
        if isinstance(other, Sequence):
            return other in self.sequence
        if isinstance(other, str):
            return other in str(self.strand)
        # check a substrand
        if isinstance(other, Strand):
            # zip the strand and the positions
            sub = list(zip(other.strand, other.positions))
            full = list(zip(self.strand, self.positions))

            # sliding window to check if the sub is in the full
            sub_len = len(sub)
            for i in range(len(full) - sub_len + 1):
                if full[i : i + sub_len] == sub:
                    return True

            # if not found, return False
            return False

        return False

    def __hash__(self):
        return id(self)

    ###
    ### PROPERTIES
    ###

    @property
    def coords(self) -> Coords:
        """
        The oxDNA 3D coordinates of the strand.
        If the coordinates are not set, they are calculated for a strand in
        a perfect helix with the same length of the strand.
        """
        if not self._sequence or not self._coords.is_empty():
            return self._coords
        coords = Coords.compute_helix_from_nucl(
            (0, 0, 0),  # start position
            (1, 0, 0),  # base vector
            (0, 1, 0),  # normal vector
            length=len(self.sequence),
            directionality=self.directionality,
            double=False,
        )
        self._coords = coords
        return coords

    @coords.setter
    def coords(self, new_coords: Union[Coords, str, list, tuple] = None) -> None:
        """
        Set the oxDNA 3D coordinates of the strand.
        Check the Coords documentation for the format of the coordinates.
        """
        if isinstance(new_coords, np.ndarray):
            if new_coords.size == 0:
                self._coords = Coords()
                return
        elif not new_coords:
            self._coords = Coords()
            return

        # sanity check
        if len(new_coords) != len(self.sequence):
            raise MotifStructureError(
                f"The number of oxDNA coordinates "
                f"({len(new_coords)}) is different "
                f"from the number of nucleotides "
                f"({len(self.sequence)})"
            )

        if not isinstance(new_coords, Coords):
            new_coords = Coords(new_coords)

        self._coords = new_coords

    @property
    def direction(self) -> Direction:
        """
        The starting (x, y) 2D direction to build the strand in 2D.
        """
        return self._direction

    @direction.setter
    def direction(self, new_direction: Union[Direction, Position, tuple, list]) -> None:
        """
        The starting direction setter: a tuple of x,y coordinates in which one
        coordinate is 0 and the other either 1 or -1.
        """
        self._direction = self._check_position(new_direction, True)
        self._reset_positions()
        self._trigger_callbacks()

    @property
    def directionality(self) -> Literal["53", "35"]:
        """
        The directionality of the sequence ('53' or '35')
        """
        return self.sequence.directionality

    @directionality.setter
    def directionality(self, new_directionality: Literal["53", "35"]) -> None:
        """
        Set the directionality of the sequence ('53' or '35').
        """
        self.sequence._directionality = new_directionality
        self._trigger_callbacks()

    @property
    def directions(self) -> Tuple[Direction]:
        """
        A tuple of the directions of each character of the strand in 2D space
        (x,y coordinates).
        """
        if not self._directions:
            self._calculate_positions()
        return self._directions

    @property
    def end(self) -> Position:
        """
        The last position of the strand: it is the last element of the positions.
        """
        if not self._positions:
            self._calculate_positions()
        return self._positions[-1]

    @property
    def end_direction(self) -> Direction:
        """
        The last direction of the strand.
        """
        if not self._positions:
            self._calculate_positions()
        return self._directions[-1]

    @property
    def max_pos(self) -> Position:
        """
        The maximum x, y coordinates of the strand in 2D space.
        """
        if not self._positions:
            self._calculate_positions()
        return self._max_pos

    @property
    def minimal_dimensions(self) -> Position:
        """
        The minimum start Positions required to draw the strand in 2D
        without reaching negative positions.
        Calculating the minimal dimensions doesn't take into
        acount structural errors.
        """
        # initialize the start position
        start = Position.zero()
        pos = Position.zero()
        direction = self.direction

        for sym in self:
            # check that the structure does not reach negative coordinates,
            # in case add 1 to all negative positions
            negative = Position((int(x < 0) for x in pos))

            # add +1 to the element that is negative
            start += negative
            pos += negative

            ### UPDATE DIRECTIONS ###
            if sym in "╰╮\\":
                direction = direction.swap_xy()
            elif sym in "╭╯/":
                direction = direction.swap_change_sign_xy()
            elif sym == "⊗":
                direction += Position((0, 0, 1))
            elif sym == "⊙":
                direction -= Position((0, 0, 1))
            # else the direction is the same

            # calculate new position
            pos += direction

        return start

    @property
    def min_pos(self) -> Position:
        """
        The minimum x, y coordinates of the strand in 2D space.
        """
        if not self._positions:
            self._calculate_positions()
        return self._min_pos

    @property
    def next_pos(self) -> Position:
        """
        The next position of the strand: it is the last position plus the
        end_direction. Useful for joining two strands.
        """
        if not self._positions:
            self._calculate_positions()
        return self._next_pos

    @property
    def positions(self) -> Tuple[Position]:
        """
        A tuple of the positions of each character of the strand in 2D space
        (x,y coordinates).
        """
        if not self._positions:
            self._calculate_positions()
        return self._positions

    @property
    def prec_pos(self) -> Position:
        """
        The preceding position of the strand: it is the start position
        minus the direction. Useful for joining two strands.
        """
        # i could calculate this as self.start - self.direction
        # but if there are no positions, they will be likely needed soon
        if not self._positions:
            self._calculate_positions()
        return self._prec_pos

    @property
    def sequence(self) -> Sequence:
        """
        The complete sequence of the strand.
        """
        return self._sequence

    @sequence.setter
    def sequence(self, new_seq: Union[str, Sequence]) -> None:
        """
        Set the sequence of the strand.
        If the sequence is a string, it is split into the sequences of the strand.
        If the new sequence is longer than the current sequence, the sequence is
        added at the end of the strand.
        The strands keeps the current directionality, indipendently from the
        directionality of the new sequence.
        """
        new_seq = self._check_line(new_seq)
        self._updated_sequence(new_seq)

    @property
    def sequence_list(self) -> List[Sequence]:
        """
        The list of consecutive sequences in the strand.
        """
        seq_list = []
        for sl in self._seq_slice:  # iterate over the slices of the sequence
            seq_list.append(Sequence(self.strand[sl], self.directionality))
        return seq_list

    @property
    def sequence_slice(self) -> List[slice]:
        """
        The list of slices of the sequence in the strand.
        """
        return self._seq_slice

    @property
    def seq_positions(self) -> Tuple[Position]:
        """
        The positions of each nucleotide in the strand sequence (x,y coordinates).
        """
        if not self._positions:
            self._calculate_positions()
        return self._seq_positions

    @property
    def strand(self) -> str:
        """
        The strand text characters.
        """
        return self._strand

    @strand.setter
    def strand(self, new_strand: str) -> None:
        """
        Set the strand text characters.
        """
        # check the strand
        new_strand = self._check_line(new_strand)
        # update the strand and the sequence
        self._update_sequence_insertion(new_strand)
        # update the positions
        self._reset_positions()
        # notify the upper class that the strand has changed
        self._trigger_callbacks()

    @property
    def start(self):
        """
        Start position of the strand in the 2D representation.
        """
        return self._start

    @start.setter
    def start(self, new_start: Union[Position, tuple, list]) -> None:
        """
        Start position setter: check that the input is
        a Position/tuple of x,y coordinates.
        """
        self._start = self._check_position(new_start)
        # reset the 2D maps
        self._reset_positions()
        self._trigger_callbacks()

    @property
    def strands_block(self) -> "StrandsBlock":
        """
        The block of strands that has locked 3D coordinates
        transformation (their relative position in the 3D space
        doesn't change when joined to other strands).
        """
        return self._strands_block

    @strands_block.setter
    def strands_block(self, new_block: Optional["StrandsBlock"] = None) -> None:
        """
        Assign a new strands block to the strand.
        """
        # unlock this strand from the previous block
        if new_block is None:
            self._strands_block = StrandsBlock(self)
            return

        # check if the new block is a StrandsBlock
        if not isinstance(new_block, StrandsBlock):
            raise ValueError(
                f"The strands block must be a StrandsBlock object. "
                f"Got {type(new_block)} instead."
            )

        # assign the new block to the strand
        self._strands_block = new_block
        self._strands_block.add(self)

    ###
    ### STATIC METHODS
    ###

    @staticmethod
    def join_strands(strand1: "Strand", strand2: "Strand") -> Optional["Strand"]:
        """
        Join two strands together if they match in the right direction and have
        a compatible directionality. This also combines the pseudoknot information
        and the 3D coordinates of the strands (transforming the coordinates of the
        second strand to match the first strand extension). The directionality of
        the combined strand is taken from the first strand if it has a sequence,
        otherwise it is taken from the second strand. The callbacks of the first
        strand are used for the new strand. The strands block of the joined strand
        includes the strands blocks of both strands.

        Parameters
        ----------
        strand1 : Strand
            The first strand to join.
        strand2 : Strand
            The second strand to join.

        Returns
        -------
        Optional[Strand]
            The joined strand if possible, otherwise None.
        """
        if not isinstance(strand2, Strand) or not isinstance(strand1, Strand):
            raise TypeError(
                f"The objects to join are not a Strand object, "
                f"got {type(strand2)}, {type(strand2)} instead."
            )

        ### LOAD THE STRAND EDGES
        s1_start_prev = strand1.prec_pos
        s1_end_next = strand1.next_pos
        # check also strand2 or we could have problems if the strand has one symbol
        s2_start_prev = strand2.prec_pos
        s2_end_next = strand2.next_pos

        ### There are 4 case:
        #       in 2 case we have to invert the second strand
        #       in 2 case we can just add them

        ### CHECK THE SECOND STRAND INVERSION
        s2_inverted = False

        # case 1: (1)-->,<--(2)
        if s1_end_next == strand2.end and s2_end_next == strand1.end:
            strand2.invert()
            s2_inverted = True

        # case 2: <--(1),-->(2)
        elif s1_start_prev == strand2.start and s2_start_prev == strand1.start:
            strand2.invert()
            s2_inverted = True

        ### CHECK ADDITION AND JOINING
        join_order = None
        joined = None

        # case 3: (1)-->,-->(2)
        if s1_end_next == strand2.start:
            join_order = 1

        # case 4: <--(1),<--(2)
        elif s1_start_prev == strand2.end:
            join_order = -1

        # join them in the order found
        if join_order is not None:
            strand1._check_addition(strand2, copy=False)

            if join_order == 1:
                first, second = strand1, strand2
            else:
                first, second = strand2, strand1

            # join the coordinates
            coords = Coords.combine_coords(first, second)

            # set the directionality of the strand with the sequence
            # always priority to strand1
            directionality = strand1.directionality
            if not strand1._sequence:
                directionality = strand2.directionality

            # create a new joined strand
            joined = Strand(
                first.strand + second.strand,
                directionality=directionality,
                start=first.start,
                direction=first.direction,
                coords=coords,
                # always priority to strand1 callbacks
                callbacks=strand1.callbacks,
            )

            # update the positional parameters
            if strand1._positions and strand2._positions:
                joined._positions = first._positions + second._positions
                joined._directions = first._directions + second._directions
                joined._seq_positions = first._seq_positions + second._seq_positions
                joined._prec_pos = first._prec_pos
                joined._next_pos = second._next_pos
                joined._max_pos = Position(
                    (
                        max(first._max_pos[0], second._max_pos[0]),
                        max(first._max_pos[1], second._max_pos[1]),
                    )
                )
                joined._min_pos = Position(
                    (
                        min(first._min_pos[0], second._min_pos[0]),
                        min(first._min_pos[1], second._min_pos[1]),
                    )
                )

            # update the pseudoknots information
            joined.pk_info = first._combine_pk_info(second)

        # revert back the strand if it was inverted
        if s2_inverted:
            strand2.invert()

        # check if the strand is empty
        if joined is None:
            return

        # join the strands block of the combined and first strand
        joined.strands_block.update(
            strand1.strands_block,
            strand2.strands_block,
            avoid_strands={strand1, strand2},
        )
        return joined

    @staticmethod
    def _check_position(input_pos_dir, direction: bool = False) -> Position:
        """
        Check that the input is a valid position or direction.

        Parameters
        ----------
        input_pos_dir : Union[Position, tuple, list]
            The input position or direction to check.
        direction : bool, default False
            If True, check if the input is a direction (a tuple with one coordinate
            equal to 0 and the other either 1 or -1).

        Returns
        -------
        Position
            The input position or direction as a Position object.
        """
        if direction and len(set(input_pos_dir) & {-1, 0, 1}) != 2:
            raise ValueError(
                f"2D direction not allowed. The allowed values are:\n"
                f"RIGHT (1, 0); DOWN (0, 1); LEFT (-1, 0); UP (0, -1).\n"
                f"Got {input_pos_dir} instead."
            )

        if isinstance(input_pos_dir, (Direction, Position)):
            return input_pos_dir

        if (
            not isinstance(input_pos_dir, (tuple, list))
            or not isinstance(input_pos_dir[0], (int, np.int64))
            or not isinstance(input_pos_dir[1], (int, np.int64))
        ):
            raise ValueError(
                f"The 2D coordinates must be a tuple/list of (x,y) "
                f"integer values. Got {input_pos_dir} instead."
            )

        return Position(input_pos_dir)

    ###
    ### PROTECTED METHODS
    ###

    def _calculate_positions(self) -> None:
        """
        Trace the strand in 2D space and create the 2D properties of the strand.
        It cleans the symbols, replacing them with the corresponding nice ASCII
        representation. It also checks for structural errors in the strand
        drawing.
        This function calculates the properties:
            - _prec_pos
            - _positions
            - _directions
            - _next_pos
            - _max_pos
            - _min_pos

        Raises
        ------
        MotifStructureError
            If the strand has structural errors in the 2D representation.
            (e.g. negative coordinates, wrong symbols, etc.)
        """
        ### initialize all the variables
        pos = self.start
        max_pos = list(pos)
        min_pos = list(pos)
        direction = self.direction
        self._prec_pos = pos - direction
        positions = []
        # pos_set = set() # to check weird crossings
        directions = [direction]
        seq_positions = []
        # Use list for efficient string concatenation
        new_strand = []
        # first conversion of the symbols that are easy to interpret
        translated = self.strand.translate(symb_to_road)

        # initialize the error message
        xy_error_msg = (
            "The {xy_axis} direction of the strand (d{xy_axis}: {dir_xy})"
            " is not compatible with the next symbol ({cur_sym}).\n"
            "\tCurrent strand: {cur_strand}.\n"
            "\tFull strand: {full_strand}."
        )

        # traverse the strand building the 2D map
        for i, sym in enumerate(translated):

            # Convert the turns to the corresponding symbols
            if sym == "/":
                if direction[0] == -1 or direction[1] == -1:
                    sym = "╭"
                elif direction[0] == 1 or direction[1] == 1:
                    sym = "╯"
            elif sym == "\\":
                if direction[0] == 1 or direction[1] == -1:
                    sym = "╮"
                elif direction[0] == -1 or direction[1] == 1:
                    sym = "╰"
            # conver the up and down symbols to match the directionality
            elif sym in ("↑", "↓"):
                if (
                    direction[1] == -1
                    and self.directionality == "53"
                    or direction[1] == 1
                    and self.directionality == "35"
                ):
                    sym = "↑"
                else:
                    sym = "↓"
            elif sym in nucleotides:
                # add the base positions
                seq_positions.append(pos)

            # Append the new symbol to the strand list
            new_strand.append(sym)

            # Check for invalid positions and raise error early
            if any(c < 0 for c in pos):
                raise MotifStructureError(
                    f"The strand reaches negative coordinates: "
                    f"{pos}. The current positions are: "
                    f"{positions}."
                    f" The last direction is: {direction}. "
                    f"Current strand: {''.join(new_strand)}"
                )

            # Check symbols and directions errors
            elif direction[0] and any(
                (
                    sym == "│",
                    direction[0] == 1 and sym in "╰╭",
                    direction[0] == -1 and sym in "╮╯",
                )
            ):
                cur_strand = "".join(new_strand)
                raise MotifStructureError(
                    xy_error_msg.format(
                        xy_axis="x",
                        dir_xy=direction[0],
                        cur_sym=sym,
                        cur_strand=cur_strand,
                        full_strand=self.strand,
                    )
                )

            elif direction[1] and any(
                (
                    sym == "─",
                    direction[1] == 1 and sym in "╭╮",
                    direction[1] == -1 and sym in "╰╯",
                )
            ):
                cur_strand = "".join(new_strand)
                raise MotifStructureError(
                    xy_error_msg.format(
                        xy_axis="y",
                        dir_xy=direction[1],
                        cur_sym=sym,
                        cur_strand=cur_strand,
                        full_strand=self.strand,
                    )
                )

            ### DON'T NAME THE CROSSINGS IF THEY WORKS
            # allowed crossing, signal it
            elif pos in positions and sym not in "┼":
                warnings.warn(
                    f"The strand is doing a crossing not allowed, "
                    f"'{sym}' is trying to overwrite the symbol at position"
                    f"'{pos}'. The symbol will be overwritten "
                    f"with '┼'.",
                    AmbiguosStructure,
                    stacklevel=3,
                )
                sym = "┼"
                new_strand[-1] = sym

            # Add symbol to the 2D map
            positions.append(pos)

            if pos[0] > max_pos[0]:
                max_pos[0] = pos[0]
            if pos[1] > max_pos[1]:
                max_pos[1] = pos[1]
            if pos[0] < min_pos[0]:
                min_pos[0] = pos[0]
            if pos[1] < min_pos[1]:
                min_pos[1] = pos[1]

            # Update directions
            if sym in "╰╮\\":
                direction = direction.swap_xy()
            elif sym in "╭╯/":
                direction = direction.swap_change_sign_xy()
            # EXPERIMENTAL
            elif sym == "⊗":
                direction = direction + Position((0, 0, 1))
            elif sym == "⊙":
                direction = direction - Position((0, 0, 1))

            directions.append(direction)

            # Update positions
            pos = pos + direction

        # Update the strand symbols
        new_strand = "".join(new_strand)
        if self._strand != new_strand:
            self._strand = new_strand

        # Update the 2D properties
        self._positions = tuple(positions)
        self._directions = tuple(directions)
        self._seq_positions = tuple(seq_positions)
        self._next_pos = pos
        self._max_pos = Position(max_pos)
        self._min_pos = Position(min_pos)

    def _check_addition(
        self, other: Union[str, Sequence, "Strand"], copy: bool = True
    ) -> "Strand":
        """
        Check a second object for addition to the strand.
        If the object is a string or Sequence, it is converted to a Strand object.

        Parameters
        ----------
        other : Union[str, Sequence, Strand]
            The object to check for addition.
        copy : bool, default True
            If True, return a copy of the object.
            If False, return the object itself.

        Returns
        -------
        Strand
            The object to add to the strand.

        Raises
        -------
        TypeError
            If the object is not a valid type for addition.
        MotifStructureError
            If the strand has structural errors in the 2D representation.
            (e.g. different directionality, etc.)
        """
        # + str
        if isinstance(other, str):
            return Strand(other, self.directionality)

        # + Sequence error
        elif (
            isinstance(other, Sequence)
            and (self.sequence or "5" in self._strand or "3" in self._strand)
            and self.directionality != other.directionality
        ):
            raise MotifStructureError(
                f"Cannot add a strand with a Sequence with "
                f"different directionality.\n"
                f"\tStrand: {self._strand}. \n"
                f"\tStrand directionality: "
                f"{self.directionality}\n"
                f"\tSequence: {other}\n"
                f"\tSequence directionality: "
                f"{other.directionality}"
            )
        # + Sequence no error
        elif isinstance(other, Sequence):
            return Strand(other, other.directionality)

        # + Not a strand
        elif not isinstance(other, Strand):
            raise TypeError(f"{other} is not a valid type for addition")

        # + Strand error
        if (
            (self.sequence or "5" in self._strand or "3" in self._strand)
            and (other.sequence or "5" in other._strand or "3" in other._strand)
            and self.directionality != other.directionality
        ):
            raise MotifStructureError(
                f"Cannot add two strands with different "
                f"directionality. \n"
                f"\tFirst strand: {self._strand} \n"
                f"\tFirst strand directionality: "
                f"{self.directionality}\n"
                f"\tSecond strand {other._strand}\n"
                f"\tSecond strand directionality: "
                f"{other.directionality}"
            )

        if copy:
            return other.copy()

        return other

    def _check_line(
        self, line: Union[str, Sequence], translator: Dict = symb_to_none
    ) -> str:
        """
        Check that the input is a valid strand or sequence.
        The translator is used to check the symbols in the line
        are valid for the strand (translator: symb_to_none) or
        valid for the sequence (translator: nucl_to_none).

        Parameters
        ----------
        line : Union[str, Sequence]
            The input strand to check.
        translator : dict, default symb_to_none
            The translator dictionary to check the strand symbols.

        Raises
        ------
        TypeError
            If the input is not a valid strand or sequence.
        MotifStructureError
            If the strand has structural errors in the 2D representation.
            (e.g. terminal symbols in the middle of the strand, etc.)
        ValueError
            If the input contains invalid symbols for the strand.
        """
        # all good for Sequences
        if isinstance(line, Sequence):
            return line._sequence

        # Not a string
        if not isinstance(line, str):
            raise TypeError(
                f"The strand must be a string or a sequence object. "
                f"Got {type(line)} instead."
            )

        line = line.upper()  # convert to uppercase
        # Check that the terminal symbols are not in the middle of the strand
        for term_sym in ("3", "5"):
            if term_sym not in line:
                continue
            if line.count(term_sym) > 1:
                raise MotifStructureError(
                    f"The strand can have only one start and one "
                    f"end symbol ('5' and '3'). "
                    f"Got {line.count(term_sym)} "
                    f"'{term_sym}' symbols."
                )

            if term_sym not in (line[0], line[-1]):
                raise MotifStructureError(
                    f"The '5' and '3' symbols are terminal "
                    f"symbols. Got '{term_sym}' end at index "
                    f"'{line.index(term_sym)}' instead.\n"
                    f"\tFull strand: {line}."
                )

        allowed_symbols = accept_symbol
        if translator == nucl_to_none:
            allowed_symbols = nucleotides

        if line.translate(translator):
            raise ValueError(
                f"The string '{line}' contains symbols not allowed in "
                f"ROAD. The symbols allowed are: {allowed_symbols}."
            )

        return line

    def _combine_pk_info(self, other: "Strand") -> Optional[Dict]:
        """
        Combine the pseudoknot information of two strands.
        If the strands have no pseudoknot information, return None.
        The pk_info is a dictionary with:
            - id: list of pseudoknot ids (with hyphen for complementarity)
            - ind_fwd: list of tuples with the first and last index
                        of the pseudoknot in the sequence
            - E: list of energies of the pseudoknots
            - dE: list of dE values of the pseudoknots

        Parameters
        ----------
        other : Strand
            The other strand to combine with.

        Returns
        -------
        Optional[Dict]
            The combined pseudoknot information if available, otherwise None.
        """
        # check if the strands have pseudoknot information
        if not self.pk_info and not other.pk_info:
            return None

        new_pk_info = {"id": [], "ind_fwd": [], "E": [], "dE": []}
        if self.pk_info:
            new_pk_info = copy.deepcopy(self.pk_info)
        if other.pk_info:
            offset = len(self.sequence)
            new_pk_info["id"] += other.pk_info["id"]
            new_pk_info["ind_fwd"] += [
                (x[0] + offset, x[1] + offset) for x in other.pk_info["ind_fwd"]
            ]
            new_pk_info["E"] += other.pk_info["E"]
            new_pk_info["dE"] += other.pk_info["dE"]
        return new_pk_info

    def _reset_positions(self) -> None:
        """
        Reset the strand 2D properties: positions, directions and base positions.
        """
        # Main positional parameters
        self._positions = None
        self._seq_positions = None
        self._directions = None

    def _updated_sequence(
        self, new_sequence: Optional[Union[str, Sequence]] = None, **kwargs
    ) -> None:
        """
        Callback function to update the sequence of the strand.
        It updates the sequence properties and the strand string.

        Parameters
        ----------
        new_sequence : str or Sequence, default None
            The new sequence to update.
        kwargs : dict
            Arbitrary keyword arguments passed to the Callback class.
        """
        # cerate a new sequence object and register the callback
        if isinstance(new_sequence, str):
            new_sequence = Sequence(
                new_sequence,
                self._sequence.directionality,
                callback=self._updated_sequence,
            )

        new_seq_str = new_sequence._sequence

        # reset the coordinates and positions if the sequence changed length
        if len(new_seq_str) != len(self._sequence):
            self._coords = Coords(())  # reset the oxDNA coordinates
            self._reset_positions()

        # assign the new sequence
        self._sequence = new_sequence
        build_strand = self._strand.translate(nucl_to_none)
        seq_ind = 0
        # iterate over the slices of the sequence
        for sl in self._seq_slice:
            # the length of the sequence slice
            subseq_len = sl.stop - sl.start
            # replace the sequence slice with the new sequence
            build_strand = (
                build_strand[: sl.start]
                + new_seq_str[seq_ind : seq_ind + subseq_len]
                + build_strand[sl.start :]
            )
            seq_ind += subseq_len

        # add the rest of the sequence if any
        if new_seq_str[seq_ind:]:
            build_strand += new_seq_str[seq_ind:]

        # update the strand
        self._strand = build_strand
        self._trigger_callbacks(**kwargs)

    def _update_sequence_insertion(
        self, strand: str, directionality: str = None
    ) -> None:
        """
        This function takes a strands string and updates the properties of the strand,
        sequence and sequence slice.

        Parameters
        ----------
        strand : str
            The strand string to update.
        directionality : str, default None
            The directionality of the strand. If None, it is taken from the
            current sequence.
        """
        # initialize the variables
        strand_str = str(strand)
        new_sequence = ""
        current_sequence = ""
        seq_slice = []

        # iterate over the strand
        for ind, sym in enumerate(strand_str + " "):
            # if the symbol is a nucleotide, add it to the current sequence
            if sym in nucleotides:
                current_sequence += sym

            else:
                # if the current sequence is not empty
                # and the current symbol is not a nucleotide,
                # add it to the sequence slice and reset the current sequence
                if current_sequence:
                    new_sequence += current_sequence
                    seq_slice.append(slice(ind - len(current_sequence), ind))
                    current_sequence = ""

        # reset the oxDNA coordinates if the sequence is not the same
        if len(new_sequence) != len(self._coords):
            self._coords = Coords(())

        # pick the directionality
        if directionality is None:
            directionality = self._sequence.directionality

        # assign the new sequence
        self._sequence = Sequence(
            new_sequence, directionality, callback=self._updated_sequence
        )

        self._seq_slice = seq_slice
        self._strand = strand_str

    ###
    ### METHODS
    ###

    def copy(self, callback=None) -> "Strand":
        """
        Create a copy of the current strand, including coordinates and attributes.

        Parameters
        ----------
        callback : callable, default None
            A callback function to be registered in the copied strand.

        Returns
        -------
        Strand
            A deep copy of the strand.
        """
        ### Initialize a new strand object
        new_strand = Strand.__new__(Strand)

        # attributes that are calculated fresh or immutable that
        # can be just reassigned
        attr_to_copy = {
            "_strand",  # strand properties
            "_seq_slice",
            "_start",  # 2D positions
            "_direction",
            "_positions",
            "_directions",
            "_seq_positions",
            "_prec_pos",
            "_next_pos",
            "_max_pos",
            "_min_pos",
        }
        for attr in attr_to_copy:
            setattr(new_strand, attr, getattr(self, attr))

        # sequence attributes
        new_seq = self._sequence.copy(callback=new_strand._updated_sequence)
        new_strand._sequence = new_seq

        new_strand._coords = self._coords.copy()
        new_strand._strands_block = StrandsBlock(new_strand)

        ### callbacks
        new_strand._callbacks = [callback] if callback else []

        ### pk_info
        new_strand.pk_info = copy.deepcopy(self.pk_info)

        return new_strand

    def draw(
        self, canvas: list = None, return_string: bool = True, plane: int = 0
    ) -> Optional[str]:
        """
        Draw the strand in 2D space using a canvas.
        The canvas is a list of strings, where each string is a line of the canvas.
        The strand is drawn using the characters and positions of the strand.

        Parameters
        ----------
        canvas : list, default None
            The canvas to draw the strand on. If None, a new canvas is created.
        return_string : bool, default True
            If True, return the drawn canvas as a string (for printing).
        plane : int, default 0
            The plane to draw the strand on (experimental). If 0, the strand
            is drawn in the default plane.

        Returns
        -------
        str or None
            The drawn canvas as a string if return_string is True, otherwise None.

        Raises
        -------
        MotifStructureError
            If the canvas is not large enough to draw the strand or a position is
            trying to overwrite a symbol in the canvas.
        """
        max_x = max([pos[0] for pos in self.positions], default=-1) + 1
        max_y = max([pos[1] for pos in self.positions], default=-1) + 1

        ### CHECK THE CANVAS OR CREATE IT ACCORDING TO THE MAXIMUM POSITIONS
        if not canvas:
            canvas = [" " * (max_x)] * max_y
        elif len(canvas) < max_y:
            raise MotifStructureError(
                f"Error while drawing the strands. The number "
                f"of lines in the canvas ({len(canvas)}) is "
                f"lower than the line required by the strand "
                f"({max_y})"
            )
        elif any([len(line) < max_x for line in canvas]):
            raise MotifStructureError(
                f"Error while drawing the strands. The number "
                f"of line characters in the canvas "
                f"({[len(line) for line in canvas]}) "
                f"is lower than characters required by the "
                f"strand ({max_x})"
            )

        ### ADD A SYMBOL AND CHECK THE CANVAS
        for pos, sym in zip(self.positions, self.strand):

            # experimental plane
            if plane > 0 and len(pos) < 3:
                continue

            canv_sym = canvas[pos[1]][pos[0]]

            # experimental plane
            if len(pos) > 2 and pos[2] != plane:
                continue

            # the strand is crossing a symbol in the canvas
            if canv_sym != " " and canv_sym not in "┼+" and sym not in "┼+":
                current_canvas = "\n".join(canvas)
                raise MotifStructureError(
                    f"Error while drawing the strands. "
                    f"'{sym}' at position {pos} is trying"
                    f" to overwrite the symbol '{canv_sym}'. "
                    f"Current drawing: \n{current_canvas}.\n"
                    f"Current strand: {self.strand}"
                )

            # the symbol replace the space in the canvas
            else:
                canvas[pos[1]] = (
                    canvas[pos[1]][: pos[0]] + sym + canvas[pos[1]][pos[0] + 1 :]
                )

        if return_string:
            return "\n".join(canvas)

    def extend(
        self,
        direction: Union[Position, Tuple[int, int]],
        until: Union[Position, Tuple[int, int]] = None,
        check: bool = True,
    ) -> None:
        """
        Extend the strand in the specified direction by a given length.

        Parameters
        ----------
        direction : Union[Position, Tuple[int, int]], default (1, 0)
            The direction to extend the strand in. It can be a Position object or
            a tuple of (x, y) coordinates. It must contain coordinates that are either
            1, -1 or 0, where 0 is the coordinate that does not change.
        until : Union[Position, Tuple[int, int]], default None
            The position until which to extend the strand. If None, the strand is
            extended until the origin position (0, 0).
        check : bool, default True
            If True, check the direction and until position for validity.
        """
        if check:
            direction = self._check_position(direction, direction=True)
        if until is None:
            until = Position.zero()
        elif check:
            until = self._check_position(until)

        dir_symb = ["─", "│", "•"]

        # don't extend terminal symbols
        if self._strand[0] in "35":
            start_mask = Position.zero()
        else:
            start_mask = Position(
                d1 if d1 == d2 else 0 for d1, d2 in zip(direction, -self._direction)
            )
        if self._strand[-1] in "35":
            end_mask = Position.zero()
        else:
            end_mask = Position(
                d1 if d1 == d2 else 0 for d1, d2 in zip(direction, self.end_direction)
            )

        # if 1 not in start_mask and 1 not in end_mask:
        #     return self

        # normalize the direction according to until
        start_until = Position(
            u if u is not None else s for u, s in zip(until, self._start)
        )
        end_until = Position(u if u is not None else e for u, e in zip(until, self.end))

        # calculate the start and end amounts to extend
        start_amount = (start_until - self._start) * start_mask
        end_amount = (end_until - self.end) * end_mask

        for i, d in enumerate(direction):
            if d == 0:
                continue
            if start_amount[i] > 0:
                self._start += start_amount * direction
                self._strand = dir_symb[i] * start_amount[i] + self._strand
            if end_amount[i] > 0:
                self._strand += dir_symb[i] * end_amount[i]

        # update the strand and the sequence
        self._update_sequence_insertion(self._strand)
        # update the positions
        self._reset_positions()
        # notify the upper class that the strand has changed
        self._trigger_callbacks()

        return self

    def flip(
        self,
        horizontally: bool = True,
        vertically: bool = False,
        flip_start: bool = True,
    ) -> None:
        """
        Flip the 2D structure of the strand. This affects the orientation on the
        canvas but does not revert the sequence.

        Parameters
        ----------
        horizontally : bool, default=True
            Whether to flip the strand horizontally (mirror along the y-axis).
        vertically : bool, default=False
            Whether to flip the strand vertically (mirror along the x-axis).
        flip_start : bool, default=True
            Whether to also flip the starting point of the strand.
        """

        if flip_start:
            positions = self.positions

            # CALCULATE NEW START POSITION
            new_x, new_y = self._start

            if horizontally:
                new_x = max(positions, key=lambda x: x[0])[0] - self._start[0]
            if vertically:
                new_y = max(positions, key=lambda x: x[1])[1] - self._start[1]

            self._start = self._start.replace(x=new_x, y=new_y)

        # CALCULATE NEW DIRECTION and update the strand sybols
        if horizontally:
            self._direction = self._direction.replace(
                x=-self._direction[0], y=self._direction[1]
            )
            self._strand = self._strand.translate(horiz_flip)

        if vertically:
            self._direction = self._direction.replace(
                x=self._direction[0], y=-self._direction[1]
            )
            self._strand = self._strand.translate(verti_flip)

        self._reset_positions()
        # trigger the callbacks only once
        self._trigger_callbacks()

    def get_base_at_pos(self, pos: Union[Position, Tuple[int, int]]) -> str:
        """
        Get the base at the specified 2D position in the strand.

        Parameters
        ----------
        pos : Position or tuple of int
            The position to get the base from.

        Returns
        -------
        str
            The base at the specified position.
        """
        pos = self._check_position(pos)
        return str(self._sequence[self.seq_positions.index(pos)])

    def get_char_at_pos(self, pos: Union[Position, Tuple[int, int]]) -> str:
        """
        Get the character at the specified 2D position in the strand.

        Parameters
        ----------
        pos : Position or tuple of int
            The position to get the character from.

        Returns
        -------
        str
            The character at the specified position.
        """
        pos = self._check_position(pos)
        return self.strand[self._positions.index(pos)]

    def get_position_map(self) -> Dict[Tuple[int, int], str]:
        """
        Get a dictionary of positions as keys and the corresponding
        characters as values.

        Returns
        -------
        Dict[Tuple[int, int], str]
            A dictionary mapping positions to characters.
        """
        return {pos: char for pos, char in zip(self.positions, self.strand)}

    def insert(self, idx: int, val: str) -> None:
        """
        Insert a character into the strand at the specified index.

        Parameters
        ----------
        idx : int
            The index at which to insert the character.
        val : str
            The character(s) to insert into the strand.
        """
        val = self._check_line(val)
        strand_line = list(self.strand)
        strand_line.insert(idx, val)
        # this trigger the callbacks
        self.strand = "".join(strand_line)

    def invert(self) -> "Strand":
        """
        Invert the start and end positions of the strand,
        effectively reversing the sequence directionalty
        and 3D coordinates orientation.

        Returns
        -------
        Strand
            The updated strand instance.
        """
        # don't trigger the callbacks, the strand is the same for the Motif
        seq_len = len(self.sequence)
        strand_len = len(self.strand)

        # update the start positions and direction
        end = self.end
        self._direction = -self.end_direction
        self._start = end

        # update the strand and sequence
        self._strand = self._strand[::-1]
        self._sequence._directionality = self._sequence._directionality[::-1]
        self._sequence._sequence = self._sequence._sequence[::-1]
        self._seq_slice = [
            slice(strand_len - sl.stop, strand_len - sl.start)
            for sl in self._seq_slice[::-1]
        ]

        # invert the positions if already built
        if self._positions:
            self._positions = self._positions[::-1]
            self._directions = tuple(-d for d in self._directions[::-1])
            self._seq_positions = self._seq_positions[::-1]
            self._prec_pos, self._next_pos = self._next_pos, self._prec_pos
            # max pos and min pos don't change

        # adjust the pseudoknots information
        if self.pk_info:
            new_ind_fwd = []
            for start, end in self.pk_info["ind_fwd"]:
                new_ind_fwd.append((seq_len - end, seq_len - start))
            self.pk_info["ind_fwd"] = new_ind_fwd

        # update the 3D coordinates
        self._coords.reverse_in_place()

        return self

    def join(
        self, other: "Strand", trigger_callbacks: bool = True
    ) -> Optional["Strand"]:
        """
        Join this strand with another strand, aligning their start and end
        positions and direction. This behaves like the @staticmethod join_strands,
        but it modifies the current strand instead of creating a new one.
        Check the join_strands method for more details.

        Parameters
        ----------
        other : Strand
            The other strand to join.
        trigger_callbacks : bool, default True
            If True, trigger the callbacks after joining the strands.

        Returns
        -------
        Strand or None
            The joined strand, or None if the strands cannot be joined.
        """
        if not isinstance(other, Strand):
            raise ValueError(
                f"The object to join is not a Strand object,"
                f" got {type(other)} instead."
            )
        # load the edges of the strands
        s1_start_prev = self.prec_pos
        s1_end_next = self.next_pos
        s2_start_prev = other.prec_pos
        s2_end_next = other.next_pos

        ### cases to invert the second strand
        s2_inverted = False

        # case 1: (1)-->,<--(2)
        if s1_end_next == other.end and s2_end_next == self.end:
            other.invert()
            s2_inverted = True

        # case 2: <--(1),-->(2)
        elif s1_start_prev == other.start and s2_start_prev == self.start:
            other.invert()
            s2_inverted = True

        # now the strands are in the right position to be joined
        ### case in which I can just add
        join_order = None

        # case 3: (1)-->,-->(2)
        if s1_end_next == other.start:
            join_order = 1

        # case 4: <--(1),<--(2)
        elif s1_start_prev == other.end:
            join_order = -1

        if join_order is not None:
            self._check_addition(other, copy=False)

            if join_order == 1:
                first, second = self, other
            else:
                first, second = other, self

            if not self._sequence:
                self._sequence._directionality = other._sequence.directionality

            self._start = first._start
            self._direction = first._direction

            self._coords = Coords.combine_coords(first, second)

            # this calls the reset of all the positions
            self._update_sequence_insertion(first.strand + second.strand)

            if self.sequence:
                self.strands_block.update(
                    self.strands_block, other.strands_block, avoid_strands={other}
                )

            if first._positions and second._positions:
                self._positions = first._positions + second._positions
                self._directions = first._directions + second._directions
                self._seq_positions = first._seq_positions + second._seq_positions
                self._prec_pos = first._prec_pos
                self._next_pos = second._next_pos
                self._max_pos = Position(
                    (
                        max(first._max_pos[0], second._max_pos[0]),
                        max(first._max_pos[1], second._max_pos[1]),
                    )
                )
                self._min_pos = Position(
                    (
                        min(first._min_pos[0], second._min_pos[0]),
                        min(first._min_pos[1], second._min_pos[1]),
                    )
                )

            # calculte the new pseudoknots information
            self.pk_info = first._combine_pk_info(second)

        # revert the strand to the original position
        if s2_inverted:
            other.invert()

        if join_order is None:
            return

        # trigger the callbacks if needed
        if trigger_callbacks:
            self._trigger_callbacks()
        return self

    def pop(self, idx: int) -> str:
        """
        Remove and return a character from the strand at the specified index.

        Parameters
        ----------
        idx : int
            The index of the character to remove.

        Returns
        -------
        str
            The removed character.
        """
        strand_line = list(self.strand)
        popped_val = strand_line.pop(idx)
        self.strand = "".join(strand_line)  # this trigger the callbacks
        return popped_val

    def reverse(self) -> "Strand":
        """
        Reverse the directionality of the strand sequence
        (5' to 3' becomes 3' to 5' and vice versa).

        Returns
        -------
        Strand
            The updated strand instance.
        """

        self._sequence.reverse(inplace=True)
        return self

    def save_3d_model(
        self,
        filename: str = "strand",
        return_text: bool = False,
        pdb: bool = False,
        **kwargs,
    ) -> Optional[Tuple[str, str]]:
        """
        Save the 3D structure of the strand in oxDNA format (.dat, .top)
        and optionally as PDB.

        Parameters
        ----------
        filename : str, default='strand'
            The base name for the output files.
        return_text : bool, default=False
            If True, returns the configuration and topology text instead
            of saving to file (used for real-time visualization).
        pdb : bool, default=False
            If True, exports the structure as a PDB file (requires
            `oxDNA_analysis_tools`).
        **kwargs : dict
            Additional keyword arguments passed to the PDB export.

        Returns
        -------
        tuple of str or None
            Returns (conf_text, top_text) if `return_text` is True, otherwise None.
        """
        ### check for the sequence and directionality
        if self.directionality == "53":
            seq = str(self.sequence)
            coords = self.coords
        else:
            seq = str(self.sequence[::-1])
            coords = self.coords.reverse()

        ### initialize the strand and nucleotide length
        seq_len = len(seq)
        n_strands = 1

        ### initialize the configuration and topology text
        conf_text = "t = 0\nb = 100 100 100\nE = 0 0 0\n"
        top_text = seq + " type=RNA circular=false \n"

        ### Buil the configuration text
        for pos, a1, a3 in coords:
            conf_text += (
                f"{pos[0]} {pos[1]} {pos[2]} "
                f"{a1[0]} {a1[1]} {a1[2]} "
                f"{a3[0]} {a3[1]} {a3[2]}\n"
            )

        ### ADD THE PROTEINS
        for protein in coords.proteins:
            for pos, a1, a3 in protein.coords:
                conf_text += (
                    f"{pos[0]} {pos[1]} {pos[2]} "
                    f"{a1[0]} {a1[1]} {a1[2]} "
                    f"{a3[0]} {a3[1]} {a3[2]}\n"
                )

            ### ADD THE PROTEIN TO THE TOPOLOGY
            top_text += protein.sequence + " type=peptide circular=false \n"
            seq_len += len(protein)
            n_strands += 1

        top_text = f"{seq_len} {n_strands} 5->3\n" + top_text
        if return_text:
            return conf_text, top_text

        ### write the oxDNA file
        filename = filename.split(".")[0]  # remove the extension from the filename
        conf_file = f"{filename}.dat"
        with open(conf_file, "w", encoding="utf-8") as f:
            f.write(conf_text)

        ### write the top file
        top_file = f"{filename}.top"
        with open(top_file, "w", encoding="utf-8") as f:
            f.write(top_text)

        ### write the pdb file
        if pdb:
            if not oat_installed:
                warnings.warn(
                    "oxDNA_analysis_tools is not installed. " "Skipping PDB export.",
                    UserWarning,
                )
                return

            # Read oxDNA configuration
            system, _ = strand_describe(top_file)
            ti, di = describe(top_file, conf_file)
            conf = get_confs(ti, di, 0, 1)[0]
            conf = inbox(conf, center=True)

            oxDNA_PDB(conf, system, filename, **kwargs)

    def shift(
        self,
        shift_position: Union["Position", "Direction", Tuple[int, int]],
        check: bool = True,
    ) -> "Strand":
        """
        Shift the strand position in the 2D layout.

        Parameters
        ----------
        shift_position : Position or Direction or tuple of int
            The shift to apply in x and y direction.
        check : bool, default True
            If True, check if the shift is valid (i.e. not negative coordinates).

        Returns
        -------
        Strand
            The updated strand instance.
        """
        shift = self._check_position(shift_position)
        self._start += shift

        # don't update the positions if the strand is not built
        if self._positions is None:
            self._trigger_callbacks()
            return self

        # apply the shift to all the positions
        positions = tuple(pos + shift for pos in self._positions)

        # check if the shift is valid
        if check and any(coor < 0 for pos in positions for coor in pos):
            raise ValueError(
                f"The strand reaches negative coordinates: "
                f"The positions would be: {positions}."
            )

        # directly update the positional attributes
        self._positions = positions
        self._seq_positions = tuple(pos + shift for pos in self._seq_positions)
        self._prec_pos += shift
        self._next_pos += shift
        self._max_pos += shift
        self._min_pos += shift

        self._trigger_callbacks()
        return self

    def transform(self, transformation_matrix: Any) -> None:
        """
        Apply a transformation matrix to the strand's 3D coordinates.
        Check the Coords class for the transformation matrix format.

        Parameters
        ----------
        transformation_matrix : array-like
            The transformation matrix to apply.
        """
        self.coords.transform(transformation_matrix)


class StrandsBlock(set):
    """
    A container class that holds a collection of Strand instances.

    This class inherits from Python's built-in list and is designed to manage
    a group of Strand objects. It provides methods for adding, removing,
    joining, and transforming strands in the block.
    All the strands in the block have locked 3D coordinates respect to
    each other. This means that when a strand in the block is transformed,
    all the other strands in the block are transformed as well with the same
    transformation matrix. This is foundamental for managing 3D structures
    of motifs.

    Parameters
    ----------
    *args : Strand
        Strand instances to be added to the block.
    **kwargs : Any
        Additional keyword arguments passed to the list constructor.

    Attributes
    ----------
    items : list
        The list of Strand instances stored in the block.
    """

    def __init__(self, *args: "Strand") -> None:
        """
        Initialize a new StrandsBlock with Strand instances.

        Parameters
        ----------
        *args : Strand
            Strand instances to be added to the block.

        Raises
        ------
        ValueError
            If any of the input arguments is not an instance of Strand.
        """
        for s in args:
            if not isinstance(s, Strand):
                raise ValueError("The input arguments must be Strand instance")
            else:
                self.add(s)

    def __add__(self, other: "StrandsBlock") -> "StrandsBlock":
        """
        Concatenate two StrandsBlock instances.

        Parameters
        ----------
        other : StrandsBlock
            The StrandsBlock to be added.

        Returns
        -------
        StrandsBlock
            A new StrandsBlock instance with the combined strands.
        """
        return StrandsBlock(*(array for array in super().__add__(other)))

    def add(self, item: "Strand") -> None:
        """
        Append a Strand instance to the StrandsBlock.

        Parameters
        ----------
        item : Strand
            The Strand instance to append.

        Raises
        ------
        ValueError
            If the item is not an instance of Strand.
        """
        if not isinstance(item, Strand):
            raise ValueError(f"The item must be a Strand instance, got {type(item)}")

        # append the strand to the block if it has coordinates
        if not item.coords.is_empty():
            super().add(item)
            item._strands_block = self

    def remove(self, item: "Strand") -> None:
        """
        Remove a Strand instance from the StrandsBlock.

        Parameters
        ----------
        item : Strand
            The Strand instance to remove.

        Raises
        ------
        ValueError
            If the item is not an instance of Strand.
        """
        if not isinstance(item, Strand):
            raise ValueError(f"The item must be a Strand instance, got {type(item)}")

        super().remove(item)
        item.strands_block = None

    def transform(self, T: Any) -> None:
        """
        Apply a transformation matrix to all the strands in the StrandsBlock.
        For the transformation matrix format, check the Coords class and its
        transform method.

        Parameters
        ----------
        T : np.ndarray
            The transformation matrix to apply to each strand.

        Notes
        -----
        This method assumes that the transformation matrix is compatible
        with the strand's coordinate system.
        """
        for strand in self:
            strand.transform(T)

    def update(
        self,
        *others: List["StrandsBlock"],
        avoid_strands: Optional[Set["Strand"]] = None,
    ) -> None:
        """
        Update the strands block with the strands in the input iterable.
        This method is similar to the add method, but it allows for
        adding multiple strands at once. It also allows for avoiding
        specific strands in the update process.

        Parameters
        ----------
        *other : StrandsBlock
            The StrandsBlock instances to update with.
        avoid_strands : set of Strand, default None
            A set of strands to avoid in the update process.

        Raises
        ValueError
            If any of the input arguments is not a StrandsBlock, or
            if avoid_strand is not a Strand.
        """
        if not all(isinstance(s, StrandsBlock) for s in others):
            raise ValueError("The input arguments must be StrandsBlock instances")

        # check if the avoid_strand is a Strand instance
        if avoid_strands is not None and any(
            not isinstance(s, Strand) for s in avoid_strands
        ):
            raise ValueError("The avoid_strands must be a set of Strand instances")

        if avoid_strands is None:
            avoid_strands = set()

        # add the strands to the block
        for sb in others:
            for s in sb:
                if s not in avoid_strands:
                    self.add(s)
