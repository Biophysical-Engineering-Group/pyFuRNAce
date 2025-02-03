import warnings
import copy
from .symbols import *
from .callback import Callback
from .position import Position, Direction
from .sequence import Sequence
from .coordinates_3d import Coords
try:
    from oxDNA_analysis_tools.UTILS.RyeReader import get_confs, describe, strand_describe, inbox
    from oxDNA_analysis_tools.oxDNA_PDB import oxDNA_PDB
    oat_installed = True
except ImportError:
    oat_installed = False

class Strand(Callback):

    ### 
    ### MAGIC METHODS
    ###

    def __init__(self, strand: str = '', directionality: str = '53', start: tuple = (0,0), direction: Direction=Direction.RIGHT, coords: Coords=None, strands_block = None, **kwargs):
        super().__init__(**kwargs) # register the callback
        ### strings properties
        self._sequence = Sequence(sequence='', directionality=directionality, callback=self._updated_sequence)
        strand = strand.upper()
        if not isinstance(strand, Sequence):
            self._check_line(strand)
        self._update_sequence_insertion(strand)

        ### check the 2D direction 
        self._check_position(direction, True)
        self._direction = direction.value if isinstance(direction, Direction) else tuple(direction)

        ### check the start position
        if start == "min":
            self._start = self.minimal_dimensions[1]
        else:
            self._check_position(start)
            self._start = tuple(start)
        ### map properties
        self._reset_maps()
        if coords is None:
            self._coords = Coords()
        else:
            self.coords = coords
        if strands_block is not None:
            self.strands_block = strands_block
        else:
            self.strands_block = StrandsBlock(self)

    def __repr__(self):
        return self.strand
    
    def __str__(self):
        return self.strand
        
    def __len__(self):
        # return len(self.strand)
        return len(self.strand)
    
    def __getitem__(self, idx):
        return self.strand[idx]

    def __setitem__(self, idx, val):
        val = val.upper()
        # convert to list to accept slice indexing and negative indexing
        new_strand_line = list(self.strand)
        new_strand_line[idx] = val
        new_strand_str = "".join(new_strand_line)
        self.strand = new_strand_str ### update the strand

    def __add__(self, other):
        other = self._check_addition(other)
        coords = Coords.combine_coords(self, self.coords, other, other.coords)
        new_strand = Strand(str(self)+ str(other), 
                            self.sequence.directionality, 
                            self.start, self.direction, 
                            callbacks=self._callbacks,
                            coords=coords)
        new_pk_info = self._combine_pk_info(other)
        if new_pk_info:
            new_strand.pk_info = new_pk_info
        return new_strand
    
    def __iadd__(self, other):
        other = self._check_addition(other)
        self.strand = self.strand + other.strand
        coords = Coords.combine_coords(self, self.coords, other, other.coords)
        self._coords = coords
        self._trigger_callbacks()
        new_pk_info = self._combine_pk_info(other)
        if new_pk_info:
            self.pk_info = new_pk_info
        return self

    def __radd__(self, other):
        if other == 0:
            return self
        elif isinstance(other, str):
            return Strand(other + str(self), directionality=self.directionality, start=self.start, direction=self.direction, callbacks=self._callbacks)
        elif isinstance(other, Sequence):
            return Strand(str(other) + str(self), directionality=self.directionality, start=self.start, direction=self.direction, callbacks=self._callbacks)
        else:
            return self.__add__(other)
        
    def __eq__(self, other):
        """ Check that the two strands have same map"""
        if isinstance(other, Strand):
            return self.map == other.map and self.directionality == other.directionality
        return False
    
    def __bool__(self):
        """ Return true if the strand has strand characters"""
        return bool(self.strand)
    
    def __contains__(self, other):
        """Check if a position or a substrand in included in the strand."""
        # check a position (a tuple or list containing two integers)
        if isinstance(other, (tuple, list)) and len(other) == 2 and isinstance(other[0], int) and isinstance(other[1], int):
            return other in self.map
        if isinstance(other, Sequence):
            return other in self.sequence
        if isinstance(other, str):
            return other in str(self.strand)
        # check a substrand
        if isinstance(other, Strand):
            for key, value in other.map.items():
                if key not in self.map or self.map[key] != value:
                    return False
            return True
        return False

    ### 
    ### PROPERTIES
    ###

    @property
    def strand(self):
        """ The strand build from the symbols and sequences"""
        return self._strand
    
    @strand.setter
    def strand(self, new_strand):
        new_strand = new_strand.upper()
        self._check_line(new_strand)
        self._update_sequence_insertion(new_strand)
        self._reset_maps()
        self._trigger_callbacks()

    @property
    def sequence(self):
        """ The complete sequence of the strand"""
        return self._sequence
    
    @sequence.setter
    def sequence(self, new_seq):
        """ Set the sequence of the strand. If the sequence is a string, it is split into the sequences of the strand.
        If the new sequence is longer than the current sequence, the sequence is added at the end of the strand."""
        #takes the sequence and splits it up into the sequences of specific length making up the sequence list
        new_seq = new_seq.upper()
        self._check_line(new_seq)
        self._updated_sequence(new_seq)

    @property
    def directionality(self):
        """ The directionality of the sequence"""
        return self.sequence.directionality
    
    @directionality.setter
    def directionality(self, new_directionality):
        """ Set the directionality of the sequence"""
        self.sequence.directionality = new_directionality
        self._trigger_callbacks()

    @property
    def sequence_list(self):
        """ The list of sequences of the strand"""
        seq_ind = 0
        seq_list = []
        for sl in self._seq_slice: # iterate over the slices of the sequence
            sequence_increment = sl.stop - sl.start
            seq_list.append(self.sequence[seq_ind: seq_ind+sequence_increment])
            seq_ind += sequence_increment
        return seq_list
        
    @property
    def symbols(self):
        """ The symbols of the strand"""
        return self.strand.translate(nucl_to_none)
    
    @property 
    def sequence_slice(self):
        """ The list of slices of the sequence in the strand"""
        return self._seq_slice

    @property
    def start(self):
        """ Start position of the strand in the 2d representation: a tuple of x,y coordinates. """
        return self._start

    @start.setter
    def start(self, new_start):
        """ Start position setter: check that the inpur is a tuple/list of x,y coordinates. """
        self._check_position(new_start)
        self._start = tuple(new_start)
        self._reset_maps()
        self._trigger_callbacks()

    @property
    def direction(self):
        """ The starting (x, y) 2D direction to build the strand in 2D. """
        return self._direction
    
    @direction.setter
    def direction(self, new_direction):
        """ The starting direction setter: a tuple of x,y coordinates in which one coordinate is 0 and the other is +- 1. """
        self._check_position(new_direction, True)
        self._direction = new_direction.value if isinstance(new_direction, Direction) else tuple(new_direction)
        self._reset_maps()
        self._trigger_callbacks()
    
    @property
    def map(self):
        """ A dictionary of x,y coordinates for each character of the strand (ordered according to the character order)"""
        if not self._map:
            self._calculate_map()
        return self._map

    @property
    def direction_map(self):
        """ The map of the strand exiting direction for each position"""
        if not self._direction_map:
            self._calculate_map()
        return self._direction_map
    
    @property
    def base_map(self):
        """ A list of x,y coordinates for each nucleotide of the strand (the map attribute filtered to include only nucleotide)"""
        if self._base_map:
            return self._base_map
        self._base_map = dict()
        for pos, sym in self.map.items():
            if sym in nucleotides:
                self._base_map[pos] = sym
        return self._base_map

    @property
    def end(self):
        """ The last position of the strand: it is the last element of the map"""
        if not self._map:
            self._calculate_map()
        return self._end # the last element is the end symbol
    
    @property
    def end_direction(self):
        """ The last direction of the strand"""
        return self.direction_map[self.end]
    
    @property
    def prev_pos(self):
        if not self._map:
            self._calculate_map()
        return self._prev_pos

    @property
    def next_pos(self):
        if not self._map:
            self._calculate_map()
        return self._next_pos
    
    @property
    def minimal_dimensions(self):
        """ Get the minimal size of the canvas and the suggested start (without taking into account structural errors)"""
        start = [0, 0]
        pos = [0,0]
        pos_max = [0, 0]
        direction = self.direction
        for sym in self:
            # check that the structure does not reach negative coordinates (x,y), in case add 1 to all positions
            for i in range(2):
                if pos[i] < 0:
                    for l in (start, pos, pos_max):
                        # add +1 to the element that is negative
                        l[i] += 1
                # check if we find a nex maximum element in x,y
                if pos[i] > pos_max[i]:
                    pos_max[i] = pos[i]
            ### UPDATE DIRECTIONS ###
            if sym in "╰╮\\":
                direction = direction[::-1]
            elif sym in "╭╯/":
                direction = [-x for x in direction[::-1]]
            # calculate new position
            pos[0] += direction[0]
            pos[1] += direction[1]
        return tuple(pos_max), tuple(start)
    
    @property
    def coords(self):
        """ The oxDNA coordinates of the strand"""
        if not self._coords.is_empty():
            return self._coords
        coords = Coords.compute_helix_from_nucl((0,0,0), # start position
                                            (1,0,0), # base vector
                                            (0,1,0), # normal vector
                                            length= len(self.sequence),
                                            directionality = self.directionality,
                                            double = False)
        self._coords = coords
        return coords
    
    @coords.setter
    def coords(self, new_coords):
        if not new_coords:
            self._coords = Coords()
        if new_coords.size > 0 and new_coords.shape and new_coords.shape[0] != len(self.sequence):
            raise MotifStructureError(f"The number of oxDNA coordinates ({new_coords.shape[0]}) is different from the number of nucleotides ({len(self.sequence)})")
        if not isinstance(new_coords, Coords):
            new_coords = Coords(new_coords)
        self._coords = new_coords

    @property
    def strands_block(self):
        """ The block of strands to which the strand belongs"""
        return self._strands_block
    
    @strands_block.setter
    def strands_block(self, new_block):
        if new_block is None:
            self._strands_block = StrandsBlock(self)
            return
        if not isinstance(new_block, StrandsBlock):
            raise ValueError(f"The strands block must be a StrandsBlock object. Got {type(new_block)} instead.")
        self._strands_block = new_block
        self._strands_block.append(self)

    ### 
    ### PROTECTED METHODS
    ###

    def _calculate_map(self):
        """ A function to calculate the map """
        pos = list(self.start)
        direction = self.direction
        self._prev_pos = (pos[0] - direction[0], pos[1] - direction[1])
        pos_dict = {}
        dir_dict = {}
        new_strand = []  # Use list for efficient string concatenation
        translated = self.strand.translate(symb_to_road)
        for sym in translated:
            # Convert the symbol to road symbols
            if sym == '/':
                if direction[0] == -1 or direction[1] == -1:
                    sym = '╭'
                elif direction[0] == 1 or direction[1] == 1:
                    sym = '╯'
            elif sym == '\\':
                if direction[0] == 1 or direction[1] == -1:
                    sym = '╮'
                elif direction[0] == -1 or direction[1] == 1:
                    sym = '╰'
            elif sym in ('↑', '↓'):
                if direction[1] == -1 and self.directionality == '53' or\
                    direction[1] == 1 and self.directionality == '35':
                    sym = '↑'
                else:
                    sym = '↓'
            
            # Append the new symbol to the strand list
            new_strand.append(sym)

            # Check for invalid positions and raise error early
            if pos[0] < 0 or pos[1] < 0:
                raise MotifStructureError(f"The strand reaches negative x,y coordinates: {pos}. The current map is: {pos_dict}. The last direction is: {direction}. Current strand: {''.join(new_strand)}")

            pos_tuple = tuple(pos)
            
            # Check symbols and directions
            if direction[0] and any((sym == "│", direction[0] == 1 and sym in "╰╭", direction[0] == -1 and sym in "╮╯")):
                raise MotifStructureError(f"The x direction of the strand (dx: {direction[0]}) is not compatible with the next symbol ({sym}). Current strand: {''.join(new_strand)}. Full strand: {self.strand}")
            elif direction[1] and any((sym == "─", direction[1] == 1 and sym in "╭╮", direction[1] == -1 and sym in "╰╯")):
                raise MotifStructureError(f"The y direction of the strand (dy: {direction[1]}) is not compatible with the next symbol ({sym}). Current strand: {''.join(new_strand)}. Full strand: {self.strand}")
            elif pos_tuple in pos_dict and sym not in "┼":
                warnings.warn(f"The strand is doing a crossing not allowed, '{sym}' is trying to overwrite the symbol '{pos_dict[pos_tuple]}'. The symbol will be overwritten with '┼'.", AmbiguosStructure, stacklevel=3)
                sym = '┼'
                new_strand[-1] = sym

            # Add symbol to the map
            pos_dict[pos_tuple] = sym

            # Update directions
            if sym in "╰╮\\":
                direction = (direction[1], direction[0])
            elif sym in "╭╯/":
                direction = (-direction[1], -direction[0])
            dir_dict[pos_tuple] = direction

            # Update positions
            pos[0] += direction[0]
            pos[1] += direction[1]

        # Update the strand symbols
        new_strand = ''.join(new_strand)
        if self._strand != new_strand:
            self._strand = new_strand
        self._map = pos_dict
        self._direction_map = dir_dict
        self._end = (pos[0] - direction[0], pos[1] - direction[1])
        self._next_pos = tuple(pos)

    @staticmethod
    def _check_position(input_pos_dir, direction=False):
        """ Check the input position or direction. If direction=True, also check the constrains for a direction"""
        if isinstance(input_pos_dir, (Direction, Position)):
            return True 
        if not (isinstance(input_pos_dir, (tuple, list)) and len(input_pos_dir) == 2 and isinstance(input_pos_dir[0], int) and isinstance(input_pos_dir[1], int)):
            raise ValueError(f'The 2D coordinates must be a tuple/list of (x,y) integer values. Got {input_pos_dir} instead.')
        if direction and len(set(input_pos_dir) & {-1, 0, 1}) != 2:
            raise ValueError(f'2D direction not allowed. The allowed values are: right (1, 0); bottom (0, 1); left (-1, 0); top (0, -1). Got {input_pos_dir} instead.')
        return True

    def _check_line(self, line, translator=symb_to_none):
        if isinstance(line, Sequence):
            return True
        if not isinstance(line, str):
            raise ValueError(f"The strand must be a string or a sequence object. Got {type(line)} instead.")
        for term_sym in ("3", "5"):
            if term_sym not in line:
                continue
            if line.count(term_sym) > 1:
                raise MotifStructureError(f"The strand can have only one start and one end symbol ('5' and '3'). Got {line.count(term_sym)} '{term_sym}' symbols.")
            if term_sym not in (line[0], line[-1]):
                raise MotifStructureError(f"The '5' and '3' symbols are terminal symbols. Got '{term_sym}' end at index '{line.index(term_sym)}' instead. Full strand: {line}.")
        allowed_symbols = road_symbols
        if translator == nucl_to_none:
            allowed_symbols = nucleotides
        if line.translate(translator):
            raise ValueError(f"The string '{line}' contains symbols not allowed in ROAD. The symbols allowed are: {allowed_symbols}.")
        #direction of the strand is saved in the first Seqeunce in the sequence list
        return True

    def _check_addition(self, other, copy=True):
        if isinstance(other, str):
            return Strand(other, self.sequence.directionality)
        elif isinstance(other, Sequence) and self.sequence and self.sequence.directionality != other.directionality:
            raise MotifStructureError(f'Cannot add a strand with a Sequence with different directionality. \n\tStrand: {self._strand}. \n\tStrand directionality: {self.directionality}\n\tSequence: {other}\n\tSequence directionality: {other.directionality}')
        elif isinstance(other, Sequence):
            return Strand(Sequence, self.directionality)
        elif not isinstance(other, Strand):
            raise ValueError(f'{other} is not a valid type for addition')
        elif self.sequence and other.sequence and self.directionality != other.directionality:
            raise MotifStructureError(f'Cannot add two strands with different directionality. \n\tFirst strand: {self._strand}. \n\tFirst strand directionality: {self.directionality}\n\tSecond strand {other._strand}\n\tSecond strand directionality: {other.directionality}')
        if copy:
            return other.copy()

    def _update_sequence_insertion(self, strand: str):
        """ Update the sequence insertion according to the symbols in the strand"""
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
                # if the current sequence is not empty and the current symbol is not a nucleotide, add it to the list of sequences
                if current_sequence:
                    new_sequence += current_sequence
                    seq_slice.append(slice(ind - len(current_sequence), ind))
                    current_sequence = "" # reset the current sequence
        if len(new_sequence) != len(self._sequence):
            self._coords = Coords(()) # reset the oxDNA coordinates
        self._sequence = Sequence(new_sequence, self._sequence.directionality, callback=self._updated_sequence)
        self._seq_slice = seq_slice
        self._strand = strand_str
        if not self._sequence:
            return
        # check the start and end symbols
        for term_sym in ("3", "5"):
            # the terminal symbol is not in the strand
            if term_sym not in self._strand:
                continue
            # identified a terminal symbol, check the directionality and revert it if necessary
            if term_sym == self._strand[0] and term_sym != self._sequence.directionality[0]:
                self._sequence.directionality = self._sequence.directionality[::-1]
                # warnings.warn(f"The strand starts with the '{term_sym}' symbol but the sequence directionality is '{self._sequence.directionality}'.", AmbiguosStructure, stacklevel=3)
            elif term_sym == self._strand[-1] and term_sym != self._sequence.directionality[-1]:
                self._sequence.directionality = self._sequence.directionality[::-1]
                # warnings.warn(f"The strand ends with the '{term_sym}' symbol but the sequence directionality is '{self._sequence.directionality}'.", AmbiguosStructure, stacklevel=3)

    def _updated_sequence(self, new_sequence, **kwargs):
        if isinstance(new_sequence, str):
            new_sequence = Sequence(new_sequence, self._sequence.directionality, callback=self._updated_sequence)
        if len(new_sequence) != len(self._sequence):
            self._coords = Coords(()) # reset the oxDNA coordinates
        self._sequence = new_sequence
        build_strand = self._strand.translate(nucl_to_none)
        seq_ind = 0
        for sl in self._seq_slice: # iterate over the slices of the sequence
            sequence_increment = sl.stop - sl.start # the length of the sequence slice
            # replace the sequence slice with the new sequence
            build_strand = build_strand[:sl.start] + str(new_sequence)[seq_ind: seq_ind+sequence_increment] + build_strand[sl.start:]
            seq_ind += sequence_increment
        # add the rest of the sequence
        build_strand += str(new_sequence)[seq_ind:]
        self._strand = build_strand
        self._reset_maps()
        self._trigger_callbacks(**kwargs)

    def _reset_maps(self):
        """ Reset the strand map, base map and direction map"""
        # Main map parameters
        self._map = None
        self._base_map = None
        self._direction_map = None
        # # Additional parameters calculated from the map
        # self._end = None
        # self._prev_pos = None
        # self._next_pos = None

    def _combine_pk_info(self, other):
        if not hasattr(self, 'pk_info') and not hasattr(other, 'pk_info'):
            return None
        new_pk_info = {"id": [], 'ind_fwd': [], 'E': [], 'dE': []}
        if hasattr(self, 'pk_info'):
            new_pk_info = copy.deepcopy(self.pk_info)
        if hasattr(other, 'pk_info'):
            offset = len(self.sequence)
            new_pk_info['id'] += other.pk_info['id']
            new_pk_info['ind_fwd'] += [(x[0] + offset, x[1] + offset) for x in other.pk_info['ind_fwd']]
            new_pk_info['E'] += other.pk_info['E']
            new_pk_info['dE'] += other.pk_info['dE']
        return new_pk_info

    ### 
    ### STATIC METHODS
    ###


    @staticmethod
    def join_strands(strand1, strand2):
        """ Join two consecutive strands together considering their end and start position.
        --------------------------------------------------------------------------------------
        strand2: Strand
            the strand2 strand to join
        """
        if not isinstance(strand2, Strand) and not isinstance(strand1, Strand):
            raise ValueError(f'The objects to join are not a Strand object, got {type(strand2)}, {type(strand2)} instead.')
        # there are four cases to check, the general approach is keeping the main start and end direction of the main strand
        s1_start_prev = strand1.prev_pos
        s1_end_next = strand1.next_pos
        # using the same approach for the strand2 strand to double check, otherwise we could have problems for joinins a strand with one symbol
        s2_start_prev = strand2.prev_pos
        s2_end_next = strand2.next_pos

        ### case in which I have to invert the second stand
        s2_inverted = False
        # case 1: (1)-->,<--(2)
        if s1_end_next == strand2.end and s2_end_next == strand1.end:
            strand2.invert()
            s2_inverted = True
        # case 2: <--(1),-->(2)
        elif s1_start_prev == strand2.start and s2_start_prev == strand1.start:
            strand2.invert()
            s2_inverted = True
        # now the strands are in the right position to be joined

        ### case in which I can just add
        # case 3: (1)-->,-->(2)
        if s1_end_next == strand2.start:
            strand1._check_addition(strand2, copy=False) # check right direction for addition
            # join the coordinates
            coords = Coords.combine_coords(strand1, strand1.coords, strand2, strand2.coords)
            # set the directionality of the strand with the sequence
            directionality = strand1.directionality if strand1.sequence else strand2.directionality
            # create a new joined strand
            joined_strand = Strand(strand1.strand + strand2.strand, directionality=directionality, start=strand1.start, direction=strand1.direction, coords=coords, callbacks = strand1.callbacks)
            
            # update the pseudoknots information
            new_pk_info = strand1._combine_pk_info(strand2)
        # case 4: <--(1),<--(2)
        elif s1_start_prev == strand2.end:
            strand1._check_addition(strand2, copy=False) # check right direction for addition
            # join the coordinates
            coords = Coords.combine_coords(strand2, strand2.coords, strand1, strand1.coords)
            # setting the strand set the callbacks
            directionality = strand1.directionality if strand1.sequence else strand2.directionality
            # create a new joined strand
            joined_strand = Strand(strand2.strand + strand1.strand, directionality=directionality, start=strand2.start, direction=strand2.direction, coords=coords, callbacks = strand1.callbacks)

            # update the pseudoknots information
            new_pk_info = strand2._combine_pk_info(strand1)
        else: # no joining possible
            if s2_inverted: # revert the strand to the original position
                strand2.invert()
            return None
        
        # update the pseudoknots information
        if new_pk_info:
            joined_strand.pk_info = new_pk_info

        if s2_inverted:
            strand2.invert()
        StrandsBlock.join_blocks(joined_strand.strands_block, strand1.strands_block, avoid_strand=strand1)
        if id(strand1.strands_block) != id(strand2.strands_block):
            StrandsBlock.join_blocks(joined_strand.strands_block, strand2.strands_block, avoid_strand=strand2)
        return joined_strand

    ### 
    ### METHODS
    ###

    def draw_strand(self, canvas: list = None, return_draw=True):
        """ Draw the 2D structure of the strand and draw it on a canvas (a list of strings). 
        If the canvas is not given, it is calculated considering the maximum x,y positions.
        --------------------------------------------------------------------------------------
        canvas: list
            The list of string on which write the strand 2D structure.
        --------------------------------------------------------------------------------------
        Return: list
            the modified canvas
        """
        map2d = self.map
        max_x = max([pos[0] for pos in map2d], default=-1) + 1
        max_y = max([pos[1] for pos in map2d], default=-1) + 1

        ### CHECK THE CANVAS OR CREATE IT ACCORDING TO THE MAXIMUM POSITIONS ###
        if not canvas:
            canvas = [" " * (max_x)] * max_y
        elif len(canvas) < max_y:
            raise MotifStructureError(f'Error while drawing the strands. The number of lines in the canvas ({len(canvas)}) is lower than the line required by the strand ({max_y})')
        elif any([len(line) < max_x for line in canvas]):
            raise MotifStructureError(f'Error while drawing the strands. The number of line characters in the canvas ({[len(line) for line in canvas]}) is lower than characters required by the strand ({max_x})')
        
        ### ADD A SYMBOL AND CHECK THE CANVAS ###
        for pos, sym in map2d.items():
            canv_sym = canvas[pos[1]][pos[0]]
            # the strand is crossing a symbol in the canvas
            if canv_sym != ' ' and canv_sym not in "┼+" and sym not in "┼+":
                current_canvas = '\n'.join(canvas)
                raise MotifStructureError(f"Error while drawing the strands. '{sym}' at position {pos} is trying to overwrite the symbol '{canv_sym}'. Current drawing: \n{current_canvas}.\nCurrent strand: {self.strand}")
            # the symbol replace the space in the canvas
            else:
                canvas[pos[1]] = canvas[pos[1]][:pos[0]] + sym + canvas[pos[1]][pos[0] + 1:] 
        
        if return_draw: 
            return '\n'.join(canvas)
    
    def flip(self, horizontally: bool = True, vertically: bool = False, flip_start: bool = True):
        """ Flip the structure 2D representation (doesn't reverse the sequence).
        --------------------------------------------------------------------------------------
        horizontally: bool, default True
            Flip the structure horizontally
        vertically: bool, default True
            Flip the structure vertically
        start: bool, default False
            Flip the start position too
        """
        if flip_start:
            map2d = self.map
            # calculate the new start position
            new_x, new_y = self._start
            if horizontally:
                new_x = max(map2d.keys(), key=lambda x: x[0])[0] - self._start[0]
            if vertically:
                new_y = max(map2d.keys(), key=lambda x: x[1])[1] - self._start[1]
            self._start = (new_x, new_y)
        if horizontally:
            self._direction = (-self._direction[0], self._direction[1])
            self._strand = self._strand.translate(horiz_flip) 
        if vertically:
            self._direction = (self._direction[0], -self._direction[1])
            self._strand = self._strand.translate(verti_flip)
        self._trigger_callbacks() # trigger the callbacks only once

    def invert(self):
        """ Invert the start/end sides to build the structure"""
        end = self.end
        self._direction = (- self.end_direction[0], - self.end_direction[1])
        self._start = end
        self._sequence._directionality = self._sequence._directionality[::-1]
        self._coords.reverse_in_place()
        self._update_sequence_insertion(self._strand[::-1])
        self._reset_maps()
        # don't trigger the callbacks, the strand is the same for the Motif perspective
        seq_len = len(self.sequence)

        # adjust the pseudoknots information
        if hasattr(self, 'pk_info'):
            new_ind_fwd = []
            for start, end in self.pk_info['ind_fwd']:
                new_ind_fwd.append((seq_len - end, seq_len - start))
            self.pk_info['ind_fwd'] = new_ind_fwd
                
        return self
    
    def reverse(self):
        """ Reverse the direction of the sequence of the strand"""
        self._sequence.reverse()
        return self
    
    def shift(self, x_y_tuple):
        """ Shift the strand in the 2D representation.
        --------------------------------------------------------------------------------------
        x_y_tuple: tuple
            the shift in the x and y direction
        """
        self._check_position(x_y_tuple)
        x, y = x_y_tuple.value if isinstance(x_y_tuple, Direction) else tuple(x_y_tuple)
        self._start = (self.start[0] + x, self.start[1] + y)
        self._reset_maps()
        self._trigger_callbacks()
        return self
    

    def insert(self, idx, val):
        """ Insert a sequence at index
        --------------------------------------------------------------------------------------
        idx: int
            the index at which adding the character
        val: str
            the characters to add
        """     
        val = val.upper()   
        self._check_line(val)
        strand_line = list(self.strand)
        strand_line.insert(idx, val)
        self.strand = "".join(strand_line) # this trigger the callbacks

    def pop(self, idx):
        """ Pop the element at index
        --------------------------------------------------------------------------------------
        idx: int
            the index at which adding the character
        val: str
            the characters to add
        """        
        strand_line = list(self.strand)
        popped_val = strand_line.pop(idx)
        self.strand = "".join(strand_line) # this trigger the callbacks
        return popped_val

    def join(self, other):
        """ Join two consecutive strands together considering their end and start position.
        --------------------------------------------------------------------------------------
        other: Strand
            the other strand to join
        """
        if not isinstance(other, Strand):
            raise ValueError(f'The object to join is not a Strand object, got {type(other)} instead.')
        # there are four cases to check, the general approach is keeping the main start and end direction of the main strand
        s1_start_prev = self.prev_pos
        s1_end_next = self.next_pos
        # using the same approach for the other strand to double check, otherwise we could have problems for joinins a strand with one symbol
        s2_start_prev = other.prev_pos
        s2_end_next = other.next_pos

        ### case in which I have to invert the second stand
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
        # case 3: (1)-->,-->(2)
        if s1_end_next == other.start:
            self._check_addition(other, copy=False) # check right direction for addition
            # join the coordinates
            coords = Coords.combine_coords(self, self.coords, other, other.coords)
            # setting the strand set the callbacks
            if not self._sequence: # if the sequence is empty, set the sequence of the other strand, so the direction is set correctly
                self._sequence = other.sequence
            self.strand = self.strand + other.strand
            self._coords = coords
            if self.strands_block != other.strands_block:
                StrandsBlock.join_blocks(self.strands_block, other.strands_block, avoid_strand=other)
            else:
                self.strands_block.remove(other)

            # calculte the new pseudoknots information
            new_pk_info = self._combine_pk_info(other)

        # case 4: <--(1),<--(2)
        elif s1_start_prev == other.end:
            self._check_addition(other, copy=False) # check right direction for addition
            # join the coordinates
            coords = Coords.combine_coords(other, other.coords, self, self.coords)
            # setting the strand set the callbacks
            self._start = other.start
            self._direction = other.direction
            if not self._sequence: # if the sequence is empty, set the sequence of the other strand, so the direction is set correctly
                self._sequence = other.sequence
            self.strand = other.strand + self.strand 
            self._coords = coords
            if self.strands_block != other.strands_block:
                StrandsBlock.join_blocks(self.strands_block, other.strands_block, avoid_strand=other)
            else:
                self.strands_block.remove(other)

            # calculte the new pseudoknots information
            new_pk_info = other._combine_pk_info(self)

        else: # no joining possible
            if s2_inverted: # revert the strand to the original position
                other.invert()
            return None

        # update the pseudoknots information
        if new_pk_info:
            self.pk_info = new_pk_info

        return self
    
    def go_through(self, pos):
        """Take a start or end position, and return the strand (in the direction that is walking 
        across the strand), and the next position in the strand."""
        if pos == self.start:
            return self.strand, (self.end[0] + self.end_direction[0], self.end[1] + self.end_direction[1])
        elif pos == self.end:
            return self.strand[::-1], (self.start[0] - self.direction[0], self.start[1] - self.direction[1])
        return None, None
    
    def isAtStartOrEnd(self, position):
        """
        Check if the given position is at the start or the end of the strand.

        Parameters:
        - position: A tuple representing the (x, y) coordinates of the position to check.

        Returns:
        - None if the position is neither at the start nor at the end.
        - 'start' if the position is at the start of the strand.
        - 'end' if the position is at the end of the strand.
        """
        self._check_position(position)
        if position == self.start:
            return 'start'
        elif position == self.end:
            return 'end'
        else:
            return None
        
    def transform(self, transformation_matrix):
        """ Apply a transformation matrix to the strand coordinates.
        """
        self.coords.transform(transformation_matrix)
        
    def save_3d_model(self, filename: str = 'strand', return_text: bool = False, pdb: bool = False, **kwargs):
        """ Save the oxDNA coordinates of the strand in a file.
        --------------------------------------------------------------------------------------
        filename: str
            The name of the file where to save the oxDNA coordinates
        """
        ### check for the sequence and directionality
        if self.directionality == '53':
            seq = str(self.sequence)
            coords = self.coords
        else:
            seq = str(self.sequence[::-1])
            coords = self.coords[::-1]

        ### initialize the strand and nucleotide length
        seq_len = len(seq)
        n_strands = 1

        ### initialize the configuration and topology text
        conf_text = f't = 0\nb = 100 100 100\nE = 0 0 0\n'
        top_text = seq + ' type=RNA circular=false \n'

        ### Buil the configuration text
        for pos, a1, a3 in coords:
            conf_text += f'{pos[0]} {pos[1]} {pos[2]} {a1[0]} {a1[1]} {a1[2]} {a3[0]} {a3[1]} {a3[2]}\n'
        ### ADD THE PROTEINS
        for protein in coords.proteins:
            for pos, a1, a3 in protein.coords:
                conf_text += f'{pos[0]} {pos[1]} {pos[2]} {a1[0]} {a1[1]} {a1[2]} {a3[0]} {a3[1]} {a3[2]}\n'
            ### ADD THE PROTEIN TO THE TOPOLOGY
            top_text += protein.sequence + ' type=peptide circular=false \n'
            seq_len += len(protein)
            n_strands += 1

        if return_text:
            return conf_text, top_text

        ### write the oxDNA file
        filename = filename.split('.')[0] # remove the extension from the filename
        conf_file = f'{filename}.dat'
        with open(conf_file, 'w') as f:
            f.write(conf_text)

        ### write the top file
        top_file = f'{filename}.top'
        with open(top_file, 'w') as f:
            f.write(f'{seq_len} {n_strands} 5->3\n' + top_text)

        ### write the pdb file
        if pdb:
            if not oat_installed:
                warnings.warn("oxDNA_analysis_tools is not installed. Skipping PDB export.", UserWarning)
                return
            # Read oxDNA configuration
            system, _ = strand_describe(top_file)
            ti, di = describe(top_file, conf_file)
            conf = get_confs(ti, di, 0, 1)[0]
            conf = inbox(conf, center=True)

            oxDNA_PDB(conf, system, filename, **kwargs)


    def copy(self, **kwargs):
        """ Create a fast copy of the strand"""
        new_stand = Strand(self.strand, self.directionality, start=self.start, direction=self.direction, coords=self.coords.copy(), **kwargs)
        for attr in self.__dict__.keys():
            if attr not in new_stand.__dict__.keys():
                setattr(new_stand, attr, copy.deepcopy(getattr(self, attr)))
            # don't copy the list because this copies the callbacks
            elif type(getattr(self, attr)) in (int, float, str, bool, tuple, dict):
                setattr(new_stand, attr, copy.copy(getattr(self, attr)))
        return new_stand


class StrandsBlock(list):

    def __init__(self, *args, **kwargs):
        for s in args:
            if not isinstance(s, Strand):
                raise ValueError("The input arguments must be Strand instance")   
            else:
                self.append(s)

    def __add__(self, other):
        return StrandsBlock(*(array for array in super().__add__(other)))

    def __setitem__(self, key, value):
        if not isinstance(value, Strand):
            raise ValueError("The item must be a Strand instance") 
        super().__setitem__(key, Coords(value))
        self[key].strands_block = self

    @staticmethod
    def join_blocks(s_block1, s_block2, avoid_strand=None):
        """ Join two strands blocks together"""
        if not isinstance(s_block1, StrandsBlock) or not isinstance(s_block2, StrandsBlock):
            raise ValueError("The input arguments must be StrandsBlock instances")
        if avoid_strand is not None and not isinstance(avoid_strand, Strand):
            raise ValueError("The avoid_strand must be a Strand instance")
        # avoid the strand to join
        s_block1 = [strand for strand in s_block1 if id(strand) != id(avoid_strand)]
        s_block2 = [strand for strand in s_block2 if id(strand) != id(avoid_strand)]
        return StrandsBlock(*s_block1, *s_block2)

    def append(self, item):
        if not isinstance(item, Strand):
            raise ValueError(f"The item must be a Strand instance, got {type(item)}")
        # append the strand to the block if it is not already in the block
        if id(item) not in [id(strand) for strand in self]:
            super().append(item)
            self[-1].strands_block = self

    def remove(self, item):
        if not isinstance(item, Strand):
            raise ValueError(f"The item must be a Strand instance, got {type(item)}")
        else:
            for i, strand in enumerate(self):
                if id(strand) == id(item):
                    del self[i]
                    strand.strands_block = None

    def transform(self, T):
        """Apply the transformation matrix T to all the strands in the strands_block"""
        for strand in self:
            strand.transform(T)