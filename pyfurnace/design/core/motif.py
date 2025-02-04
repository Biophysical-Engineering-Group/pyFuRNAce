from math import copysign
import warnings
import numpy as np
from inspect import signature
try:
    from oxDNA_analysis_tools.external_force_utils.force_reader import write_force_file
    from oxDNA_analysis_tools.external_force_utils.forces import mutual_trap
    from oxDNA_analysis_tools.UTILS.RyeReader import get_confs, describe, strand_describe, inbox
    from oxDNA_analysis_tools.oxDNA_PDB import oxDNA_PDB
    oat_installed = True
except ImportError:
    oat_installed = False
from .symbols import *
from .callback import Callback
from .sequence import Sequence
from .strand import Strand, StrandsBlock
from .basepair import BasePair
from .position import Position

class Motif(Callback):
    """
    Parent-class of all of the different motifs (eg. isdovetail, kissing loop,...) making up an RNA-Origami
    """

    def __init__(self, 
                 strands: list = None, 
                 *args, 
                 basepair: dict = None, 
                 structure: str = None, 
                 autopairing: bool = True, 
                 lock_coords: bool = True, 
                 join: bool = True, 
                 copy: bool = False, 
                 **kwargs):
        '''
        Attribute of elements of the class Motif.
        ---------------------------------------------------------------------
        basepair: dict
            dictionary with positions as key and symbols of base pair as value
        strands: list strands
            strands to add to the structure
        '''
        super().__init__(**kwargs)

        ### Initialize the attributes
        self._junctions = {}
        self._sequence = ''
        self._map = {}
        self._base_map = {}
        self._sequence_index_map = {}
        self._pair_map = {}
        self._structure = None
        self._strands = []
        self.lock_coords = lock_coords
        self._strands_block = StrandsBlock() # arrays of the 3d coordinates of the motif

        ### Add the strands to the motif
        if strands is None:
            strands = []
        elif isinstance(strands, Strand):
            strands = [strands]
        all_strands = list(strands) + list(args)
        if not all(isinstance(s, Strand) for s in all_strands):
            raise ValueError(f"All the elements in the strands input must be of type Strand. Got types {[type(s) for s in all_strands]}.")

        if copy:
            all_strands = self.copy_strands_preserve_blocks(all_strands, self)
        else:
            for s in all_strands:
                s.register_callback(self._updated_strands)

        # the consecutive strands are joined and copied into the motif _strand list
        if join:
            all_strands = self.join_strands(strands = all_strands)
        self._strands = all_strands # set the strands

        if self.lock_coords:
            # IMPORTANT: look at the _coords to avoid locking strands that don't have coordinates (therefore should not be locked)
            self._strands_block = StrandsBlock(*[s for s in self._strands if not s._coords.is_empty()])

        ### Add the basepair to the motif
        ### if the basepair dictionary is given, the autopairing is off and viceversa
        if basepair and structure:
            raise ValueError("Both basepair and structure cannot be given at the same time.")
        if structure:
            self.structure = structure
            autopairing = False
        else:
            if basepair:
                # convert the basepair dictionary to a dictionary of Position objects
                autopairing = False
            else:
                autopairing = autopairing
                basepair = BasePair()
            # initialize basepair property
            self.basepair = BasePair(basepair, callback = self._updated_basepair)
        # set autopairing
        self.autopairing = autopairing


    ### 
    ### CLASS METHODS (construct an object without using __init__())
    ###

    @classmethod
    def from_list(cls, motif_list, **kwargs):
        """ Obtain a Motif object from a list of strings with a motif text.
        Each Strand is read starting from the 5."""
        ### COLLECT THE START POSITION ALONG THE LIST OF THE MOTIF ###
        strand_list = [] # list of strands
        mapped_pos = set() # dictionary with the position of the strands
        # Get the maximum x and y positions
        max_x = max([len(line) for line in motif_list], default=1) 
        max_y = len(motif_list) - 1
        ### Uniform the line length
        for index, line in enumerate(motif_list):
            diff = max_x - (len(line) - 1)
            if diff:
                motif_list[index] = line + ' ' * diff
        for y, line in enumerate(motif_list):
            for x, char in enumerate(line):
                if char in "^v><": # Special direction characters
                    if char == '^':
                        direction = (0, -1)
                    elif char == 'v':
                        direction = (0, 1)
                    elif char == '>':
                        direction = (1, 0)
                    elif char == '<':
                        direction = (-1, 0)
                    start_pos = (x + direction[0], y + direction[1])
                    ### TRACE THE STRANDS ### 
                    # trace the strand to obtain the character of strand
                    strand_chars = Motif.trace(motif_list, start_pos, direction, limits = (max_x, max_y))
                    # create the Strand object
                    strand = Strand(strand_chars, directionality='53', start=start_pos, direction=direction)
                    # add the Strand to the list
                    strand_list.append(strand)
                    mapped_pos.update(set(strand.map.keys()))
                elif char == "5" and (x, y) not in mapped_pos: # got a start_position, now find the direction of the strand
                    right_sym = None
                    bot_sym = None
                    left_sym = None
                    top_sym = None
                    if x + 1 <= max_x:
                        right_sym = motif_list[y][x + 1].capitalize()
                    if y + 1 <= max_y:
                        bot_sym = motif_list[y + 1][x].capitalize()
                    if x > 0:
                        left_sym = motif_list[y][x - 1].capitalize()
                    if y > 0:
                        top_sym = motif_list[y - 1][x].capitalize()
                    direction = None
                    if right_sym in road_symbols and motif_list[y][x+1] not in ' 3':
                        direction = (1, 0)
                    elif bot_sym in road_symbols and bot_sym not in ' 3':
                        direction = (0, 1)
                    elif left_sym in road_symbols and left_sym not in ' 3':
                        direction = (-1, 0)
                    elif top_sym in road_symbols and top_sym not in ' 3':
                        direction = (0, -1)
                    if direction is None:
                        raise MotifStructureError(f"Strand at position {x, y} has no direction.")
                    start_pos = (x + direction[0], y + direction[1])

                    ### TRACE THE STRANDS ### 
                    # trace the strand to obtain the character of strand
                    strand_chars = Motif.trace(motif_list, start_pos, direction, limits = (max_x, max_y))
                    # create the Strand object
                    strand = Strand(strand_chars, directionality='53', start=start_pos, direction=direction)
                    # add the Strand to the list
                    strand_list.append(strand)
                    mapped_pos.update(set(strand.map.keys()))

        ### TRACE THE BASEPAIRS ###
        basepair = BasePair()
        for y, line in enumerate(motif_list):
            for x, sym in enumerate(line):
                pos1 = None
                # if position corresponds to a basepair symbol, add it to the dictionary
                if sym in "┊!*:":
                    pos1 = (x,y-1) 
                    pos2 = (x,y+1)
                elif sym in "=":
                    pos1 = (x-1,y)
                    pos2 = (x+1,y)
                if pos1 is not None:
                    basepair[Position(pos1)] = Position(pos2)
        return cls(strands=strand_list, basepair=basepair, **kwargs)

    @classmethod
    def from_text(cls, motif_text, **kwargs):
        """ Obtain a Motif object from a text with the motif sketcehd"""
        motif_list = [line for line in motif_text.split('\n')]
        return cls.from_list(motif_list, **kwargs)
    
    @classmethod
    def from_file(cls, file_path, **kwargs):
        """ Build a motif object from a textfile with the motif sketcehd"""
        with open(file_path, 'r') as f:
            motif_text = f.read()
        return cls.from_text(motif_text, **kwargs)
    
    @classmethod
    def concat(cls, motifs, *args, axis: int = 1, extend = False, copy=True, align=True, position_based=False, align_junctions: list=None, unlock_strands=False, **kwargs):
        """Concatenate the motifs (horizotally or vertically) by stacking them.
        extend: bool, default False
            If True, the junctions are extended to account for the shift.
        align: bool, default True
            If True, the motifs are aligned before concatenation.
        align_junctions: list of tuple, default None
            List of tuples with the index of junctions to align for each couple of motif. The junctions are given as tuples with the index of the junction in the first motif and the index of the junction in the second motif.
        position_based: bool, default False
            The motifs are shifted based on the junctions position if set to False. If set to True, the motifs are shifted based on the minimum and maximum positions of the motifs.
            If no junctions are detected, the motifs are shifted based on the positions.
        """
        if isinstance(motifs, Motif):
            motifs = [motifs]
        # put the motifs in a list and remove the empty motifs
        motifs = [m for m in motifs if m] + [m for m in args if m]
        if copy:
            motifs = [m.copy() for m in motifs]
        if align:
            aligned = Motif.align(motifs, axis=axis, align_junctions=align_junctions)
        else :
            aligned = motifs
        x_pos = axis
        y_pos = int(not axis)
        for ind, m1 in enumerate(aligned[:-1]):
            m2 = aligned[ind+1]
            if extend:
                if ind == 0:
                    m1.extend_junctions(skip_axis=y_pos)
                m2.extend_junctions(skip_axis=y_pos)
            if position_based: # if position based, calculate the shift based on the positions
                max_pos_m1 = m1.max_pos[y_pos]
                min_pos_m2 = m2.min_pos[y_pos]
            else:
                m1_junct = m1.junctions[(x_pos,y_pos)]
                m2_junct = m2.junctions[(-x_pos,-y_pos)]
                if m1_junct and m2_junct: # if there are junctions, calculate the shift based on the junctions
                    ind1, ind2 = 0, 0
                    if align_junctions and align_junctions[ind]:
                        ind1, ind2 = align_junctions[ind]
                    max_pos_m1 = m1_junct[ind1][y_pos]
                    min_pos_m2 = m2_junct[ind2][y_pos]
                else: # just concatenate the motifs
                    max_pos_m1 = m1.max_pos[y_pos]
                    min_pos_m2 = 0
                    if m2.min_pos[y_pos] > max_pos_m1: 
                        # if the second motif is already shifted, don't shift
                        max_pos_m1 = -1

            shift_pos = max_pos_m1 - min_pos_m2 + 1
            m2.shift((shift_pos * axis, shift_pos * int(not axis))) # move motif to the right of the first motif
        basepair = BasePair()
        autopairing = all([m.autopairing for m in aligned])
        # if all the motifs have autopairing on, the new motif has autopairing on
        if not autopairing:
            for m in aligned:
                    basepair.update(m.basepair)
        # if unlock strands is True, unlock all the strands
        if unlock_strands:
            for m in aligned:
                for s in m:
                    s.strands_block = None
        new_motif = cls(**kwargs)
        new_motif.autopairing = autopairing
        new_motif._basepair = basepair
        new_motif.replace_all_strands([strand for m in aligned for strand in m], copy=False, join=True)
        return new_motif

    ### 
    ### STATIC METHODS
    ###

    @staticmethod
    def trace(motif_list, pos, direction, limits):
        """ Trace the strands in a motif_list """
        # pos is a tuple (x,y) of the character we are considering
        # so it is a tuple (char_ind, line_ind)
        x = pos[0]
        y = pos[1]
        # we reach the edges of the list
        if pos[0] > limits[0] or pos[0] == -1 or pos[1] > limits[1] or pos[1] == -1:
            return ""
        # we have reached a new character
        sym = motif_list[y][x]
        
        # you reach a 5' or 3' end
        if sym == '3':
            return sym
        elif sym in bp_symbols:
            raise MotifStructureError("The strand leads to a base pairing")
        elif sym in "╰╮\\": # reverse the tuple
            direction = direction[::-1]
        elif sym in "╭╯/": # reverse the tuple and change sign
            direction = (- direction[1], - direction[0])
        elif sym == ' ' or sym not in road_symbols:
            return ""
        # if sym in 'AUCGNKSW─-|│+┼': keep the same direction
        pos = (pos[0] + direction[0], pos[1] + direction[1])
        return sym + Motif.trace(motif_list, pos, direction, limits)

    @staticmethod
    def join_strands(strands):
        """ Try to join consecutive strands. Return the list of joined strands. """
        joined_strands = set()
        strands = [s for s in strands if s] # remove the empty strands
        # ordering strand is:
        # - the strands with 5' end
        # - strand with the lowest y start position
        # - strand with the lowest x start position
        strands = sorted(strands, key=lambda s: (-int('5' in s.strand), *s.start[::-1]))
    
        ind1 = 0
        for ind1, s1 in enumerate(strands[:-1]):
            if ind1 in joined_strands: # if the strand is empty, or is already joined skip it
                continue 
            # join the strand with the consecutive strands until no more strands are joined

            for ind2, s2 in enumerate(strands):
                # skip the second strand if: it's empty, it's the same strand, the strands are already joined
                if ind1 == ind2 or ind2 in joined_strands: 
                    continue 

                joined = Strand.join_strands(s1, s2) # try to join the strands
                if joined is not None: # if the strands are joined
                    joined_strands.add(ind2) # add the index of the second strand to the joined strands
                    strands[ind1] = joined # replace the first strand with the joined strand
                    s1 = joined # update the first strand

        return [s for i, s in enumerate(strands) if i not in joined_strands] # return the strands that are not joined

    @staticmethod
    def align(motifs, *args, axis: int = 1, extend: bool = False, align_junctions:list=None, align_to:str="first", **kwargs):
        """Align the motif horizontally (axis = 1) or vertically (axis = 0) by shifting the motifs.
        If extend is True, the junctions are extended to account for the shift.
        --------------------------------------------------------------------------------------
            motifs: list of Motif
                list of motifs to align
            axis: int, default 1
                The axis along which the motifs are aligned. 1 for horizontal alignment, 0 for vertical alignment.
            extend: bool, default False
                If True, the junctions are extended to account for the shift.
            align_junctions: list of tuple, default None
                List of tuples with the index of junctions to align for each couple of motif. The junctions are given as tuples with the index of the junction in the first motif and the index of the junction in the second motif.
            align_to: str, default "first"
                How to aling the motif, the options are "first" to align to the first two junctions, "last" to align to the last two junctions, "center" to align to the center of the motifs.
                If no junctions are detected, the motifs are aligned to "center"
                """
        ### get the alignment direction
        if axis in (1, 0):
            direction = (axis, int(not axis))
        else:
            raise ValueError(f'{axis} is not a valid value for the axis parameter. The axis parameter must be 0 or 1')
        # the direction to shift is the opposite of the direction of the axis
        shift_direction = (int(not axis), axis)

        ### check alignment type
        if align_to not in ("first", "last", "center"):
            raise ValueError(f'{align_to} is not a valid value for the align_to parameter. The align_to parameter must be "first", "last" or "center"')

        ### check if the motifs are of the right type
        if issubclass(type(motifs), Motif):
            motifs = [motifs]
        motifs = motifs + list(args)

        ### aling all the motifs
        ind1 = 0
        ind2 = 1
        n_motifs = len(motifs)
        connected_motifs = [True] * n_motifs # start with the assumption that all the motifs are connected
        while ind2 < n_motifs: # stop when the second motif is the last one
            m1 = motifs[ind1]
            m2 = motifs[ind2]
            # skip the empty motifs
            if not m2:
                ind2 += 1
                continue
            if not m1:
                ind1 += 1
                continue
            if align_junctions and align_junctions[ind1]:
                # align the junction with the given index
                junct_ind1, junct_ind2 = align_junctions[ind1]
            else:
                # aling the first junction of the first motif with the first junction of the second motif
                junct_ind1, junct_ind2 = 0, 0 # by default align the first junction
                if align_to == "last":
                    junct_ind1, junct_ind2 = -1, -1
            # m1._check_addition(m2, direction)
            junctions1 = m1.junctions[direction]
            junctions2 = m2.junctions[(-axis, -int(not axis))]
            ### calculate the shift of the two motifs and shift the motif with the lower junction
            if align_to == "center":
                # if align_to is center align to the center of the motifs
                shift = int((m1.max_pos[axis] + m1.min_pos[axis]) / 2  - (m2.max_pos[axis] + m2.min_pos[axis]) / 2)
            elif junctions1 and junctions2 and junctions1[junct_ind1] and junctions2[junct_ind2]:
                # if junctions are detected, calculate the shift based on the junctions
                shift = junctions1[junct_ind1][axis] - junctions2[junct_ind2][axis]
            else: # nothing to align, the motifs are not connected
                connected_motifs[ind1] = False # mark the first motif as not connected
                shift = 0

            if not shift:
                pass
            elif shift > 0:
                # shift only the newly appended motif
                m2.shift((shift * shift_direction[0], shift * shift_direction[1]), extend=extend)
            else:
                # shift the left motif and all the previous one if they are connected
                for j in range(ind1, -1, -1):
                    if not connected_motifs[j]:
                        break
                    if motifs[j]:
                        motifs[j].shift((-shift * shift_direction[0], -shift * shift_direction[1]), extend=extend)
            # update the indices
            ind1 += 1
            ind2 += 1
        return motifs
    
    @staticmethod
    def get_sequential_shift(motifs, *args, axis: int = 1, position_based=True, **kwargs): 
        """Calculate the shift needed to align the motifs sequentially, based on the junctions or the positions of the motifs.
        Example:
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
        motifs += list(args)
        x_pos = axis
        y_pos = int(not axis)
        shifts = [0]
        current_shift = 0
        for ind, m1 in enumerate(motifs[:-1]):
            m2 = motifs[ind+1]
            ### calculate the shift based on the junctions or the positions of the motifs
            m1_junct = m1.junctions
            m2_junct = m2.junctions
            if position_based: # calculate the shift based on the positions 
                max_pos_m1 = m1.max_pos[y_pos]
                min_pos_m2 = m2.min_pos[y_pos]
            else:
                # if the junctions are detected, calculate the shift based on the junctions
                if m1_junct[(x_pos,y_pos)] and m2_junct[(-x_pos,-y_pos)]:
                    max_pos_m1 = max((pos[y_pos] for pos in m1.junctions[(x_pos,y_pos)]))
                    min_pos_m2 = min((pos[y_pos] for pos in m2.junctions[(-x_pos,-y_pos)]))
                else: # just concatenate the motifs
                    max_pos_m1 = m1.max_pos[y_pos]
                    min_pos_m2 = 0

            current_shift += max_pos_m1 - min_pos_m2 + 1
            shifts.append(current_shift)
        return shifts

    ### 
    ### MAGIC METHODS
    ###
    
    def __getitem__(self, idx):
        """ Get the strand at index """
        return self._strands[idx]

    def __setitem__(self, idx, strand):
        """ Set the strand at index """
        if not isinstance(strand, Strand):
            raise ValueError(f"{strand} is not a Strand object.")
        strand.register_callback(self._updated_strands)
        self._strands[idx] = strand
        self._updated_strands()
        self._strands = self.join_strands(self._strands)
        if self.lock_coords:
            self._strands_block = StrandsBlock(*[s for s in self._strands if not s._coords.is_empty()])

    def __len__(self):
        """ Get the number of strands"""
        return len(self._strands)
    
    def __repr__(self):
        """ Return the list of strands"""
        return str(list(self))
    
    def __str__(self):
        if not self:
            return ''
        # produce the canvas
        canvas_repr = [' ' * self.num_char] * self.num_lines 
        # draw each strand to the canvas
        for strand in self:
            strand.draw_strand(canvas_repr, return_draw=False)
        # add base pairing at position (if the base pairing position is free)
        for pos1, pos2 in self.basepair.items():
            if pos1[1] - pos2[1] == 0 and abs(pos1[0] - pos2[0]) == 2: # select the horizontal basepair symbol
                pos = ((pos1[0] + pos2[0])//2, pos1[1])
                sym = '='
            elif  pos1[0] - pos2[0] == 0 and abs(pos1[1] - pos2[1]) == 2: # select the vertical basepair symbol
                pos = (pos1[0], (pos1[1] + pos2[1])//2)
                sym = '┊'
            else:
                continue
            if pos[1] < len(canvas_repr) and pos[0] < len(canvas_repr[0]) and canvas_repr[pos[1]][pos[0]] == ' ': # if the position is free, add the basepair symbol
                canvas_repr[pos[1]] = canvas_repr[pos[1]][:pos[0]] + sym + canvas_repr[pos[1]][pos[0]+1:] 
            # else: # send a warning if the position is already occupied
            #     warnings.warn(f'Hidden basepair at position {pos}. Trying to pair {pos1} with {pos2}.', stacklevel=2)
        return '\n'.join(canvas_repr)

    def __add__(self, other):
        '''
        Addition of two elements of class Motif, were a new element of class Motif is created. 
        Elements are added by stacking horizontally.
        --------------------------------------------------------------------------------------
        other_motif: Motif
            element which is added (if it is of other type error occurs)
        --------------------------------------------------------------------------------------
        return: Element of type Motif, in which the the two elements have been added.
        '''
        self._check_addition(other)
        if not self and not other:
            return Motif()
        elif not self:
            return other.copy()
        elif not other:
            return self.copy()
        # create a copy of the motif to add. All operations will be performed on the copies
        self_copy = self.copy() 
        other_copy = other.copy()
        self.align([self_copy, other_copy])  # align the motifs horizontally
        max_right_x = max((pos[0] for pos in self_copy.junctions[(1,0)]))
        min_left_x = min((pos[0] for pos in other_copy.junctions[(-1,0)]))
        x_shift = max_right_x - min_left_x + 1
        other_copy.shift((x_shift, 0)) # move motif to the right of the first motif
        # if one of the two motifs doesn't use autopairing, set it false and update the basepair dictionary
        if self_copy.autopairing and other_copy.autopairing:
            new_basepair = BasePair()
        else:
            new_basepair = self_copy.basepair
            new_basepair.update(other_copy.basepair)
        # return the new motif collecting the strands
        return Motif(strands = list(self_copy)+list(other_copy), basepair=new_basepair) 
    
    def __iadd__(self, other):
        self._check_addition(other)
        if not other:
            return self
        # copy the other motif to prevent modifying it 
        other_copy = other.copy()
        if self:
            self.align([self, other_copy])  # align the motifs horizontally
            max_right_x = max((pos[0] for pos in self.junctions[(1,0)]))
            min_left_x = min((pos[0] for pos in other_copy.junctions[(-1,0)]))
            other_copy.shift((max_right_x - min_left_x + 1, 0)) # move motif to the right of the first motif
        # if other has autopairing off, set it off
        if not other_copy.autopairing:
            self.autopairing = False
            self._basepair.update(other_copy.basepair)
        # add the other strands to self
        for strand in other_copy:
            self._strands.append(strand)
        # join the strands
        self._strands = self.join_strands(self._strands)
        if self.lock_coords:
            self._strands_block = StrandsBlock(*[s for s in self._strands if not s._coords.is_empty()])
        for s in self:
            s.register_callback(self._updated_strands)
        self._updated_strands()
        return self
     
    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)
        
    def __bool__(self):
        """ Return False there are no strands or all the strands are False"""
        if not self._strands:
            return False
        for s in self._strands:
            if s:
                return True
        return False
    
    def __eq__(self, other):
        """ Check if the motifs are equal """
        if not isinstance(other, Motif):
            return False
        if len(self) != len(other):
            return False
        for s1, s2 in zip(self, other):
            if s1 != s2:
                return False
        if self.basepair != other.basepair:
            return False
        return True
    
    # def copy(self, **kwargs):
    #     """ Return a copy of the motif """
    #     copied = super().copy(**kwargs)
    #     if copied.lock_coords:
    #         copied._strands_block = StrandsBlock(*[s for s in self._strands if not s._coords.is_empty()])
    #     return copied
        
    ### 
    ### PROPERTIES
    ###

    @property
    def strands(self):
        """ Return the list of strands """
        return self._strands

    @property
    def num_lines(self):
        """
        The number of lines in the motif structure.
        """
        if not self._map:
            self.map
        return self._max_pos[1] + 1
        
    @property
    def num_char(self):
        """
        The number of characters in the structure lines
        """
        if not self._map:
            self.map # calculate the map
        return self._max_pos[0] + 1
    
    @property
    def max_pos(self):
        """
        The maximum position of the motif
        """
        if not self._map:
            self.map # calculate the map
        return self._max_pos
    
    @property
    def min_pos(self):
        """The minimum position occupied of the motif"""
        if not self._map:
            self.map # calculate the map
        return self._min_pos
    
    @property
    def map(self):
        """ 
        A dictionary with the position occupied as key and the strand index as value.
        """
        ### the map is already calculated and the strands aren't changed, return the map
        if self._map:
            return self._map

        ### update the map
        max_pos = []
        min_pos = []
        map_dict = {}
        for ind, strand in enumerate(self):
            for pos in strand.map:
                map_dict[pos] = ind
                ### initialize the maximum and minimum position 
                if not max_pos:
                    max_pos = list(pos)
                    min_pos = list(pos)
                ### update the maximum and minimum position
                if pos[0] > max_pos[0]:
                    max_pos[0] = pos[0]
                elif pos[0] < min_pos[0]:
                    min_pos[0] = pos[0]
                if pos[1] > max_pos[1]:
                    max_pos[1] = pos[1]
                elif pos[1] < min_pos[1]:
                    min_pos[1] = pos[1]
        self._map = map_dict
        if not max_pos:
            self._max_pos = (-1,-1)
        else:
            self._max_pos = tuple(max_pos)
        if not min_pos:
            self._min_pos = (0,0)
        else:
            self._min_pos = tuple(min_pos)
        return map_dict
    
    @property
    def base_map(self):
        """
        A dictionary with the position occupied as key and the strand index as value.
        """
        if self._base_map:
            return self._base_map
        
        base_map = {}
        for i, s in enumerate(self._strands):
            base_map.update({pos: i for pos in s.base_map})
        self._base_map = base_map
        return base_map
    
    @property
    def sequence_index_map(self):
        """
        A dictionary with the position occupied as key and the index of the sequence as value.
        """
        if self._sequence_index_map:
            return self._sequence_index_map
        
        pos_to_index = {} # convert the positions to indexes
        strand_offset = 0

        for s in self:
            # check the direction of the strand
            search_dir = 1 if s.directionality == '53' else -1

            # create the map from the positions to the indexes
            for ind, pos in enumerate(list(s.base_map)[::search_dir]):
                pos_to_index[pos] = strand_offset + ind

            # store the break points of the strands
            strand_offset += len(s.base_map)
        
        self._sequence_index_map = pos_to_index
        return pos_to_index
    
    @property
    def pair_map(self):
        """ Tha map of the paired indexes (alternative to the dot bracket notation) """
        if self._pair_map:
            return self._pair_map
        
        base_pair_dict = self.basepair
        ind_pair_map = {}
        for s in self:
            search_dir = 1 if s.directionality == '53' else -1
            for pos in list(s.base_map)[::search_dir]:
                paired_pos = base_pair_dict.get(pos, default=None)
                ind1 = self.sequence_index_map[pos]
                if paired_pos is not None: # ignore unpaired positions
                    ind2 = self.sequence_index_map[paired_pos]
                    ind_pair_map[ind1] = ind2
                else: # unpaired position
                    ind_pair_map[ind1] = None

        self._pair_map = ind_pair_map
        return ind_pair_map
        

    @property
    def junctions(self):
        '''
        A dictionary with the four junctions directions as key (right: (1,0), bottom: (0,1), left: (-1,0), top: (0,1)) and an ordered list of
        junction positon as value
        '''
        ### if the junctions are already calculated and the strands aren't changed, return them
        if self._junctions:
            return self._junctions
        
        ### update the junctions
        junctions_dict = {(1,0): [], (0,1): [], (-1,0): [], (0,-1): []} # initialize the dictionary
        #                  right,    bottom,     left,      top

        ### collect the junctions
        for s in self:
            if not s: # skip the empty strands
                continue
            # skip if the the start is '5' or '3' and the strand is not just the symbol
            if s.map[s.start] not in '35' or len(s.map) == 1: 
                junction_dir = (-s.direction[0], -s.direction[1]) # for the direction, invert the signes to have the direction of the junction
                junctions_dict[junction_dir].append(s.start) # add the position to the junction
            # skip if the the end is '5' or '3' and the strand is not just the symbol
            if s.map[s.end] not in '35' or len(s.map) == 1: 
                junctions_dict[s.end_direction].append(s.end)

        ### order the junctions
        for key, val in junctions_dict.items():
            # order the bottom/top junctions
            if key in ((0,1), (0,-1)):
                val.sort(key=lambda a : a[0]) # order the junctions from left to right
            # order the left/right junctions
            elif key in ((-1,0), (1,0)):
                val.sort(key=lambda a : a[1]) # order the junctions from top to bottom
        self._junctions = junctions_dict
        return junctions_dict

    @property
    def basepair(self):
        """
        A dictionary with positions as keys and basepairs symbols as value.
        """
        if self.autopairing and (not self._map or not self._sequence or not self._basepair):
            self._calculate_basepair() # calculate the basepair dictionary
            self._structure = None  # reset the dot bracket
        return self._basepair 
    
    @basepair.setter
    def basepair(self, basepair_dict):
        """
        Set the basepair dictionary. If the basepair dictionary is set, autopairing is turned off.
        """
        if not isinstance(basepair_dict, (dict, BasePair)):
            raise ValueError(f"{basepair_dict} must be a dictionary or a BasePair object. Got {type(basepair_dict)} instead.")
        try:
            basepair_dict = {Position(k): Position(v) for k, v in basepair_dict.items()}
        except Exception as e:
            raise ValueError(f"Error converting the basepair dictionary to a dictionary of Position objects: {e}")
        self.autopairing = False
        self._basepair = BasePair(basepair_dict, callback = self._updated_basepair)
        self._trigger_callbacks()
            

    @property
    def structure(self):
        """ Return the dot bracket representation of the motif """
        ### if the dot bracket is already calculated and the strands aren't changed, return it
        if self._structure:
            return self._structure
        
        # initialize variables to build the pair map

        ### create the pair_map from the positions to the index

        # calculate the break points of the strands
        break_points = []
        last_ind = 0
        for i, s in enumerate(self):
            if not s.sequence: # skip the strands without sequence
                continue
            last_ind += len(s.sequence)
            break_points.append(last_ind) # store the break points of the strands

        ### CREATE THE DOT BRACKET NOTATION ###
        dotbracket = pair_map_to_dot_bracket(self.pair_map, last_ind) # get the dot bracket notation

        for i, bp in enumerate(break_points[:-1]):
            dotbracket = dotbracket[:bp + i] + '&' + dotbracket[bp + i:]
        self._structure = dotbracket

        return self._structure
    
    @structure.setter
    def structure(self, structure: str):
        """ Set the dot bracket notation and the basepair dictionary"""
        if not isinstance(structure, str):
            raise ValueError(f"{structure} must be a string. Got {type(structure)} instead.")
        if len(structure) != len(self.sequence):
            raise ValueError(f"The length of the dot bracket must be equal to the length of the sequence. Got {len(structure)} for {len(self.sequence)}")
        basepair_dict = BasePair()
        pair_map = dot_bracket_to_pair_map(structure.replace('&', '')) # get the dictionary of connected indexes without  the separation symbol
        base_positions = list(self.base_map.keys()) # get the list of base positions in simple order
        for index, pos in enumerate(base_positions): # iterate over the base positions
            paired = pair_map[index]
            if paired is not None: # in the index of the base position is paired to something
                basepair_dict[Position(pos)] = Position(base_positions[paired])
        self.basepair = basepair_dict # set the basepair dictionary
        self._structure = structure
    
    @property
    def sequence(self):
        """ Return the sequence of the motif """
        ### if the sequence is already calculated and the strands aren't changed, return it
        if self._sequence:
            self._sequence 
        
        ### update the sequence
        tot_seq = ''
        for s in self:
            if not s.sequence: # skip the strands without sequence
                continue
            if s.directionality == '35':
                tot_seq += str(s.sequence)[::-1]
            else:
                tot_seq += str(s.sequence)
            tot_seq += '&'
        self._sequence = tot_seq.strip('&') # remove separator at the end
        return self._sequence
    
    @sequence.setter
    def sequence(self, seq_list: str = None):
        """ Set the sequence of each strand the motif """
        if not isinstance(seq_list, (str, Sequence)) and not isinstance(seq_list, (tuple, list)) and all(isinstance(s, (str, Sequence)) for s in seq_list):
            raise ValueError(f"{seq_list} must be a string, a Sequence object or a list of strings or Sequence objects. Got {type(seq_list)} instead.")
        if isinstance(seq_list, (str, Sequence)):
            seq_list = seq_list.split('&')
        if len(seq_list) != len(self):
            raise ValueError(f"The number of sequences must be equal to the number of strands. Got {len(seq_list)} sequences for {len(self)} strands.")
        seq_list = seq_list[:] # copy the list
        for s in self:
            s.sequence = seq_list.pop(0)
        self._sequence = None

    @property
    def lock_coords(self):
        """ Return the boolean indicating if the coordinates are locked """
        return self._lock_coords
    
    @lock_coords.setter
    def lock_coords(self, lock_coords):
        """ Set wether the coordinates of the motif must be locked """
        self._lock_coords = bool(lock_coords)
    
    ### 
    ### METHODS
    ###
    
    def flip(self, horizontally: bool = True, vertically: bool = True, strand_index: list = None):
        """ Flip the strands . Flip them horizontally (horizontally=True), vertically (vertically=True)
        or diagonally (horizontally= True, vertically= Ture)"""
        # save the initial index of character and lines, which change every time you change a strand
        idx_char = self.num_char - 1
        idx_lines = self.num_lines - 1
        # create new basepair dictionary in case autopairing is off
        new_basepair_dict = BasePair()
        for pos1, pos2 in self.basepair.items():
            if horizontally:
                pos1 = (idx_char - pos1[0], pos1[1]) # flip the horizontal position
                pos2 = (idx_char - pos2[0], pos2[1]) # flip the horizontal position
            if vertically:
                pos1 = (pos1[0], idx_lines - pos1[1]) # flip the vertical position
                pos2 = (pos2[0], idx_lines - pos2[1]) # flip the vertical position
            new_basepair_dict[Position(pos1)] = Position(pos2) # save the new basepair
        self._basepair = new_basepair_dict # save the new basepair
        for ind, strand in enumerate(self):
            if strand_index and ind not in strand_index: # if we want to flip only speficif strands, don't flip the unspecified strands
                continue
            # flip the start position of the strands: the new start is the border - the previous position
            new_start = list(strand.start)
            if horizontally:
                new_start[0] = idx_char - strand.start[0]
            if vertically:
                new_start[1] = idx_lines - strand.start[1]
            strand.start = new_start
            strand.flip(horizontally, vertically, flip_start=False) # flip the strand symbol
        return self
    
    def append(self, strand, join=True, copy=False):
        """ Add a strand to the motif. 
        --------------------------------------------------------------------------------------
        strand: Strand
            strand to add to the motif
        """
        if not isinstance(strand, Strand):
            raise ValueError(f"{strand} is not a Strand object.")
        if copy:
            strand = strand.copy()
        strand.register_callback(self._updated_strands)
        self._strands.append(strand)
        self._updated_strands()
        if join:
            self._strands = self.join_strands(self._strands)
        if self.lock_coords:
            self._strands_block = StrandsBlock(*[s for s in self._strands if not s._coords.is_empty()])
        return self
    
    def pop(self, index=-1):
        """ Remove the strand at index from the motif. 
        --------------------------------------------------------------------------------------
        index: int, default -1
            The index of the strand to remove from the motif.
        """
        strand = self._strands.pop(index)
        # REMOVE THE PAIRS IN WHICH THE STRAND IS INVOLVED
        new_basepair = BasePair()
        for k, v in self._basepair.items():
            if k in strand.map or v in strand.map:
                continue
            new_basepair[Position(k)] = Position(v)
        self._basepair = new_basepair
        strand._clear_callbacks()
        self._updated_strands()
        self._updated_basepair()
        if self.lock_coords:
            self._strands_block = StrandsBlock(*[s for s in self._strands if not s._coords.is_empty()])
        return strand
    
    def insert(self, index, strand, join=True, copy=False):
        """ Insert a strand at index in the motif. 
        --------------------------------------------------------------------------------------
        index: int
            The index where to insert the strand.
        strand: Strand
            The strand to insert in the motif.
        """
        if not isinstance(strand, Strand):
            raise ValueError(f"{strand} is not a Strand object.")
        if copy:
            strand = strand.copy()
        strand.register_callback(self._updated_strands)
        self._strands.insert(index, strand)
        self._updated_strands()
        if join:
            self._strands = self.join_strands(self._strands)
        if self.lock_coords:
            self._strands_block = StrandsBlock(*[s for s in self._strands if not s._coords.is_empty()])
        return self
    
    def replace_all_strands(self, strands, copy=True, join=False):
        """ Replace all the strand with a new list of strands
        --------------------------------------------------------------------------------------
        strands: list of Strand
            The new list of strands to insert in the motif.
        join: bool, default False
            If True, try to join the strands that are consecutive.
        """
        # copy the strands if needed
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
        # register the callback for each strand
        self._strands = new_strands
        if self.lock_coords:
            self._strands_block = StrandsBlock(*[s for s in self._strands if not s._coords.is_empty()])
        self._updated_strands()
        return self

    def rotate(self, times=1):
        """ Rotate the structure of 90 degree clockwise.
        --------------------------------------------------------------------------------------
        times: int, defualt 1
            The number of times to rotate the structure of 90 degrees.
        """
        for _ in range(times%4): 
            num_lines = self.num_lines
            for s in self:
                sign = +1
                if s._direction[1]: # when rotating from vertical to horizontal: change sign
                    sign = -1
                s._start = (num_lines - 1 - s.start[1], s.start[0])
                s._direction = (sign * s._direction[1], sign * s._direction[0]) 
                s.strand = s.strand.translate(rotate_90)   
            if not self.autopairing and self._basepair:
                new_bp = BasePair()
                for k, v in self._basepair.items():
                    new_k = (num_lines - 1 - k[1], k[0])
                    new_v = (num_lines - 1 - v[1], v[0])
                    new_bp[Position(new_k)] = Position(new_v)
                self._basepair = new_bp
        return self

    def shift(self, direction, extend=False):
        """ Shift the motif of (x,y) position. 
        --------------------------------------------------------------------------------------
        direction: tuple of int
            The number of positions to shift the motif in the x and y direction.
        """
        Strand._check_position(input_pos_dir = direction)
        try:
            direction = Position(direction)
        except Exception as e:
            raise ValueError(f"Error converting the direction to a Position object: {e}")
        min_pos = self.min_pos
        # check if the shift will bring the strands to negative positions
        if min_pos[0] + direction[0] < 0 or min_pos[1] + direction[1] < 0:
            raise ValueError(f'The motif cannot be shifed. The strands cannot be drawn at negative positons. Attempt to draw the motif at position {min_pos[0] + direction[0], min_pos[1] + direction[1]}')
        # check if there is a vertical shift
        junc_vert = []
        junc_hor = []
        if direction[1] > 0 and extend:
            junc_vert = self.junctions[(0, - int(copysign(1, direction[1])))] # take the junctions in the opposite direction of the vertical shift
        if direction[0] > 0 and extend:
            junc_hor = self.junctions[(- int(copysign(1, direction[1])), 0)] # take the junctions in the opposite direction of the horizontal shift
                
        for s in self:
            extend_strand = []
            # check if the strand has a junction in the opposite direction of the shift
            for pos in junc_vert:
                if pos in s.map: # the strand has a junction in the opposite direction of the vertical shift
                    extend_strand.append(Strand('│'*abs(direction[1]), start=(pos[0] + direction[0], pos[1]), direction=(0, int(copysign(1, direction[1])))))
            for pos in junc_hor:
                if pos in s.map: # the strand has a junction in the opposite direction of the horizontal shift
                    extend_strand.append(Strand('─'*abs(direction[0]), start=(pos[0], pos[1] + direction[1]), direction=(int(copysign(1, direction[1])), 0)))
            # for each strand shift the starting position
            # s.start = (s.start[0] + direction[0], s.start[1] + direction[1])
            s.shift(direction)
            #loop through the strands forming the extensions
            for s2 in extend_strand:
                s.join(s2) # join the strand to the vertical extension
            
        # shift every basepair too
        if not self.autopairing:
            self._basepair = self.basepair.shift(direction)
        return self
        
    def extend_junctions(self, skip_axis: int = None, until: tuple = (None, None)):
        """ Extend the junctions in order to reach the border of the motif. """
        juncts = self.junctions
        skip_direction = None
        if skip_axis in (0, 1):
            skip_direction = (skip_axis, int(not skip_axis))
        elif skip_axis:
            raise ValueError(f'{skip_axis} is not a valid value for the axis parameter. You can skip the axis 0 or 1')
        try:
            for axis, sym in ((1, '─'), (0, '│')): # consider the two axis
                naxis = int(not axis)
                pos_direction = (axis, naxis) 
                neg_direction = (- axis, - naxis)
                if pos_direction == skip_direction: 
                    continue # if the direction is the one to skip, skip it
                for pos in juncts[neg_direction]: # consider all the junction positions in the negative direction
                    main_dir = pos_direction
                    strand_at_pos = self[self.map[pos]] # get the strand at the position
                    if pos[naxis] != 0: # if the position is not at the minimum border
                        # add a strand of the symbol, starting from the position and going to the border, in the opposite direction
                        extend_strand = Strand(sym * pos[naxis], start=(pos[0] * naxis, pos[1] * axis), direction=pos_direction)
                        if not (set(extend_strand.map) & set(self.map)): # check that the extension doesn't overlap with other strands
                            strand_at_pos.join(extend_strand)
                for pos in juncts[pos_direction]:
                    main_dir = neg_direction
                    strand_at_pos = self[self.map[pos]] # get the strand at the position
                    max_pos = self.max_pos[naxis]
                    if until[naxis]: # if the until parameter is set to a specific position use it
                        max_pos = until[naxis] # set the maximum position to the until parameter
                    if pos[naxis] < max_pos: # if the position is not at the maximum border but can be extended to
                        # add a strand of the symbol, starting from the position and going to the border, in the opposite direction
                        extend_strand = Strand(sym * (max_pos - pos[naxis]), start=(pos[0] + 1 * axis, pos[1] + 1 * naxis), direction=pos_direction)
                        if not (set(extend_strand.map) & set(self.map)): # check that the extension doesn't overlap with other strands
                            strand_at_pos.join(extend_strand)
        except MotifStructureError as e:
            print(f"The junctions cannot be extended in the direction {main_dir}. The structure is not valid. Error:")
            print("\t", e)
        return self
    
    def strip(self, skip_axis=None):
        """ Remove the empty lines/columns in the motif structure. """
        min_pos = self.min_pos
        shift = (-min_pos[0], -min_pos[1])
        if skip_axis == 0:
            shift = (-min_pos[0], 0)
        elif skip_axis == 1:
            shift = (0, -min_pos[1])
        self.shift(shift)
        return self
    
    def sort(self, key=None, reverse=False):
        """ Sort the strands in the motif. """
        if not key:
            # sort the strand according to the lowest start position
            key = lambda s: s.start[::-1]
        self._strands.sort(key=key, reverse=reverse)
        return self
    
    def folding_barriers(self, kl_delay: int = 150, structure: str = None):
        ### kl_delay is the number of nts of delay before KLs snap closed in the Topology Checker section
	    ###  A 'realistic' setting might be around 350 ~ 1second.  A more conservative setting of 150 is set as the default.
        if structure is None:
            structure = self.structure
        ss_len = len(structure)
        s_map = dot_bracket_to_pair_map(structure)
        barriers = [''] * len(structure)

        # don't count the nucleotides in the 3' terminal position
        terminal_3_bonus = 16

        current_barr_count = 0
        for ind in range(ss_len):
            # set the default barrier   
            barriers[ind] = 'N'

            if ind > ss_len - terminal_3_bonus:
                if s_map[ind]:
                    barriers[s_map[ind]] = 'N'
                continue

            # get the bracket at the position
            bra = structure[ind]

            ### unpaired nts are not barriers and reset the current barrier count
            if bra == '.':
                barriers[ind] = 'N'
                current_barr_count = 0

            ### opening pairs may be barriers or not, indicate them with 'w', reset the barrier count
            elif bra == '(':
                barriers[ind] = 'w'
                current_barr_count = 0

            ### closing pairs may be barriers or not, depending on the topology count
            elif bra == ')':
                # if the opening pair is not blocked, we close and reset the current barrier count
                if barriers[s_map[ind]] == 'w':
                    barriers[ind] = 'N'
                    barriers[s_map[ind]] = 'N'
                    current_barr_count = 0

                # if the closing pair is already blocked, then we mark the barrier reached and start counting
                elif barriers[s_map[ind]] == 'S':
                    # if the current barrier count is greater than 5, wa mark it with 'S'
                    if current_barr_count > 5:
                        barriers[ind] = 'S'
                        barriers[s_map[ind]] = 'W'
                    # if the current barrier count is less than 5, we mark it with 'W'
                    else:
                        barriers[ind] = 'W'
                        barriers[s_map[ind]] = 'W'
                    current_barr_count += 1

            # if the index is bigger than the kl_delay, we check if the internal KLs are closed
            if ind > kl_delay:
                
                # we check if the KLs are closed, if yes, we mark them with 'N'
                close_sym = structure[ind - kl_delay]
                if close_sym not in '(.)' and close_sym in db_pairs.values():
                    barriers[ind - kl_delay] == 'N'
                    barriers[s_map[ind - kl_delay]] == 'N'

                    # for each nt in the delay, we check if there are any barriers, if yes, we mark them with 'S'
                    for k in range(s_map[ind - kl_delay], ind-kl_delay):
                        if barriers[k] == 'w':
                            barriers[k] = 'S'
                            
        penalty = {'S': 2, 'W': 1, 'w': 1, 'N': 0}
        repl = [penalty[i] for i in barriers]
        score = sum(repl)
        # print("FIX THIS FUNCTIONS!!!", end='\r')
        return ''.join(barriers), score

    
    def save_3d_model(self, filename: str = 'motif', 
                      config: bool = True,
                      topology: bool = True,
                      forces: bool = False, 
                      pk_forces: bool = False, 
                      return_text: bool = False, 
                      pdb=False, **kwargs):
        """ Save the motif in a oxDNA file format. """
        def get_kwargs_names(func):
            sig = signature(func)
            # Extract parameters that have default values (kwargs)
            kwargs = [param.name for param in sig.parameters.values() if param.default != param.empty]
            return kwargs

        filename = filename.split('.')[0] # remove the extension from the filename
        strands = [s for s in self if s.sequence]
        n_nucleotides = sum([len(s.sequence) for s in strands])
        n_strands = len(strands)
        # create the conformation and topology text
        conf_text = 't = 0\nb = 1000 1000 1000\nE = 0 0 0\n'
        topology_text = ''

        ### ADD THE STRANDS TO THE CONFORMATION AND TOPOLOGY TEXTS ###
        for s in strands:
            # check for the sequence direction
            if s.directionality == '53':
                seq = str(s.sequence)
                coord_array = s.coords.array
            else:
                seq = str(s.sequence[::-1])
                coord_array = s.coords[::-1]
            # add the coordinates to the conformations text
            if config:
                for pos, a1, a3 in coord_array:
                    conf_text += f'{pos[0]} {pos[1]} {pos[2]} {a1[0]} {a1[1]} {a1[2]} {a3[0]} {a3[1]} {a3[2]}\n'

            # add the sequence to the topology text
            if topology:
                topology_text += seq + ' type=RNA circular=false \n'

            # add the proteins to the conformation text
            try:
                for protein in s.coords.proteins:
                    if config:
                        for pos, a1, a3 in protein.coords:
                            conf_text += f'{pos[0]} {pos[1]} {pos[2]} {a1[0]} {a1[1]} {a1[2]} {a3[0]} {a3[1]} {a3[2]}\n'
                    # add the proteins to the topology text
                    if topology:
                        topology_text += f'{protein.sequence} type=peptide circular=false \n'
                    n_nucleotides += len(protein)
                    n_strands += 1
            except Exception as e:
                print('Problem with proteins', s, e)

        topology_text = f'{n_nucleotides} {n_strands} 5->3\n' + topology_text

        if return_text:
            return conf_text, topology_text

        # save the files
        conf_file = f'{filename}.dat'
        if config:
            with open(conf_file, 'w') as f:
                f.write(conf_text)
        top_file = f'{filename}.top'
        if topology:
            with open(top_file, 'w') as f:
                f.write(topology_text)

        ### save the external forces
        if not forces and not pk_forces:
            pass
        elif not oat_installed:
            warnings.warn("oxDNA_analysis_tools is not installed. Skipping force writing.", UserWarning)
        else:
            trap_kw_names = get_kwargs_names(mutual_trap)
            trap_kwargs = {k: v for k, v in kwargs.items() if k in trap_kw_names}
            pair_map = dot_bracket_to_pair_map(self.structure.replace('&', ''))
            trap_kwargs.setdefault('stiff', 0.09)
            trap_kwargs.setdefault('r0', 1.2)
            trap_kwargs.setdefault('PBC', True)
            trap_kwargs.setdefault('rate', 0)
            trap_kwargs.setdefault('stiff_rate', 0)
            force_list = []
            pk_force_list = []
            ss = self.structure.replace('&', '')
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
                if pk_forces and ss[n1] not in '.()':
                    pk_force_list.append(trap1)
                    pk_force_list.append(trap2)
            if forces:
                write_force_file(force_list, f'{filename}_forces.txt')
            if pk_forces:
                write_force_file(pk_force_list, f'{filename}_pk_forces.txt')

        if pdb:
            if not oat_installed:
                warnings.warn("oxDNA_analysis_tools is not installed. Skipping PDB export.", UserWarning)
            else:
                # Read oxDNA configuration
                system, _ = strand_describe(top_file)
                ti, di = describe(top_file, conf_file)
                conf = get_confs(ti, di, 0, 1)[0]
                conf = inbox(conf, center=True)
                # remove the proteins from the configuration if no pdb files are provided
                if not kwargs.get('protein_pdb_files'):
                    strand_offset = 0
                    to_pop = []
                    conf_to_keep = []
                    for i, strand in enumerate(system.strands):
                        strand_end = strand_offset + strand.get_length()
                        if strand.type == 'peptide':
                            to_pop.append(i)
                        else:
                            conf_to_keep.append((strand_offset, strand_end))
                        strand_offset = strand_end
                    for i in to_pop[::-1]:
                        system.strands.pop(i)

                oxdna_pdb_kw_names = get_kwargs_names(oxDNA_PDB)
                oxdna_pdb_kwargs = {k: v for k, v in kwargs.items() if k in oxdna_pdb_kw_names}
                oxdna_pdb_kwargs.setdefault('uniform_residue_names', True)
                oxDNA_PDB(conf, system, filename, **oxdna_pdb_kwargs)


    def copy(self, **kwargs):
        """ Return a copy of the motif. """
        return Motif(Motif.copy_strands_preserve_blocks(self._strands), basepair=self._basepair, autopairing=self.autopairing, copy=False, join=False, lock_coords=self.lock_coords, **kwargs)
                    

    ### 
    ### PROTECTED METHODS
    ###

    def _updated_strands(self, **kwargs):
        """ Update the motif when the strands are changed. 
        This is a callback function that is called when the strands are changed. 
        The callback has to contains *args and **kwargs to be compatible with the Callback class."""
        self._sequence = ''
        self._map = {}
        self._base_map = {}
        self._sequence_index_map = {}
        self._pair_map = {}
        self._max_pos = None
        self._min_pos = None
        self._junctions = {}
        if self.autopairing:
            self._basepair = BasePair()
            self._structure = None
        self._trigger_callbacks(**kwargs)

    def _updated_basepair(self, **kwargs):
        """ Update the basepair dictionary when the sequence is changed. 
        The callback has to contains *args and **kwargs to be compatible with the Callback class."""
        self._structure = None # reset the dot bracket
        self._trigger_callbacks(**kwargs)
        
    def _calculate_basepair(self):
        """
        Calculate a dictionary with positions as keys and basepairs symbols as value and store it in the _basepair attribute.
        The basepair is calculate simply considering complementary bases one position away! It's not very accurate!
        """
        basepair = BasePair(callback=self._updated_basepair)
        for ind1, strand1 in enumerate(self): # take first strand
            bmap1 = strand1.base_map # get the nucleotide map of the first strand
            for strand2 in self[ind1:]: # take second strand (use 'self[ind1:]' to avoid strands already calculated)
                bmap2 = strand2.base_map # get the nucleotide map of the second strand
                for pos1, base1 in bmap1.items(): # check each base in the first strand
                    if pos1 in basepair: 
                        continue # if the base is already occupied in a base pairing, skip it
                    for direction in ((2,0), (0,2), (-2,0), (0,-2)):
                        # makig a for cycle to control if there is a complementary base in all the direction  
                        pos2 = (pos1[0] + direction[0], pos1[1] + direction[1]) # calculate the position of the second base
                        if pos2 in basepair: continue # if the position is already occupied, skip it
                        bp_pos = (pos1[0] + direction[0]//2, pos1[1] + direction[1]//2) # calculate the position of the basepair symbol
                        if pos2 in bmap2 and bp_pos not in strand2.map and bp_pos not in strand1.map: 
                            # the strand has a base in the pairing position and the basepair symbol position is free
                            base2 = bmap2[pos2]
                            if base2 in base_pairing[base1]: # check the base pairing works
                                # append position and connected symbol 
                                basepair[Position(pos2)] = Position(pos1)
                                break # don't test other directions
        # update the basepair
        self._basepair = basepair

    def _check_addition(self, other, direction = (1,0)):
            if not self or not other:
                return # if one of the two motifs is empty, return without error
            if not isinstance(other, Motif):
                raise ValueError(f'{other} is not a valid type for addition')
            # take the junctions of the left and right side of the motifs
            Strand._check_position(input_pos_dir = direction, direction=True)
            junction1 = self.junctions[direction]
            junction2 = other.junctions[(- direction[0], - direction[1])]
            if not junction1 or not junction2:
                raise MotifStructureError(f"The motifs cannot be added in the direction {direction}, the junctions are missing. Junctions motif1: {junction1}, Junctions motif2: {junction2}. If you want to concatenate the motifs, use pf.Motif.concat() method.")
            elif len(junction1) != len(junction2):
                raise MotifStructureError(f"The number of connecting strands is different in the direction {direction}. Junctions motif1: {junction1}, Junctions motif2: {junction2}")
            for i in range(len(junction1)-1):
                # calculate the y difference between strands junctions, it must match
                if abs(junction1[i][1] - junction1[i+1][1]) != abs(junction2[i][1] - junction2[i+1][1]):
                    raise MotifStructureError(f"The junctions have a different y distance; junctions1: {junction1}, junctions2: {junction2}")  
                # # calculate the x difference between strands junctions, it must match
                if abs(junction1[i][0] - junction1[i+1][0]) != abs(junction2[i][0] - junction2[i+1][0]):
                    raise MotifStructureError(f"The junctions have a different x distance; junctions1: {junction1}, junctions2: {junction2}")
    
    @staticmethod
    def copy_strands_preserve_blocks(strands, motif = None):
        # IMPORTANT: keep strand that are part of the same motif linked
        motifs_id = {id(s.strands_block) for s in strands} # get the set of strands block id
        new_motifs_dict = {key: StrandsBlock() for key in motifs_id} # for each strands block, make a new one and link to the old id
        
        strands_copy = []
        # copy all strands
        for strand in strands:
            # check that the strands are in the strand class
            if not isinstance(strand, Strand):
                raise ValueError(f"{strand} is not a Strand object.")
            new_strands_block = new_motifs_dict[id(strand.strands_block)] # collect the new strands block
            # copy the strand and add the callback
            if motif:
                copied = strand.copy(callback = motif._updated_strands)
            else:
                copied = strand.copy()
            strands_copy.append(copied)
            new_strands_block.append(copied) # add the origami to the new strands block
        return strands_copy
    
    def save_text(self, filename: str = 'motif'):
        """ Save the motif in a text file. """
        with open(f'{filename}.txt', 'w') as f:
            f.write(str(self))

    def save_structure_text(self, filename: str = 'motif'):
        """ Save the structure of the motif in a text file. """
        # split the structure for each strand
        splitted_db = self.structure.split('&')
        # create a new motif with the structure
        db_motif = self.copy()

        # in each strand, replace the bases with the structure symbols
        for s, db in zip(db_motif, splitted_db):
            seq_count = 0
            new_strand = list(s.strand)
            for i in range(len(s)):
                if new_strand[i] in iupac_code:
                    new_strand[i] = db[seq_count]
                    seq_count += 1
            # update the strand with the new structure
            s._strand = ''.join(new_strand)
            s._map = None
        db_motif.save_text(filename)