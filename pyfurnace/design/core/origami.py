from functools import wraps
from typing import List, Tuple, Union, Literal, Callable
from .symbols import *
from .position import Position, Direction
from .callback import Callback
from .sequence import Sequence
from .strand import Strand
from .motif import Motif


class Origami(Callback):
    """
    A class for building and manipulating RNA Origami structures.

    The Origami class organizes RNA `Motif` objects into a 2D matrix to
    represent spatial arrangements of strands. It supports stacking motifs
    horizontally and vertically, calculating vertical junctions and connections,
    and exporting the full structure and sequence.

    Parameters
    ----------
    matrix : Motif or list of Motif or list of list of Motif
        A motif or 1D/2D list of Motifs to initialize the origami layout.
    *args : Motif or list of Motif
        Additional motifs or rows of motifs to add to the matrix.
    align : {'left', 'first', 'center'}, default='left'
        How motifs should be aligned vertically.
    copy : bool, default=False
        Whether to create a copy of the motifs before adding them to the origami.
    ss_assembly : bool, default=False
        Wether to assemble the 3D structure of the origami without locking the 
        coordinates of the motifs.
    **kwargs : dict
        Additional keyword arguments passed to the Callback base class.

    Attributes
    ----------
    motif : Motif
        Combined Motif representing the full assembled Origami.
    sequence : str
        Full nucleotide sequence of the Origami.
    structure : str
        Dot-bracket notation of the RNA secondary structure.
    strands : list of Strand
        List of individual strands in the Origami.
    junctions : list of dict
        Connection junctions (left, right, top, bottom) for each row.
    base_map : dict
        Map from nucleotide symbol positions (x, y) to matrix slice (y, x).
    map : dict
        Map from symbol positions (x, y) to motif matrix slice (y, x).
    shift_map : dict
        Map from motif matrix slices (y, x) to spatial shifts (x, y) of 
        the motifs to shift them from their original position to the
        assembled position.
    motif_map : dict
        Map from nucleotide symbol positions (x, y) to the corresponding Motif.
    sequence_index_map : dict
        Map from nucleodited symbol positions (x, y) to sequence index.
    pseudoknots : list of dict
        Parsed pseudoknot metadata including indices and energetics.
    num_motifs : int
        Total number of motifs in the Origami.
    num_lines : int
        Number of horizontal lines (rows) in the Origami.
    num_char : list of int
        Number of characters per line, used for alignment.
    align : str
        Current vertical alignment mode.
    ss_assembly : bool
        Wether to assemble the 3D structure of the origami without locking the 
        coordinates of the motifs.

    See Also
    --------
    Motif, Strand, Sequence
    """


    def __init__(self, 
                 matrix: Union[Motif, List[Motif], List[List[Motif]]] = None, 
                 *args: Union[Motif, List[Motif]], 
                 align: Literal['left', 'first', 'center'] = 'left', 
                 copy: bool = False, 
                 ss_assembly: bool = False, 
                 **kwargs) -> None:
        """
        Initialize an Origami object with a 2D list of motifs.

        Parameters
        ----------
        matrix : Motif or list of Motif or list of list of Motif
            A motif or 1D/2D list of Motifs to initialize the origami layout.
        *args : Motif or list of Motif
            Additional motifs or rows of motifs to add to the matrix.
        align : {'left', 'first', 'center'}, default='left'
            How motifs should be aligned vertically.
        copy : bool, default=False
            Whether to create a copy of the motifs before adding them to the origami.
        ss_assembly : bool, default=False
            Wether to assemble the 3D structure of the origami without locking the 
            coordinates of the motifs.
        **kwargs : dict
            Additional keyword arguments passed to the Callback base class.
        """
        # initialize the callback
        Callback.__init__(self, **kwargs)

        # initialize the protected atrributes
        self._junctions = None
        self._motif = None
        self._map = None
        self._base_map = None
        self._shift_map = None
        self._assembled = None
        self._ss_assembly = bool(ss_assembly)
        self._pseudoknots = None

        # initialize the matrix
        if not matrix:
            matrix = []

        ### CHECK THE MATRXI
        # the matrix is a proper 2D list
        if (isinstance(matrix, (list, tuple)) 
                and all(isinstance(row, (list, tuple)) for row in matrix)):
            pass

        # the matrix is a list of motif
        elif (isinstance(matrix, (list, tuple)) 
                and any(isinstance(row, Motif) for row in matrix)):
            matrix = [matrix] # make it a 2D list

        # the matrix is a motif
        elif isinstance(matrix, Motif):
            matrix = [[matrix]]

        # unsupported type for matrix
        else:
            raise ValueError(f'The matrix variable may only contain lists of motifs'
                             f' or motifs, but it contains {type(matrix)}.')
        
        ### check the type of the args variable
        if args:
            # the args contanes lines (lists of motifs), so args it's a matrix
            if all(isinstance(item, (list, tuple)) for item in args): 
                # add the lines to the matrix
                matrix.extend(args)
                # the args contains only motifs, so it's a line
            elif all(isinstance(item, Motif) for item in args): 
                matrix[-1].extend(args) ### add the motifs to the last line
            else:
                raise ValueError(f'The args variable may only contain lists of'
                                 f' motifs or motifs, but it contains {type(args)}.')

        ### add the matrix to the object
        self._matrix = matrix
        # prepare the motifs
        self._prepare_matrix(copy=copy) 

        ### set the alignment type
        if align not in ('left', 'first', 'center'):
            raise ValueError(f'"{align}" is not an accepted value for the'
                              'align_type variable. The align_type variable'
                              ' can only be: "left", "first" or "center".')
        else:
            self._align = align

    ### 
    ### MAGIC METHODS
    ###
        
    def __add__(self, other: 'Origami') -> 'Origami':
        """
        Horizontally add another Origami to this Origami.

        Parameters
        ----------
        other : Origami
            The origami to stack horizontally.

        Returns
        -------
        Origami
            A new Origami object with horizontally concatenated motifs.
        """
        new_matrix = [[m.copy() for m in row] for row in self._matrix]

        if not isinstance(other, Origami):
            raise TypeError(f"Unsupported operand type(s) for +: "
                            "'Origami' and '{type(other).__name__}'")

        # add extra rows to the new matrix
        diff_len = len(other._matrix) - len(new_matrix)
        if diff_len > 0:
            new_matrix.extend([[] for _ in range(diff_len)])

        for i, row in enumerate(other._matrix):
            new_matrix[i].extend([m.copy() for m in row])

        return Origami(new_matrix, 
                       align=self.align, 
                       ss_assembly=self.ss_assembly, 
                       copy=False)

    
    def __bool__(self) -> bool:
        """ Return False there are no motifs or all motifs are empty"""
        if not self._matrix:
            return False
        for row in self:
            for motif in row:
                if motif:
                    return True
        return False
    
    def __getitem__(self, 
                    key: Union[int, 
                               slice, 
                               Tuple[int, int], 
                               Tuple[slice, slice],
                               Callable[[Motif], bool]]
                    ) -> Union[Motif, 
                               List[Motif], 
                               List[List[Motif]]]:

        """
        Get motifs from the matrix using slicing or filtering.

        Parameters
        ----------
        key : int, slice, tuple or callable
            If int or slice, returns the corresponding row(s) of motifs.
            If tuple of two ints, returns the motif at that position.
            If tuple of two slices, returns a submatrix of motifs.
            If a function, returns a submatrix of motifs that satisfy the function.

        Returns
        -------
        Union[Motif, List[Motif], List[List[Motif]]]
            Retrieved motif(s).

        Raises
        ------
        TypeError
            If the key is of unsupported type.
        """
        ### 2D slice
        if isinstance(key, (tuple, list)):
            y, x = key

            # two slices, return a 2D list, a sub-origami
            if isinstance(x, slice) and isinstance(y, slice): 
                return [line[x] for line in self._matrix[y]]
            
            # two integers, return a single motif
            if all(isinstance(i, int) for i in key): 
                return self._matrix[y][x]
            
            # convert any index to a slice, then return a 2D list
            if isinstance(y, int):
                y = y % len(self._matrix)
                y = slice(y, y + 1)
            if isinstance(x, int):
                x = x % len(self._matrix[y])
                x = slice(x, x+1)
            return [m for row in self._matrix[y] for m in row[x]]
        
        ### 1D slice, return the row
        elif isinstance(key, slice) or isinstance(key, int):
            return self._matrix[key]
        
        ### Function, return a 2D list of motifs that satisfy the function
        elif hasattr(key, '__call__'): 
            return [[m for m in row if key(m)] for row in self._matrix]
        
        else:
            raise TypeError("Index must be a single int/slice, "
                            "a tuple of (row, col) of int/slice, "
                            "or a function to screen the motifs, got: "
                            f"{key}, of type: {type(key)}"
                            )
    
    def __setitem__(self,
                    key: Union[int,
                            slice,
                            Tuple[int, int],
                            Tuple[slice, slice],
                            Callable[[Motif], bool]],
                    value: Union[Motif,
                                List[Motif],
                                List[List[Motif]]]) -> None:
        """
        Set motif(s) at the specified position in the matrix, trying to
        match the shape of the value to the shape of the key. The value
        is always copied when added to the matrix to avoid references
        problem when repeating motifs in the matrix.
        If the value is a single of motifs or a list of motifs, it will be 
        set in the selected row(s). If the value is a 2D list of motifs, it will be set in the
        selected region of the matrix only if the selected region is a 2D region.

        Parameters
        ----------
        key : int, slice, tuple, or callable
            If int or slice, sets the entire row(s).
            If tuple of two ints, sets a single motif at that position.
            If tuple of slices, sets a 2D region in the matrix.
            If a function, replaces motifs that satisfy the condition.

        value : Motif, list of Motif, or 2D list of Motif
            The motif(s) to insert. Must match the shape implied by `key`.

        Raises
        ------
        ValueError
            If the value does not match the expected dimensions 
            or contains invalid types.
        TypeError
            If the key is of unsupported type.
        """

        ### CHECK THE DIMENSIONALITY OF THE VALUE

        # the value is a single motif
        if isinstance(value, Motif):
            value_dimension = 0
            value = [value]

        # value is a 1D list of motifs
        elif (isinstance(value, list) and 
                all(isinstance(item, Motif) for item in value)):
            value_dimension = 1

        # value is a 2D list of motifs
        elif (isinstance(value, list)  \
                and all(isinstance(item, (list, tuple)) for item in value)
                and all(isinstance(m, Motif) for sublist in value for m in sublist)):
            value_dimension = 2

        else:
            raise ValueError(f'Only motifs, lists of motifs, or 2D lists '
                             f'of motifs can be added to the Origami, but '
                             f'the object {value} was added.')
        
        ### CHECK THE DIMENSIONALITY OF THE KEY

        mask = None
        # if the key is a function, we need to create a mask wich is a submatrix
        # with the slices of the motifs that satisfy the function
        if hasattr(key, '__call__'): 
            mask = [[(i, slice(j, j+1)) 
                        for j, m in enumerate(row) if key(m)]
                            for i, row in enumerate(self._matrix)]
            
        # the key is a single int, we select a row
        elif isinstance(key, int):
            key_dimension = 1
            y_int = key % len(self._matrix)

            y_slice = slice(y_int, y_int + 1)
            x_slice = slice(0, len(self._matrix[key]))

        # the key is a tuples of ints, we select a single motif
        elif (isinstance(key, (tuple, list))
                and all(isinstance(i, int) for i in key)):
            key_dimension = 0
            # Convert the keys to a positive integer
            y_int = key[0] % len(self._matrix)
            x_int = key[1] % len(self._matrix[y_int])

            # Convert the keys to slices
            y_slice = slice(y_int, y_int + 1)
            x_slice = slice(x_int, x_int + 1)

        # the key is a slice of a row and slice of a column
        # so this is still a 1D region
        elif isinstance(key, (tuple, list)) and\
                isinstance(key[0], int) and isinstance(key[1], slice):
            key_dimension = 1 # select a row
            y_int = key[0] % len(self._matrix)

            y_slice = slice(y_int, y_int + 1)
            x_slice = key[1] # get the slice

        # key selects a submatrix, so this is a 2D region
        elif isinstance(key, (tuple, list)) and\
                all(isinstance(i, slice) for i in key):
            key_dimension = 2 # select a 2D region
            y_slice = key[0]
            x_slice = key[1]

        # special case: vertical selection (theoretically a 1D region,
        # but for code purposes we need to treat it as a 2D region)
        elif isinstance(key, (tuple, list)) and\
                isinstance(key[0], slice) and isinstance(key[1], int):
            key_dimension = 2 # still select a 2D region, but VERTICAL
            x_int = key[1] % len(self._matrix[key[0]])
            
            y_slice = key[0]
            x_slice = slice(x_int, x_int + 1)

        else:
            raise TypeError("Origami indexes can be: \n"\
                            "\t - a function to screen the motifs, \n"\
                            "\t - an int/slice to select a row, \n"\
                            "\t - a tuple of two int/slice to select a region. \n"\
                            f"Got: {key}, of type: {type(key)}"
                            )

        ### APPLY THE MASK (if any)
        if mask is not None:
            for row_ind, row in enumerate(mask):
                for mot_ind, (i, j) in enumerate(row):

                    # just put in position for 0D/1D values
                    if value_dimension in (0, 1):
                        self._matrix[i][j] = [m.copy() for m in value]

                    # dimensionality 2: match the two submatrices
                    elif value_dimension == 2:

                        # try to set the value matching the indices
                        try:
                            self._matrix[i][j] = [value[row_ind][mot_ind]].copy()
                        except IndexError as e:
                            raise IndexError(
                                f'Error while setting the value to the Origami. '
                                f'The lists do not match. Origami indexes: y: {i}, '
                                f'x: {j}.') from e

        ### REDUCED ALL CASES TO A 2D SLICING
        else:
            for i, line in enumerate(self._matrix[y_slice]):

                # dimensionality 0 or 1: just set the value
                if value_dimension in (0, 1):
                    line[x_slice] = [m.copy() for m in value]

                # dimensionality 2 math the two submatrices
                elif key_dimension ==  2 and value_dimension == 2:
                    # try to set the value matching the indices
                    try:
                        line[x_slice] = [m.copy() for m in value[i]]
                    except IndexError as e:
                        raise IndexError(
                            f'Error while setting the value to the Origami. '
                            f'The lists do not match. Origami indexes: '
                            f'y index: {i}, x slice: {x_slice}.') from e

        ### update the motif
        self._prepare_matrix() # this also prepare the motifs
        self._updated_motif()

    def __len__(self):
        """ Get the number of rows in the origami"""
        return len(self._matrix)
    
    def __str__(self) -> str:
        """ Return a string representation of the assmebled origami 
        (the origami motif). """
        return str(self.motif)

    def __repr__(self):
        """ Return a string representation of the origami object, by
        iterating through the matrix and calling the repr method of
        each motif. """
        reprs = ''
        for line in self._matrix:
            for item in line:
                reprs += repr(item) + ', '
            reprs += ';\n'
        return reprs
    
    ###
    ### PROPERTIES
    ###

    @property
    def num_motifs(self) -> int:
        """ The number of motifs in the origami. """
        return sum(1 for line in self._matrix for item in line if item is not None)
    
    @property
    def num_char(self) -> List[int]:
        """ The number of characters in each line of the origami. """
        if not self._matrix:
            return 0
        return [sum(m.num_char for m in line) for line in self._matrix]

    @property
    def num_lines(self) -> int:
        """ The number of lines in the origami. """
        if not self._matrix:
            return 0
        return len(self._matrix)
    
    @property
    def ss_assembly(self) -> bool:
        """ Boolean indicating if the origami 3d structure
        is assembled without locking the coordinates of the motifs. """
        return bool(self._ss_assembly)
    
    @ss_assembly.setter
    def ss_assembly(self, new_ss_assembly):
        """ Set the ss_assembly attribute to True or False. """
        self._ss_assembly = bool(new_ss_assembly)
        self._updated_motif()
    
    @property
    def junctions(self) -> List[dict]:
        """ The vertical junctions for each line of the Origami.
        The junctions are a dictionary with directions as keys and positions of the 
        junctions as values.
        The motif junctions are extended, so we can just take the maximum position 
        of the motifs for connections. """

        ### check if the junctions are already calculated
        if self._junctions:
            return self._junctions
        
        junctions_list = [] # initialize the list

        ### Calculate the junctions for each line
        for line in self._matrix:
            x_shift = 0

            # initialize the dictionary
            line_junc_dict = {direct: [] for direct in Direction}
            
            # iterate through the motifs in the line
            for ind, m in enumerate(line):

                ### NOT USED
                # # add the left junctions of the first motif
                # if ind == 0: 
                #     line_junc_dict[Direction.LEFT] = m.junctions[Direction.LEFT]

                shift = Position((x_shift, 0))
                ### add the top and bottom junctions
                line_junc_dict[Direction.UP] += [shift + p 
                                                 for p in m.junctions[Direction.UP]]
                line_junc_dict[Direction.DOWN] += [shift + p 
                                                  for p in m.junctions[Direction.DOWN]]

                # just add numchar becasue the junctions are extended
                x_shift += m.num_char 

            ### NOT USED
            # # if it's the last motif, add the right junctions
            # line_junc_dict[(1, 0)] += [x_shift + p for p in m.junctions[(1, 0)]]

            # add the junctions to the list
            junctions_list.append(line_junc_dict)

        # save the junctions
        self._junctions = junctions_list
        return junctions_list

    @property
    def align(self) -> Literal['left', 'first', 'center']:
        """ The alignment type of the rows of the origami. """
        return self._align
    
    @align.setter
    def align(self, new_align):
        """ Set the alignment type of the rows of the origami. 

        Parameters
        ----------
        new_align : {'left', 'first', 'center'}
            The new alignment type for the origami.
            When set to 'left', the motif rows are aligned to the left.
            When set to 'first', the motifs rows are aligned to match the first
            vertical junction.
            When set to 'center', the motifs rows are aligned to the center.

        """

        if new_align not in ('left', 'first', 'center'):
            raise ValueError(f'"{new_align}" is not an accepted value for '
                             'the align_type variable. The align_type variable '
                             'may only a string reading "left", "first" or "center".')
        self._align = new_align
        self._updated_motif()

    @property
    def assembled(self):
        """ The matrix of the origami with the motif shifted in the correct position 
        for the assembly. The assembled matrix contains rows with the vertical 
        connection motifs. """
        if not self._assembled:
            self._assemble()
        return self._assembled
    
    @property
    def map(self) -> dict:
        """ A dictionary with the symbols position (x, y) as keys and the
        matrix index (y, x) of the motif that contains it as values. """
        if self._map is not None:
            return self._map
        
        ### build the map and the motif map
        assembled = self.assembled
        self._map = {}
        self._motif_map = {}
        for y, line in enumerate(assembled):
            # skip the connection motifs
            if y % 2 == 1:
                continue

            for x, m in enumerate(line):
                for pos in m.map:
                    # y//2 because for each line there is a connection motif
                    self._map[pos] = (y//2, x)
                    self._motif_map[pos] = m

        return self._map
    
    @property
    def shift_map(self) -> dict:
        """ A dictionary with the slice of the motif in the matrix as key (y, x) 
        and positional shift of the motif as values (y, x). The shift is the 
        difference between the position of the motif in the matrix and the position
        of the motif in the assembled origami.
        """
        if self._shift_map:
            return self._shift_map
        # the shift_map is calculated when the origami is assembled
        self._assemble()
        return self._shift_map
    
    @property
    def motif_map(self) -> dict:
        """ A dictionary with the symbols position (x, y) as keys and the
        motif that contains it as values. """
        if not self._map:
            # the motif_map is calculated when the map is calculated
            self.map
        return self._motif_map
    
    @property
    def sequence_index_map(self) -> dict:
        """ A dictionary with the nucleotide position (x, y) as keys and the
        index of the nucleotide in the sequence as values. """
        if not self._motif:
            self.motif
        return self.motif.sequence_index_map
    
    @property
    def base_map(self) -> dict:
        """ A dictionary with the nucleotide position (x, y) as keys and the
        matrix index (y, x) of the motif that contains it as values. """
        if not self._base_map:
            self._base_map = {}
            origami_map = self.map
            origami_motif_base_map = self.motif.base_map
            for pos in origami_motif_base_map.keys():
                self._base_map[pos] = origami_map[pos]
        return self._base_map
    
    @property
    def sequence(self) -> 'Sequence':
        """ The sequence of the origami, as a Sequence. """
        return self.motif.sequence
    
    @sequence.setter
    def sequence(self, new_seq):
        """Set the sequence of the origami. """

        # remove the '&' symbol
        new_seq = new_seq.replace('&', '') 
        current_seq = self.sequence.replace('&', '') 

        if (not isinstance(new_seq, (str, Sequence)) 
                or len(new_seq) != len(current_seq)):
            raise ValueError(f"The new sequence must be a string or a Sequence object"
                             f" with the same lenght of the current sequence "
                             f"({len(current_seq)}). Got type: {type(new_seq)}; with "
                             f"length: {len(new_seq)}, excluding the '&' symbols.")
        
        # adjust the offset if there are multiple strands
        offset = 0 

        # read the maps once to avoid triggering the callback and origami assembly
        origami_motif = self.motif
        motif_map = self.map
        motif_shifts = self.shift_map

        # iterate over the strands in the origami motif
        for s in origami_motif:
            
            # a tuple to identify a specific strand in a motif in the 
            # origami matrix
            strand_ID = None
            # initialize/reset the current base map
            new_strand_seq = ''

            # iterate over the nucleotides in the strand
            for ind, pos in enumerate(s.base_map):
                
                ### GET THE STRAND ID FOR THIS NUCLEOTIDE

                # get the y, x cooridnates of the motif in the matrix
                motif_yx = motif_map[pos]
                # get the x, y shift of the motif in the origami positions
                shift_yx = motif_shifts[motif_yx]
                # remove the shifts from the position of the base
                original_pos = (pos[0] - shift_yx[0], pos[1] - shift_yx[1])
                # get the motif at the position
                motif = self._matrix[motif_yx[0]][motif_yx[1]]
                # get the strand index of the motif at the base position
                strand_ind = motif.map[original_pos]
                if strand_ID is None:
                    strand_ID = (motif_yx[0], motif_yx[1], strand_ind)

                ### NEW STRAND ID? 
                # update the sequence of the previous strand before 
                # moving to the next one
            
                if strand_ID != (motif_yx[0], motif_yx[1], strand_ind):
                    
                    # get the strand and set the curent base maps
                    strand = self._matrix[strand_ID[0]][strand_ID[1]][strand_ID[2]]
                    strand.sequence = new_strand_seq

                    # reset the current base map with the new strand
                    new_strand_seq = ''

                    # update the motif and strand position to the current position
                    strand_ID = (motif_yx[0], motif_yx[1], strand_ind)

                # add the new sequence to the current base map
                new_strand_seq += new_seq[ind + offset]

            # add the last strand
            last_strand = self._matrix[strand_ID[0]][strand_ID[1]][strand_ID[2]]
            last_strand.sequence = new_strand_seq

            # update the offset 
            offset += len(s.sequence)

    @property
    def pseudoknots(self) -> List[dict]:
        """ A list of dictionaries with the pseudoknot information.
        Each dictionary contains the following keys:
            - id: a pseudoknot index to identify the pseudoknot
            - ind_fwd: a list of tuples (start, end) with the indices of the forward 
                       sequences of the pseudoknot
            - ind_rev: a list of tuples (start, end) with the indices of the reverse 
                       sequences of the pseudoknot
            - E: the energy of the pseudoknot
            - dE: the energy tolerance of the pseudoknot
        """
        if self._pseudoknots:
            return self._pseudoknots
        
        # A dictionary to store the pseudoknot information, with the pk_index as key
        # and the pk information dict as value
        pk_dict = dict()
        
        def add_pk(strand, pk_index, info_nr, shift, forward=True):
            """ Add the pseudoknot information to the pk_dict. """
            # get the pseudoknot information
            pk_info = strand.pk_info

            # add the pseudoknot information to the pk_dict
            pk_dict.setdefault(pk_index, 
                               {"id": pk_index, 
                                'ind_fwd': [], 
                                'ind_rev': [], 
                                'E': [], 
                                'dE': []})
            pk_dict[pk_index]['E'].append(pk_info['E'][info_nr])
            pk_dict[pk_index]['dE'].append(pk_info['dE'][info_nr])

            # indicate the index of the pseudoknot in the sequence
            start_pos_ind = 0 
            if strand.directionality == '35':
                start_pos_ind = -1
            pos = list(strand.base_map.keys())[start_pos_ind]

            # get the index of the sequence in the strand
            offset_ind = self.sequence_index_map[(shift[0] + pos[0], 
                                                  shift[1] + pos[1])]
            
            # get the start and end positions of the pseudoknot
            pk_start, pk_end = pk_info['ind_fwd'][info_nr]
            start_end_tuple = (offset_ind + pk_start, offset_ind + pk_end)
            # add the start and end positions to the pk_dict
            if forward:
                pk_dict[pk_index]['ind_fwd'].append(start_end_tuple)
            else:
                pk_dict[pk_index]['ind_rev'].append(start_end_tuple)

        pk_motifs = []

        ### collect all the motifs with pseudoknot information
        for i, line in enumerate(self._matrix):
            for j, m in enumerate(line):
                if any(hasattr(s, 'pk_info') for s in m):
                    pk_motifs.append((i, j))
        
        ### Iterate through the strands of the motifs with pseudoknot information
        for i, j in pk_motifs:
            m = self._matrix[i][j]
            shift = self.shift_map[(i, j)]

            # get pseudoknot IDs from the strands
            pk_strands = [s for s in m if hasattr(s, 'pk_info')]
            pk_indexes = [pk_id for s in pk_strands for pk_id in s.pk_info['id']]

            ### Adjust the pk_index for unique pseudoknots
            if any(ind[0] == '0' for ind in pk_indexes): # new 0 pseudoknot
                current_n_zero = sum(1 if key[0] == '0' else 0 for key in pk_dict)
                pk_index_0 = '0_' + str(current_n_zero + 1)            

            # add the pseudoknots
            for strand in pk_strands:
                for info_nr, pk_index in enumerate(strand.pk_info['id']):
                    reverse = pk_index[-1] == "'"
                    if pk_index[0] == '0':
                        pk_index = pk_index_0
                    elif reverse:
                        pk_index = pk_index[:-1]
                    add_pk(strand, pk_index, info_nr, shift, forward=not reverse)

        # make the average energy and average tolerance
        for pk in pk_dict.values():
            pk['E'] = sum(pk['E']) / len(pk['E'])
            pk['dE'] = sum(pk['dE']) / len(pk['dE'])
        # convert the pk_dict to a list for simplicity
        self._pseudoknots = list(pk_dict.values())
        return self._pseudoknots

    @property
    def structure(self) -> str:
        """ The dot-bracket structure of the origami. """
        return self.motif.structure
    
    @property
    def pair_map(self) -> dict:
        """The dictionary of the paired indexes (alternative to the dot bracket 
        notation). """
        if not self._motif:
            self.motif
        return self.motif.pair_map
    
    @property
    def motif(self) -> Motif:
        """ The assembled origami motif. """
        if isinstance(self._motif, Motif):
            return self._motif
        mot = None
        for line in self.assembled:
            # concatenate the motifs in the line
            mot_line = Motif.concat(line, 
                                    align=False, 
                                    unlock_strands=self._ss_assembly, 
                                    lock_coords=False)
            # add the line to the motif
            mot = Motif.concat([mot, mot_line], 
                               axis=0, align=False, 
                               lock_coords=True, 
                               unlock_strands=self._ss_assembly)
        
        self._motif = mot
        return self._motif
    
    @property
    def strands(self) -> List[Strand]:
        """ The strands of the origami. """
        return list(self.motif)

    ### 
    ###  STATIC METHODS
    ###

    @staticmethod
    def _calculate_connections(junctions1: dict,
                               junctions2: dict,
                               x_shift: Union[tuple[int, int], Position]=(0, 0), 
                               start_y: int=0) -> Tuple[Motif, int]:
        """
        Creates the connection between the rows of the origami.

        Parameters
        ----------
        junctions1: dict
            junctions of the first line
        junctions2: dict
            junctions of the second line
        x_shift: tuple
            The x shift of the junctions of the first and second line
        start_y: int
            The y position of the first line

        Returns
        -------
        Tuple[Motif, int]
            The connection motifs and height of the vertical connections
        """
        ### take the junctions of the two lines
        j1 = [pos[0] + x_shift[0] for pos in junctions1[Direction.DOWN]] 
        j2 = [pos[0] + x_shift[1] for pos in junctions2[Direction.UP]]

        # a junction is missing, then no connection
        if not j2 or not j1:
            return Motif(), 0
        
        # the number of connections is the minimum of the two junctions
        n_connect =  min((len(j1), len(j2)))
        j1 = j1[:n_connect]
        j2 = j2[:n_connect]
        
        ### CREATE THE CONNECTIONS

        # a dictionary with the connented pair index as key and 
        # a set of crossed pair indexes as value
        closed_crossings = dict() 
        # the positions that should be connected
        pairs = list(zip(j1, j2)) 

        for ind, (x1, x2) in enumerate(pairs):
            # intialize the crossed pair indexes
            closed_crossings[ind] = set() 
            # the minimum x position to connect in this pair
            x_min = min(x1, x2) 
            # the maximum x position to connect in this pair
            x_max = max(x1, x2) 

            # the crossed pairs are pairs that have at list one position between
            # the minimum and maximum positions and are not already connected
            crossed = {i for i, x12 in enumerate(pairs) 
                            if (i not in closed_crossings 
                                and (x_min <= x12[0] <= x_max 
                                     or x_min <= x12[1] <= x_max))}
            
            # update the crossed pairs for this connection
            closed_crossings[ind].update(crossed)

        ### CHECK FOR NESTED CROSSINGS
        # if pair #1 crosses pair #2; then pair #2 crosses pair #3 and pair #4
        # the shift of pair #1 has to take into account also pair #3 and pair #4

        # go through the connected pairs
        for key1 in list(closed_crossings.keys()): 
            for key2 in list(closed_crossings.keys()): 
                # if the second pair is in the crossed pairs of the first pair
                if key2 in closed_crossings[key1]: 
                    closed_crossings[key1].update(closed_crossings[key2])

        # calculate the maximum number of crossings
        max_crossing = max(len(crossed) for crossed in closed_crossings.values())

        ### MAKE THE STRANDS
        strands = [] 
        for ind, (x1, x2) in enumerate(pairs): 
            # the the number of crossings for this pair
            n_crossings = len(closed_crossings[ind])

            if x1 < x2: # the first motif is on the left
                strand = (  "│" * n_crossings 
                          + "╰" 
                          + "─" * (x2 - x1 -1) 
                          + "╮" 
                          + "│" * (max_crossing - n_crossings))

            elif x1 > x2: # the first motif is on the right
                strand = (  "│" * (max_crossing - n_crossings) 
                          + "╯" 
                          + "─" * (x1 - x2 -1) 
                          + "╭" 
                          + "│" * n_crossings)

            else: # the motifs are on the same position vertically
                strand = "│" * (max_crossing + 1)

            # can add the symbol "^" for retrocompatibility with ROAD
            strand += "↑" # if you  do this, increase the max_crossing by 1
            strands.append(Strand(strand, 
                                  start=(x1, start_y), 
                                  direction=Direction.DOWN))

        # Extra +1 to the max_crossing to add the symbol "^" or "↑"
        connection_height = max_crossing + 1 + 1 
        return Motif(strands, join=False), connection_height
    
    
    ### 
    ### PROTECTED METHODS
    ###
    
    def _prepare_matrix(self, copy: bool = False) -> None:
        """ Remove the empty motifs, add the callbackss and extend the junctions 
        (make sure all the bottom junctions reach the same y position).
        
        copy : bool, optional
            If True, the motifs are copied before being added to the origami.
        """

        ### clean the empty motifs
        prepare_matrix = []
        for y, line in enumerate(self._matrix):
            new_line = []

            for x, m in enumerate(line):
                # check if the motif is empty
                if not m:
                    continue
                    
                # check the type
                if not isinstance(m, Motif):  
                    raise ValueError(f'Only Motif or subclass of Motif can be added to '
                                     f'the Origami, but the object type {type(m)} '
                                     f'was added at position x: {x}, y: {y}.')
                
                # Copy, add the callback and save it
                if copy:
                    m = m.copy()
                m.register_callback(self._updated_motif)
                new_line.append(m)

            if new_line:
                prepare_matrix.append(new_line)

        # save the updated matrix
        self._matrix = prepare_matrix

        for ind, line in enumerate(self._matrix):
            # Align the motifs
            aligned = Motif.align(line)

            # Extend the junctions
            max_y = max([m.max_pos[1] for m in aligned], default=0)
            for m in aligned:
                m.extend_junctions(until= (None, max_y))
            
            # Add the aligned motifs to the matrix
            self._matrix[ind] = aligned

    def _assemble(self):
        """
        Return a matrix where the motifs and the jucntions of the Origami are in place
        """
        self._prepare_matrix()

        ### check the alignment shifts
        junctions = self.junctions
        ### No alignment
        if self._align == 'left':
            align_list = [0] * len(self._matrix)

        ### Align the first junctions
        elif self._align == "first":
            # take the maximum x position of the first bottom junction in all lines
            top_junctions = [0] * len(self._matrix) #initialize all the shift to 0
            bot_junctions = [0] * len(self._matrix)
            align_list = [0] * len(self._matrix)
            for ind, line in enumerate(junctions): # go through the junctions of all lines
                if line[(0, -1)]: # take the x position of the top junction
                    top_junctions[ind] = line[(0, -1)][0][0]
                if line[(0, 1)]: # take the x position of the bottom junction
                    bot_junctions[ind] = line[(0, 1)][0][0]
            # calculate the alignment: shift to the right to align the first junctions
            for ind1, bot_junct in enumerate(bot_junctions[:-1]):
                top_junct = top_junctions[ind1 + 1]
                shift = align_list[ind1] + bot_junct - top_junct # calculate the shift to align the first junctions
                if shift > 0: # the bottom line is shifted to the right, all good
                    align_list[ind1 + 1] += shift
                    continue
                for ind2 in range(ind1 + 1): # the bottom line is shifted to the left
                    align_list[ind2] -= shift # shift all top lines the right to avoid reaching negative x positions
        
        ### Align the center of the motifs
        elif self.align == "center":
            # take the maximum center position of the motifs in all lines
            max_center = max([num_char // 2 for num_char in self.num_char])
            # calculate the alignment: shift to the left to align the center of the motifs
            align_list = [max_center - num_char // 2 for num_char in self.num_char]

        ### create the assembly representation of the origami
        vert_shift = 0 # initialize the vertical shift
        assembled_matrix = []
        self._shift_map = {}
        for ind, line in enumerate(self._matrix):
            # if not line:
            #     continue # skip empty lines
            # get the x shift to align the motifs
            sequential_shift = Motif.get_sequential_shift(line, position_based=False)
            # create the shifts of the motifs by considering the shift to align the helices, the shift to place the motifs in the line, and the vertical shift
            shifts = [(align_list[ind] + shift, vert_shift) for shift in sequential_shift]
            for x, _ in enumerate(line):
                self._shift_map[(ind, x)] = shifts[x]
            # create the assembled line, copying and shifting the motifs
            assembled_matrix.append([m.copy().shift(shift) for m, shift in zip(line, shifts)])
            # update the vertical shift. If there is not vertical shift, take the last vert_shift
            vert_shift = vert_shift + max([m.num_lines for m in line], default=0)
            # create the connections between the motifs
            if ind < len(self._matrix) - 1:
               # the shift of the junctions due to the alignment
               junctions_shift = (align_list[ind], align_list[ind+1]) 
               # create the connection motif
               connection_motif, connection_height = self._calculate_connections(junctions[ind], junctions[ind + 1], x_shift = junctions_shift, start_y = vert_shift)
               # add the connection motif to the assembled matrix
               assembled_matrix.append([connection_motif])
               # update the vertical shift
               vert_shift += connection_height

        ### save the assembled matrix
        self._assembled = assembled_matrix
        return assembled_matrix

    ### 
    ### METHODS
    ###

    def index(self, condition):
        """ Find the motifs that satisfy the condition """
        if isinstance(condition, Motif):
            motif = condition
            condition = lambda x: x == motif
        elif not hasattr(condition, '__call__'): # if we submit a function, return the matrix filtered by the function
            raise ValueError(f'The condition must be a function or a Motif object, but {condition} was given.')
        return [(y, x) for y, line in enumerate(self._matrix) for x, m in enumerate(line) if condition(m)]

    @wraps(Motif.save_3d_model)
    def save_3d_model(self, *args, **kwargs):
        return self.motif.save_3d_model(*args, **kwargs)
    
    @wraps(Motif.folding_barriers) # inherit the documentation from the function
    def folding_barriers(self, kl_delay: int = 150) -> Tuple[str, int]:
        return self.motif.folding_barriers(kl_delay=kl_delay)

    def append(self, item, copy=False):
        """ Append a motif to the last line of the Origami"""
        if isinstance(item, Motif):
            if not self._matrix:
                self._matrix.append([])
            self._matrix[-1].append(item)
        elif isinstance(item, (list, tuple)) and all(isinstance(m, Motif) for m in item):
            self._matrix.append(item)
        else:
            raise ValueError(f'Only motifs or lists of motifs can be added to the Origami, but the object {item} was added.')
        self._prepare_matrix(copy=copy) # remove empty motifs and add callbacks
        self._updated_motif()
           
    def insert(self, idx, item, copy=False):
        """
        Inserts a helix at the given index
        ----------------------------------------------
        index: int or slice or tuple
            index at which the helix should be inserted, either a single index, a slice, or a tuple of slices
        item: Motif or list of Motifs
            Motif or list of motif that should be inserted
        """
        ### check the item variable
        if isinstance(item, (list, tuple)) and all(isinstance(m, Motif) for m in item):
            dimension = 2
        elif isinstance(item, Motif):
            dimension = 1
        else:
            raise ValueError(f'Only motifs or lists of motifs can be added to the Origami, but the object type {type(item)} was added.')
        
        ### add the item to the matrix according to the index
        if isinstance(idx, (int, slice)):
            if dimension == 1:
                self._matrix.insert(idx, [item])
            elif dimension == 2:
                self._matrix.insert(idx, item)
                
        elif isinstance(idx, (tuple, list)) and len(idx) == 2:
                if dimension == 2:
                    # raise ValueError(f'The index is the coordinate of a single motif: {idx}, but a list of motifs was added: {motifs}.')
                    for i, m in enumerate(item):
                        self._matrix[idx[0]].insert(idx[1] + i, m)
                if dimension == 1:
                    self._matrix[idx[0]].insert(idx[1], item)
                
        else:
            raise ValueError(f'Index must be a single index or a tuple of (row, col) or a list of two indices.')
        
        # remove empty motifs and add callbacks
        self._prepare_matrix(copy=copy)
        self._updated_motif()

    def pop(self, idx):
        """
        Removes the last line or the selected Motif
        ----------------------------------------------
        index: int, slice or tuple
            index at which the line or motif should be removed
        """
        if isinstance(idx, (int, slice)):
            if self._matrix:
                popped = self._matrix.pop(idx)
            else: return
        elif isinstance(idx, (tuple, list)) and len(idx) == 2:
            if self._matrix[idx[0]]:
                popped = self._matrix[idx[0]].pop(idx[1])
            else: return
        else:
            raise ValueError(f'Index must be a single index or a tuple of (row, col) or a list of two indices. Got {idx} instead.')
        self._prepare_matrix()
        self._updated_motif()
        return popped

    def remove(self, motif):
        """
        Removes the given motif from the Origami
        """
        for line in self._matrix:
            if motif in line:
                line.remove(motif)
                break
        self._prepare_matrix()
        self._updated_motif()

    def get_motif_at_position(self, position):
        """ Select a motif at the given position """
        if not isinstance(position, (tuple, list)) or len(position) != 2 or not all(isinstance(i, int) for i in position):
            raise ValueError(f'The position must be a tuple of two integers, but {position} was given.')
        position = tuple(position)
        if position not in self.map:
            raise ValueError(f'The position {position} is not in the map of the origami.')
        return self[self.map[(position)]]
    
    def get_motif_at_seq_index(self, index):
        """ Select a motif at the given sequence index """
        return self[self.get_slice_at_seq_index(index)]
    
    def get_slice_at_seq_index(self, index):
        """ Select a slice at the given sequence index """
        if not isinstance(index, int):
            raise ValueError(f'The sequence index must be an integer, but {index} (type: {type(index)}) was given.')
        ind_pos = None
        for pos, ind in self.sequence_index_map.items():
            if ind == index:
                ind_pos = pos
        if ind_pos is None:
            raise ValueError(f'The index {index} is not in the sequence index map of the origami.')
        return self.map[ind_pos]
    
    def get_strand_at_position(self, position):
        """ Select a strand at the given position """
        if not isinstance(position, (tuple, list)) or len(position) != 2 or not all(isinstance(i, int) for i in position):
            raise ValueError(f'The position must be a tuple of two integers, but {position} was given.')
        position = tuple(position)
        if position not in self.map:
            raise ValueError(f'The position {position} is not in the map of the origami.')
        moti_yx = self.map[position]
        motif = self[moti_yx]
        shift_motif = self.shift_map[moti_yx]
        corrected_pos = (position[0] - shift_motif[0], position[1] - shift_motif[1])
        strand_ind = motif.map[corrected_pos]
        return motif[strand_ind]
    
    def change_strand_at_position(self, position, new_strand, start=None, direction=None):
        """ Change the strand at the given position """
        strand = self.select_strand_at_position(position)
        if start:
            strand._start = start
        if direction:
            strand._direction = direction
        strand.strand = new_strand
        self._updated_motif()

    def duplicate_line(self, idx, insert_idx=None):
        """
        Duplicate a line in the origami
        ----------------------------------------------
        index: int
            index of the line that should be duplicated
        insert_idx: int
            index at which the duplicated line should be inserted
        copy_kissing_loops_seq: bool
            if True, the sequences of the kissing loops are copied to the new motifs
        """
        if not isinstance(idx, int):
            raise ValueError(f'The index must be an integer, but {idx} was given.')
        line = self._matrix[idx]
        new_line = [m.copy(callback=self._updated_motif) for m in line]
        if insert_idx is None:
            insert_idx = len(self._matrix)
        self._matrix.insert(insert_idx, new_line)
        self._prepare_matrix()
        self._updated_motif()

    def get_motif_type(self, motif_type):
        """ Return a list of motifs of the given type """
        return [m for line in self._matrix for m in line if type(m) == motif_type]
    
    def improve_folding_pathway(self, kl_delay: int = 150):
        from ..motifs import Stem
        from ..utils import start_end_stem

        # remove the motif that start the Origami
        origami = self.copy()
        start_ind = origami.index(lambda m: '5' in m)
        if start_ind:
            origami.pop(start_ind[0])

        # calculate the folding barriers
        start_barrier = origami.folding_barriers(kl_delay=kl_delay)[1]

        # Check the folding barriers of starting in each possible stem
        # of at least 5 bases
        db, stacks = dot_bracket_to_stacks(origami.structure)
        min_bar = start_barrier
        best_middle = 0
        for db, (start, end) in zip(db, stacks):
            if db in '()' and (end - start) > 4:
                middle = (start + end) // 2
                new_strucutre = rotate_dot_bracket(origami.structure, middle)
                new_bar = folding_barriers(kl_delay=kl_delay, structure=new_strucutre)[1]
                if new_bar < min_bar:
                    min_bar = new_bar
                    best_middle = middle
                    best_structure = new_strucutre

        # replace the starting motif with the new one
        for flip in range(2):
            ori_copy = origami.copy()
            start_slice = origami.get_slice_at_seq_index(best_middle)
            m = origami.get_motif_at_seq_index(best_middle)
            stem_1 = Stem(m.length // 2)
            start_end = start_end_stem()
            if flip:
                start_end.flip()
            stem2 = Stem(m.length - stem_1.length)
            ori_copy[start_slice] = [stem_1, start_end, stem2]

            if ori_copy.folding_barriers(kl_delay=kl_delay)[1] == min_bar:
                origami = ori_copy
                break

        return origami


    ### IMPLEMENT THIS!!!
    def find_best_start(self):
        return
        ### make a strand map of sequence index --> motif position
    
    def _updated_motif(self, **kwargs):
        self._motif = None
        self._junctions = None
        self._map = None
        self._base_map = None
        self._shift_map = None
        self._assembled = None
        self._pseudoknots = None
        self._trigger_callbacks(**kwargs)
        # self._trigger_callbacks()

    def to_road(self):
        ori_str = str(self)
        ori_str = ori_str.replace('↑', '^')
        ori_str = ori_str.replace('↓', '^')
        ori_str = ori_str.replace('│ ┊┊┊┊┊┊ │', '│ ****** │')
        ori_str = ori_str.replace(' ┊┊ ', ' !! ')
        ori_str = ori_str.replace(' ┊ ', ' ! ')
        return ori_str
    
    def barrier_repr(self, kl_delay: int = 150, barriers=None, return_list=False):
        motif = self.motif
        origami_lines = str(self).split('\n')
        if barriers is None:
            barriers = motif.folding_barriers(kl_delay=kl_delay)[0]
        for i, (x, y) in enumerate(motif.base_map):
            origami_lines[y] = origami_lines[y][:x] + barriers[i] + origami_lines[y][x+1:]
        if return_list:
            return origami_lines
        return '\n'.join(origami_lines)

    
    def save_text(self, filename_path, to_road=False):
        path = filename_path.split('.')[0]
        name = path.split('/')[-1].split('\\')[-1]
        with open(path + '.txt', 'w') as f:
            f.write(f'>{name}\n')
            f.write(str(self.sequence) + '\n')
            f.write(self.structure + '\n\n')
            if to_road:
                f.write(self.to_road())
            else:
                f.write(str(self))

    def save_structure_text(self, filename_path):
        self.motif.save_structure_text(filename_path)

    def copy(self):
        return Origami(self._matrix, ss_assembly=self.ss_assembly, align=self.align, copy=True)
    
    def reload(self):
        self._updated_motif()
        self._assemble()
