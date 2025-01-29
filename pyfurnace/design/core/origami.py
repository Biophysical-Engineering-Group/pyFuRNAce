from .symbols import *
from .callback import Callback
from .sequence import Sequence
from .strand import Strand
from .motif import Motif

class Origami(Callback):
    """
    Class that stacks Motifs vertically and horizontally to form the complete structure of the RNA Origami
    """

    def __init__(self, matrix: list = None, *args, align = 'left', copy=False, ss_assembly=False, **kwargs):
        """
        Attributes of class Origami
        ---------------------------------------------
        matrix: list
            2D list of all motifs making up the Origami
        """
        # initialize the callback
        Callback.__init__(self, **kwargs)

        # initialize the protected atrributes
        self._helices = None
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
        new_matrix = []

        ### check the matrix variable
        # the matrix is a proper 2D list
        if isinstance(matrix, (list, tuple)) and all(isinstance(row, (list, tuple)) for row in matrix):
            pass
        # the matrix is a list of motif
        elif isinstance(matrix, (list, tuple)) and any(isinstance(row, Motif) for row in matrix):
            matrix = [matrix] # make it a 2D list
        # the matrix is a motif
        elif isinstance(matrix, Motif):
            matrix = [[matrix]]
        # unsupported type for matrix
        else:
            raise ValueError(f'The matrix variable may only contain lists of motifs or motifs, but it contains {type(matrix)}.')
        
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
                raise ValueError(f'The args variable may only contain lists of motifs or motifs, but it contains {type(args)}.')
        
        ### prepare the motifs of the matrix
        for y, line in enumerate(matrix):
            new_line = []
            for x, m in enumerate(line):
                new_line.append(self._prepare_motif(m, y, x, copy=copy)) # prepare the motif
            new_matrix.append(new_line)

        ### add the matrix to the object
        self._matrix = new_matrix
        self._uniform_lines_len()

        ### set the alignment type
        if align not in ('left', 'first', 'center'):
            raise ValueError(f'"{align}" is not an accepted value for the align_type variable. The align_type variable may only a string reading "left", "first" or "center".')
        else:
            self._align = align

    ### 
    ### MAGIC METHODS
    ###
        
    def __add__(self, other):
        ''' 
        Addition of two elements of the class Origami or an element of class helix to an element of class origami
        --------------------------------------------------------------------------------------
        other: Origami
            element which is added (if it is of another type error occurs)
        --------------------------------------------------------------------------------------
        return: new Origami
        '''
        pass
    
    def __bool__(self):
        """ Return False there are no motifs or all motifs are empty"""
        if not self._matrix:
            return False
        for row in self:
            for motif in row:
                if motif:
                    return True
        return False
    
    def __getitem__(self, key):
        """ Get the block at index """
        ### 2D slice
        if isinstance(key, (tuple, list)):
            y, x = key
            ### if we submit two slices, return a 2D list, a sub-origami
            if isinstance(x, slice) and isinstance(y, slice): 
                return [line[x] for line in self._matrix[y]]
            if all(isinstance(i, int) for i in key): # if we submit two integers, return a single motif
                return self._matrix[y][x]
            if isinstance(y, int):
                y = y % len(self._matrix) # convert the index to a positive integer
                y = slice(y, y + 1)
            if isinstance(x, int):
                x = x % len(self._matrix[y])
                x = slice(x, x+1)
            return [m for row in self._matrix[y] for m in row[x]]
        ### 1D slice
        elif isinstance(key, slice) or isinstance(key, int):
            return self._matrix[key]
        elif hasattr(key, '__call__'): # if we submit a function, return the matrix filtered by the function
            # Return a Matrix with the motifs that satisfy the function
            return [[m for m in row if key(m)] for row in self._matrix]
        else:
            raise TypeError("Index must be a tuple of (row, col) or a int/slice to get a row.")
    
    def __setitem__(self, key, value):

        ### CHECK THE DIMENSIONALITY OF THE VALUE
        # the value is a single motif
        if isinstance(value, Motif):
            value_dimension = 0
            value = [value]
        # value is a 1D list of motifs
        elif isinstance(value, list) and \
                all(isinstance(item, Motif) for item in value):
            value_dimension = 1
        # value is a 2D list of motifs
        elif isinstance(value, list) and \
                all(isinstance(item, (list, tuple)) for item in value) and \
                all(isinstance(m, Motif) for sublist in value for m in sublist):
            value_dimension = 2
        else:
            raise ValueError(f'Only motifs, lists of motifs, or 2D lists of motifs can be added to the Origami, but the object {value} was added.')
        
        ### CHECK THE DIMENSIONALITY OF THE SLICING
        mask = None
        if hasattr(key, '__call__'): # if we submit a function, return the matrix filtered by the function
            # Select the motifs that satisfy the function
            mask = [[(i, slice(j, j+1)) 
                        for j, m in enumerate(row) if key(m)]
                            for i, row in enumerate(self._matrix)]
            
        elif isinstance(key, (tuple, list)) and\
                all(isinstance(i, int) for i in key):
            key_dimension = 0
            # Convert the keys to a positive integer
            y_int = key[0] % len(self._matrix)
            x_int = key[1] % len(self._matrix[y_int])

            # Convert the keys to slices
            y_slice = slice(y_int, y_int + 1)
            x_slice = slice(x_int, x_int + 1)

        elif isinstance(key, int):
            key_dimension = 1
            y_int = key % len(self._matrix)

            y_slice = slice(y_int, y_int + 1)
            x_slice = slice(0, len(self._matrix[key]))

        elif isinstance(key, (tuple, list)) and\
                isinstance(key[0], int) and isinstance(key[1], slice):
            key_dimension = 1 # select a row
            y_int = key[0] % len(self._matrix)

            y_slice = slice(y_int, y_int + 1)
            x_slice = key[1] # get the slice

        elif isinstance(key, (tuple, list)) and\
                all(isinstance(i, slice) for i in key):
            key_dimension = 2 # select a 2D region
            y_slice = key[0]
            x_slice = key[1]

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

        ### APPLY THE MASK
        if mask is not None:
            for row_ind, row in enumerate(mask):
                for mot_ind, (i, j) in enumerate(row):
                    if value_dimension in (0, 1):
                        self._matrix[i][j] = value
                    elif value_dimension == 2:
                        try:
                            self._matrix[i][j] = [value[row_ind][mot_ind]]
                        except IndexError as e:
                            raise ValueError(f'Error while setting the value to the Origami.\
                                            The lists do not match. Origami indexes: y: {i}, x: {j}.') from e

        ### REDUCE ALL CASES TO A 2D SLICING
        else:
            for i, line in enumerate(self._matrix[y_slice]):
                if value_dimension in (0, 1):
                    line[x_slice] = value
                elif key_dimension ==  2 and value_dimension == 2:
                    try:
                        line[x_slice] = value[i]
                    except IndexError as e:
                        raise ValueError(f'Error while setting the value to the Origami.\
                                         The lists do not match. Origami indexes: y: {i}, x: {x_slice}.') from e

        ### update the motif
        self._uniform_lines_len() # this also prepare the motifs
        self._updated_motif()

    def __len__(self):
        """ Get the number of lines in the origami"""
        return len(self._matrix)
    
    def __str__(self):
        return str(self.motif)


    def __repr__(self):
        #string in which the representation of the helix is saved
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
    def num_motifs(self):
        return sum(1 for line in self._matrix for item in line if item is not None)
    
    @property
    def num_char(self):
        if not self._matrix:
            return 0
        return [sum(m.num_char for m in line) for line in self._matrix]

    @property
    def num_lines(self):
        if not self._matrix:
            return 0
        return len(self._matrix)
    
    @property
    def num_columns(self):
        if not self._matrix:
            return 0
        return len(self._matrix[0])
    
    @property
    def ss_assembly(self):
        return self._ss_assembly
    
    @ss_assembly.setter
    def ss_assembly(self, new_ss_assembly):
        self._ss_assembly = bool(new_ss_assembly)
        self._updated_motif()
    
    @property
    def helices(self):
        if self._helices:
            return self._helices
        helix_matrix = []
        # scan the matrix line by line
        for line in self._matrix:
            new_line = []
            helix = []
            # scan the motifs in the line, avoiding empty motifs
            ind1 = 0
            ind2 = 1
            len_line = len(line)
            while ind2 < len_line:
                m1 = line[ind1]
                m2 = line[ind2]
                if not m2:
                    ind2 += 1
                    continue
                if not m1:
                    ind1 += 1
                    continue
                # check if the first and the second motif are connected, if yes, add them to the helix
                if m1.junctions[(1, 0)] and m2.junctions[(-1, 0)]:
                    if not helix:
                        helix.append(m1) # if the motif is the first in the helix, add it
                    helix.append(m2) # add the second motif to the helix
                elif helix: # if the motif are not connected and there is already an helix, save it to the line
                    new_line.append(helix)
                    helix = [] # reset the helix
                ind1 += 1
                ind2 += 1
            if helix: # if there is an helix at the end of the line, save it
                new_line.append(helix)
            helix_matrix.append(new_line)
        self._helices = helix_matrix
        return helix_matrix
    
    @property
    def junctions(self):
        """ Return the vertical junctions for each line of the Origami.
        The motif junctions are extended, there for we can just take the maximum position of the motifs."""
        ### check if the junctions are already calculated
        if self._junctions:
            return self._junctions
        junctions_list = []
        ### go through the lines and calculate the junctions
        for line in self._matrix:
            x_shift = 0
            # line = [m for m in line if m] # remove empty motifs           
            line_junc_dict = {(1,0): [], (0,1): [], (-1,0): [], (0,-1): []} # initialize the dictionary
            #                  right,    bottom,     left,      top
            for ind, m in enumerate(line):
                if ind == 0: # add the left junctions of the first motif
                    line_junc_dict[(-1,0)] = line[0].junctions[(-1,0)]
                ### add the top and bottom junctions
                line_junc_dict[(0, 1)] += [(x + x_shift, y) for x, y in m.junctions[(0, 1)]]
                line_junc_dict[(0, -1)] += [(x + x_shift, y) for x, y in m.junctions[(0, -1)]]
                x_shift += m.num_char # just add numchar becasue the junctions are extended
                # if it's the last motif, add the right junctions
                if ind == len(line) - 1:
                    line_junc_dict[(1, 0)] += [(x + x_shift, y) for x, y in m.junctions[(1, 0)]]
            junctions_list.append(line_junc_dict)
        # save the junctions
        self._junctions = junctions_list
        return junctions_list

    @property
    def align(self):
        return self._align
    
    @align.setter
    def align(self, new_align):
        if new_align not in ('left', 'first', 'center'):
            raise ValueError(f'"{new_align}" is not an accepted value for the align_type variable. The align_type variable may only a string reading "left", "first" or "center".')
        self._align = new_align
        self._updated_motif()

    @property
    def assembled(self):
        if not self._assembled:
            self._assemble()
        return self._assembled
    
    @property
    def map(self):
        """ Return a dictionary with the position of the symbols in the origami as keys and the (y, x) position in the matrix as values."""
        if self._map is not None:
            return self._map
        
        ### build the map and the motif map
        assembled = self.assembled
        self._map = {}
        self._motif_map = {}
        for y, line in enumerate(assembled):
            if y % 2 == 1: # skip the connection motifs
                continue
            for x, m in enumerate(line):
                for pos in m.map:
                    # y//2 because for each line there is a connection motif
                    self._map[pos] = (y//2, x)
                    self._motif_map[pos] = m
        return self._map
    
    @property
    def shift_map(self):
        """ Return a dictionary with the position of the (y, x) coordinates in the matrix as keys and the (x, y) shift of the motifs as values."""
        if self._shift_map:
            return self._shift_map
        self._assemble()
        return self._shift_map
    
    @property
    def motif_map(self):
        """ Motif map is a dictionary with positions as key, and motifs as values"""
        if not self._map:
            self.map
        return self._motif_map
    
    @property
    def sequence_index_map(self):
        """ Sequence map is a dictionary with positions as key, and sequence indexes as values"""
        if not self._motif:
            self.motif
        return self.motif.sequence_index_map
    
    @property
    def base_map(self):
        """ Base map is a dictionary with positions as key, and motif slices as values"""
        if not self._base_map:
            self._base_map = {}
            origami_map = self.map
            origami_motif_base_map = self.motif.base_map
            for pos in origami_motif_base_map.keys():
                self._base_map[pos] = origami_map[pos]
        return self._base_map
    
    @property
    def sequence(self):
        return self.motif.sequence
    
    @sequence.setter
    def sequence(self, new_seq):
        """Set the sequence of the origami"""
        if not isinstance(new_seq, (str, Sequence)) or len(new_seq) != len(self.sequence):
            raise ValueError(f"The new sequence must be a string or a Sequence object with the same lenght of the current sequence ({len(self.sequence)}). Got type: {type(new_seq)}; with length: {len(new_seq)}")
        new_seq = new_seq.replace('&', '') # remove the '&' symbol
        offset = 0 # adjust the offset if there are multiple strands
        origami_motif = self.motif
        motif_map = self.map
        motif_shifts = self.shift_map
        for s in origami_motif:
            mot_str_coor = None
            current_base_map = dict()
            for ind, pos in enumerate(s.base_map):
                # get the y,x cooridnates of the motif in the matrix
                motif_yx = motif_map[pos]
                # get the x,y shift of the motif in the origami positions
                shift_yx = motif_shifts[motif_yx]
                # remove the shifts from the position of the base
                shifted_pos = (pos[0] - shift_yx[0], pos[1] - shift_yx[1])
                # get the motif at the position
                motif = self._matrix[motif_yx[0]][motif_yx[1]]
                # get the strand index of the motif at the base position
                strand_ind = motif.map[shifted_pos]
                # if the strand is a new strand, add the current base map to the strand
                if mot_str_coor != (motif_yx[0], motif_yx[1], strand_ind):
                    # if there is a current base map, add it to the strand
                    if current_base_map:
                        # get the strand, at the strand index and add the curent base map
                        self._matrix[mot_str_coor[0]][mot_str_coor[1]][mot_str_coor[2]].sequence = ''.join(current_base_map.values())
                    # update the current base map with the new strand
                    current_base_map = motif[strand_ind].base_map
                    # update the motif and strand position to the current position
                    mot_str_coor = (motif_yx[0], motif_yx[1], strand_ind)
                # update the base in the base map
                current_base_map[shifted_pos] = new_seq[ind + offset]
            # add the last strand
            self._matrix[mot_str_coor[0]][mot_str_coor[1]][mot_str_coor[2]].sequence = ''.join(current_base_map.values())
            # update the offset 
            offset += len(s.sequence)

    @property
    def pseudoknots(self):
        if self._pseudoknots:
            return self._pseudoknots
        
        def add_pk(strand, pk_index, info_nr, shift, forward=True):
            pk_info = strand.pk_info
            # add the energy and tolerance of the motif
            pk_dict.setdefault(pk_index, {"id": pk_index, 'ind_fwd': [], 'ind_rev': [], 'E': [], 'dE': []})
            pk_dict[pk_index]['E'].append(pk_info['E'][info_nr])
            pk_dict[pk_index]['dE'].append(pk_info['dE'][info_nr])

            # indicate the index of the pseudoknot in the sequence
            start_pos_ind = 0 if strand.directionality == '53' else -1
            pos = list(strand.base_map.keys())[start_pos_ind]
            # get the index of the sequence in the strand
            offset_ind = self.sequence_index_map[(shift[0] + pos[0], shift[1] + pos[1])]
            pk_start, pk_end = pk_info['ind_fwd'][info_nr]
            start_end_tuple = (offset_ind + pk_start, offset_ind + pk_end)
            if forward:
                pk_dict[pk_index]['ind_fwd'].append(start_end_tuple)
            else:
                pk_dict[pk_index]['ind_rev'].append(start_end_tuple)

        pk_motifs = []
        # collect all the motifs with pk_index
        for i, line in enumerate(self._matrix):
            for j, m in enumerate(line):
                if any(hasattr(s, 'pk_info') for s in m):
                    pk_motifs.append((i, j))
        
        pk_dict = dict()
        for i, j in pk_motifs:
            m = self._matrix[i][j]
            shift = self.shift_map[(i, j)]
            pk_strands = [s for s in m if hasattr(s, 'pk_info')]
            pk_indexes = [pk_id for s in pk_strands for pk_id in s.pk_info['id']]

            ### Adjust the pk_index for the unique pseudoknots
            if any(ind[0] == '0' for ind in pk_indexes): # new 0 pseudoknot
                current_n_zero = sum(1 if key[0] == '0' else 0 for key in pk_dict)
                pk_index_0 = '0_' + str(current_n_zero + 1)            

            # add the pseudoknots
            for strand in pk_strands:
                for info_nr, pk_index in enumerate(strand.pk_info['id']):
                    reverse = pk_index[-1] == "'"
                    if pk_index[0] == '0':
                        pk_index = pk_index_0
                    add_pk(strand, pk_index, info_nr, shift, forward=not reverse)

        # make the average energy and average tolerance
        for pk in pk_dict.values():
            pk['E'] = sum(pk['E']) / len(pk['E'])
            pk['dE'] = sum(pk['dE']) / len(pk['dE'])
        # convert the pk_dict to a list for simplicity
        self._pseudoknots = list(pk_dict.values())
        return self._pseudoknots

    @property
    def structure(self):
        return self.motif.structure
    
    @property
    def pair_map(self):
        if not self._motif:
            self.motif
        return self.motif.pair_map
    
    @property
    def motif(self):
        if isinstance(self._motif, Motif):
            return self._motif
        self._motif = Motif.concat([Motif.concat(line, align=False, unlock_strands=self._ss_assembly, lock_coords=False) for line in self.assembled], axis=0, align=False, lock_coords=True, unlock_strands=self.ss_assembly)
        return self._motif
    
    @property
    def strands(self):
        return list(self.motif)

    ### 
    ###  STATIC METHODS
    ###

    @staticmethod
    def _calculate_connections(junctions1, junctions2, junctions_x_shift=(0, 0), start_y=0):
        """
        Creates the connection between the helices

        Parameters
        ----------
        junctions1: dict
            junctions of the first line
        junctions2: dict
            junctions of the second line
        junctions_x_shift: tuple
            The x shift of the bottom and top jucntions, respectively
        start_y: int
            The y position of the first line
        """
        ### take the junctions of the two lines and check them
        j1 = [pos[0] + junctions_x_shift[0] for pos in junctions1[(0, 1)]] # take bottom junctions x positions, add the shift
        j2 = [pos[0] + junctions_x_shift[1] for pos in junctions2[(0, -1)]] # take top junctions x positions, add the shift
        if not j2 or not j1:
            return Motif(), 0
        # the number of connections is the minimum of the two junctions
        n_connect =  min((len(j1), len(j2)))
        j1 = j1[:n_connect]
        j2 = j2[:n_connect]
        
        ### create the connections between the junctions j1 and j2
        closed_crossings = dict() # a dictionary with the connented pair index as key and a set of crossed pair indexes as value
        pairs = list(zip(j1, j2)) # zip the two lists
        for ind, (x1, x2) in enumerate(pairs): 
            closed_crossings[ind] = set() # intialize the crossed pair indexes
            x_min = min(x1, x2) # the minimum x position to connect
            x_max = max(x1, x2) # the maximum x position to connect
            # the crossed pairs are pairs that have at list one x position between the minimum and maximum x positions and are not already connected
            crossed = {i for i, x12 in enumerate(pairs) if i not in closed_crossings and (x_min <= x12[0] <= x_max or x_min <= x12[1] <= x_max)}
            # update the crossed pair indexes
            closed_crossings[ind].update(crossed)
        # check for nested crossing that are not directly overlapping, like j1 = [22, 23, 33, 34], j2 = [2, 3, 22, 23]
        ### update nested crossings
        for key1 in list(closed_crossings.keys()): # go through the connected pairs
            for key2 in list(closed_crossings.keys()): 
                if key2 in closed_crossings[key1]: # if the second pair is in the crossed pairs of the first pair
                    closed_crossings[key1].update(closed_crossings[key2])

        # calculate the maximum number of crossings
        max_crossing = max(len(crossed) for crossed in closed_crossings.values())

        ### make the connections
        strands = [] # a list of the strands
        for ind, (x1, x2) in enumerate(pairs): 
            n_crossings = len(closed_crossings[ind])
            if x1 < x2: # the first motif is on the left
                strand = "│" * n_crossings + "╰" + "─" * (x2 - x1 -1) + "╮" + "│" * (max_crossing - n_crossings)
            elif x1 > x2: # the first motif is on the right
                strand = "│" * (max_crossing - n_crossings) + "╯" + "─" * (x1 - x2 -1) + "╭" + "│" * n_crossings
            else: # the motifs are on the same position
                strand = "│" * (max_crossing + 1)
            # can add the symbol "^" for retrocompatibility with ROAD
            strand += "↑" # if you  do this, increase the max_crossing by 1
            strands.append(Strand(strand, start=(x1, start_y), direction=(0, 1)))
        return Motif(strands, join=False), max_crossing + 1 + 1 # Extra +1 to the max_crossing to add the symbol "^" for retrocompatibility with ROAD
    
    ### 
    ### PROTECTED METHODS
    ###

    def _prepare_motif(self, motif, y=None, x=None, copy=False):
        """
        Check if an element is of class Helix or Motif (or subclass of Motif)
        """
        if motif is None:
            return Motif()
        if not isinstance(motif, Motif):  
            error_string = f'Only Motif or subclass of Motif can be added to the Origami, but the object {motif} was added'
            if y is not None and x is not None:
                error_string += f' at position x: {x}, y: {y}'
            raise ValueError(error_string)
        if copy:
            motif = motif.copy()
        motif.register_callback(self._updated_motif)
        # motif.extend_junctions()
        return motif
    
    def _uniform_lines_len(self, copy=False):
        ### clean the empty motifs
        for ind, line in enumerate(self._matrix):
            self._matrix[ind] = [self._prepare_motif(motif, copy=copy) for motif in line if motif]

        ### align the motifs in the lines and extend the junctions
        for ind, line in enumerate(self._matrix):
            aligned = Motif.align(line) #+ [Motif()] * (max_len - len(line)))
            max_y = max([m.max_pos[1] for m in aligned], default=0)
            for m in aligned:
                # m.extend_junctions(skip_axis = 1, until= (None, max_y))
                m.extend_junctions(until= (None, max_y))
            self._matrix[ind] = aligned

    def _assemble(self):
        """
        Return a matrix where the motifs and the jucntions of the Origami are in place
        """
        self._uniform_lines_len()

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
               connection_motif, connection_height = self._calculate_connections(junctions[ind], junctions[ind + 1], junctions_x_shift = junctions_shift, start_y = vert_shift)
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


    def save_3d_model(self, filename: str = 'origami', forces: bool = False, return_text: bool = False, pdb=False, **kwargs):
        """
        Save the origami as an oxDNA file
        """
        return self.motif.save_3d_model(filename, forces=forces, return_text=return_text, pdb=pdb, **kwargs)

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
        self._uniform_lines_len(copy=copy) # remove empty motifs and add callbacks
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
        self._uniform_lines_len(copy=copy)
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
        self._uniform_lines_len()
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
        self._uniform_lines_len()
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
        self._uniform_lines_len()
        self._updated_motif()

    def get_motif_type(self, motif_type):
        """ Return a list of motifs of the given type """
        return [m for line in self._matrix for m in line if type(m) == motif_type]

    ### IMPLEMENT THIS!!!
    def find_best_start(self):
        return
        ### make a strand map of sequence index --> motif position
    
    def _updated_motif(self, *args, **kwargs):
        self._helices = None
        self._motif = None
        self._junctions = None
        self._map = None
        self._base_map = None
        self._shift_map = None
        self._assembled = None
        self._pseudoknots = None
        self._trigger_callbacks(*args, **kwargs)

    def to_road(self):
        ori_str = str(self)
        ori_str = ori_str.replace('↑', '^')
        ori_str = ori_str.replace('↓', '^')
        ori_str = ori_str.replace('│ ┊┊┊┊┊┊ │', '│ ****** │')
        ori_str = ori_str.replace(' ┊┊ ', ' !! ')
        ori_str = ori_str.replace(' ┊ ', ' ! ')
        return ori_str
    
    def save_text(self, filename_path, to_road=False):
        path = filename_path.split('.')[0]
        name = path.split('/')[-1].split('\\')[-1]
        with open(path + '.txt', 'w') as f:
            f.write(f'>{name}\n')
            f.write(self.sequence + '\n')
            f.write(self.structure + '\n\n')
            if to_road:
                f.write(self.to_road())
            else:
                f.write(str(self))

    def save_structure_text(self, filename_path):
        self.motif.save_structure_text(filename_path)

    def copy(self):
        return Origami(self._matrix, ss_assembly=self.ss_assembly, align=self.align, copy=True)
