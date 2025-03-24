import random
import warnings
from typing import Tuple
from .basepair import BasePair

###
### USEFUL SETS AND DICTIONARIES
###

nucleotides = {"A", 
               "U", 
               "C", 
               "G", # standard nucleotides
               "W", # A or U
               "M", # A or C
               "R", # A or G
               "Y", # C or U
               "K", # G or U
               "S", # G or C  
               "D", # A or G or U
               "H", # A or C or U
               "V", # A or C or G
               "B", # C or G or U
               "N", # any nucleotide
               "X", # any nucleotide for external Kissing Loops in ROAD
               '&', # separator
               }

# https://www.bioinformatics.org/sms/iupac.html
iupac_code = {"A": {"A"}, 
              "U": {"U"}, 
              "C": {"C"}, 
              "G": {"G"}, # standard nucleotides
              "W": {"A", "U"}, # A or U
              "M": {"A", "C"}, # A or C
              "R": {"A", "G"}, # A or G
              "Y": {"C", "U"}, # C or U
              "K": {"G", "U"}, # G or U
              "S": {"G", "C"}, # G or C
              "D": {"A", "G", "U"}, # A or G or U
              "H": {"A", "C", "U"}, # A or C or U
              "V": {"A", "C", "G"}, # A or C or G
              "B": {"C", "G", "U"}, # C or G or U
              "N": {"A", "U", "C", "G"}, # any nucleotide
              "X": {"A", "U", "C", "G"}, # any nucleotide for external Kissing Loops in ROAD
              "&": {"&"}, # separator
              }

base_pairing = {"A": {"B", "D", "H", "K", "N", "U", "W", "Y"},
                "U": {"A", "B", "D", "G", "H", "K", "M", "N", "R", "S", "V", "W"},
                "C": {"B", "D", "G", "K", "N", "R", "S", "V"},
                "G": {"B", "C", "D", "H", "K", "M", "N", "S", "U", "V", "W", "Y"}, # standard base pairing
                "W": {"A", "B", "D", "G", "H", "K", "M", "N", "R", "S", "U", "V", "W", "Y"},
                "M": {"B", "D", "G", "H", "K", "N", "R", "S", "U", "V", "W", "Y"},
                "R": {"B", "C", "D", "H", "K", "M", "N", "S", "U", "V", "W", "Y"},
                "Y": {"A", "B", "D", "G", "H", "K", "M", "N", "R", "S", "V", "W"},
                "K": {"A", "B", "C", "D", "G", "H", "K", "M", "N", "R", "S", "U", "V", "W", "Y"},
                "S": {"B", "C", "D", "G", "H", "K", "M", "N", "R", "S", "U", "V", "W", "Y"},
                "D": {"A", "B", "C", "D", "G", "H", "K", "M", "N", "R", "S", "U", "V", "W", "Y"},
                "H": {"A", "B", "D", "G", "H", "K", "M", "N", "R", "S", "U", "V", "W", "Y"},
                "V": {"B", "C", "D", "G", "H", "K", "M", "N", "R", "S", "U", "V", "W", "Y"},
                "B": {"A", "B", "C", "D", "G", "H", "K", "M", "N", "R", "S", "U", "V", "W", "Y"},
                "N": {"A", "B", "C", "D", "G", "H", "K", "M", "N", "R", "S", "U", "V", "W", "Y"},
                "X": {"A", "B", "C", "D", "G", "H", "K", "M", "N", "R", "S", "U", "V", "W", "Y", "X"},
                "&": {"&"} # separator
                }

db_pairs = {"(": ")", "[": "]", "{": "}", "<": ">", "A": "a", "B": "b", "C": "c", "D": "d", "E": "e", "F": "f", "G": "g", "H": "h", "I": "i", "J": "j", "K": "k", "L": "l", "M": "m", "N": "n", "O": "o", "P": "p", "Q": "q", "R": "r", "S": "s", "T": "t", "U": "u", "V": "v", "W": "w", "X": "x", "Y": "y", "Z": "z"}
all_pk_symbols =     ('[', ']', '{', '}', '<', '>', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z')
road_symbols =  nucleotides | {"." "(", ")"} | set(all_pk_symbols) | {
                "─", "│", "╭", "╮", "╰", "╯", "^", "*", "┼", "┊", '~', '◦',  # ROAD special symbols
                '↑', '↓', "⊗", "⊙", "▄", "█", # Pyfurnace special symbols
                "-", "|", "+", ":", "=", "!", " ", "/", "\\", "3", "5", "&"} 
bp_symbols = {"┊", "=", ":", "!", "*"}
# turns = {"╭", "╮", "╰", "╯", "/", "\\"} 

###
### USEFUL STRING TRANSLATIONS
###

def complement(sequence):
    return sequence.translate(nucl_to_pair)

def reverse_complement(sequence):
    return sequence.translate(nucl_to_pair)[::-1]

pseudo_to_dot = str.maketrans(''.join(all_pk_symbols), "." * len(all_pk_symbols))
pair_db_sym = str.maketrans(''.join(db_pairs.keys()) + ''.join(db_pairs.values()), 
                            ''.join(db_pairs.values()) + ''.join(db_pairs.keys()))
nucl_to_none = str.maketrans("", "", "".join(nucleotides))
symb_to_none = str.maketrans("", "", "".join(road_symbols))
nucl_to_pair = str.maketrans("AUCGWMRYKSDHVBNX", "UAGCWKYRKSHDBVNX")
bracket_to_none = str.maketrans("", "", "{[(.)]}")
# only_turns = str.maketrans("", "", "".join(road_symbols-turns))
only_nucl = str.maketrans("", "", "".join(road_symbols-nucleotides))
horiz_flip = str.maketrans("╭╮╰╯/\\", "╮╭╯╰\\/")
verti_flip = str.maketrans("╭╮╰╯/\\", "╰╯╭╮\\/")
rotate_90 = str.maketrans("╭╮╰╯│|─-", "╮╯╭╰──││")
symb_to_road = str.maketrans("-|+=:*!", "─│┼┊┊┊┊")

def pair_nucleotide(nucleotide, iupac_symbol='N'):
    """ Return the nucleotides that can pair with the given nucleotide"""
    if iupac_symbol in 'AUCG':
        return iupac_symbol
    if iupac_symbol == 'K':
        if nucleotide == 'G':
            return 'U'
        elif nucleotide == 'U':
            return 'G'
    return nucleotide.translate(nucl_to_pair)

def mutate_nucleotide(sequence, sequence_constraints, nucl_ind, pair_map):
    """ Take a nucleotide, and return a possible mutation and the corresponding pair"""
    new_set = iupac_code[sequence_constraints[nucl_ind]] - {sequence[nucl_ind]}
    if not new_set:
        return None, None
    new_nucl = random.choice(list(new_set))
    paired_nucl = pair_map[nucl_ind] # get the index of the paired nucleotide
    if paired_nucl is not None: # if the nucleotide is paired, get the paired nucleotide
        # get the possible nucleotides that can pair with the new nucleotide
        paired_nucl = pair_nucleotide(new_nucl, sequence_constraints[paired_nucl])
    return new_nucl, paired_nucl

def gc_content(seq, extended_alphabet=True):
    """Calculate the percentage of G, C, S and half K in the sequence.
    --------------------------------------------------------------------------------------
    Return: float
        The percentage of G, C, or S in the sequence.
    """
    total_count = len(seq)
    if not total_count:
        return 0
    gc_count = seq.count('G') + seq.count('C') 
    if extended_alphabet:
        gc_count += seq.count('S') + sum(map(seq.count, ['M', 'R', 'Y', 'K', 'N', 'X'])) / 2 + sum(map(seq.count, ['D', 'H'])) / 3 + sum(map(seq.count, ['V', 'B'])) * 2 / 3 
    percentage = (gc_count / total_count)
    return percentage

###
### DOT BRACKET FUNCTIONS
###

def rotate_dot_bracket(dot_bracket, shift_left):
    n = len(dot_bracket)
    
    # Step 1: Get the pair map
    pair_map = dot_bracket_to_pair_map(dot_bracket)
    
    # Step 3: Adjust the pair map
    new_pair_map = BasePair()
    for i, j in pair_map.items():
        new_i = (i - shift_left) % n
        new_j = (j - shift_left) % n if j is not None else None
            
        new_pair_map[new_i] = new_j
    
    return pair_map_to_dot_bracket(new_pair_map, n)
    
def dot_bracket_to_pair_map(dot_bracket):
    """ Convert the dot bracket notation to a double dictionary with the paired indexes"""
    pair_map = BasePair()
    stacks = [[] for _ in range(len(db_pairs))]
    for i, sym in enumerate(dot_bracket):
        if sym in db_pairs.keys():
            stack_num = list(db_pairs.keys()).index(sym)
            stacks[stack_num].append(i)
        elif sym in db_pairs.values():
            stack_num = list(db_pairs.values()).index(sym)
            if stacks[stack_num]: # pop the last element if the stack is not empty
                j = stacks[stack_num].pop()
                pair_map[j] = i
        else: # unpaired nucleotide: No pair
            pair_map[i] = None
    # If there are still elements that are not in the dictionary, add them as unpaired
    for i in range(len(dot_bracket)):
        if i not in pair_map:
            pair_map[i] = None
    return pair_map

def pair_map_to_dot_bracket(pair_map, structure_length = None):
    """ Convert the pair map to a dot bracket notation"""
    if structure_length is None:
        structure_length = max(j if j is not None else 0 for i in pair_map.items() for j in i) + 1

    ### CREATE THE DOT BRACKET NOTATION ###
    done_pairs = set()

    # Prepare the variables
    dotbracket = ['.'] * structure_length 
    bracket_count = 0
    recheck = True
    stack = []

    while recheck: # recheck the structure, for pseudoknots
        recheck = False

        for i in range(structure_length):

            paired = pair_map[i]
            if paired is None or i in done_pairs or paired in done_pairs:
                continue # ignore unpaired positions or already paired positions

            if paired > i: # the current position will close the pair later
                
                # may be a pseudknot:  Check for clash
                if stack and paired > stack[-1]:
                    recheck = True # clash detected: recheck the structure later
                    continue
                else: # no clash: add the pair to the stack
                    stack.append(paired)
                    dotbracket[i] = list(db_pairs.keys())[bracket_count]
                    dotbracket[paired] = list(db_pairs.values())[bracket_count]
                    
            else: # paired already analyzed
                if not stack: # no pair in the stack
                    continue
                if i == stack[-1]: # the current position closes the last pair
                    # Remove the pair from the map
                    done_pairs.add(i)
                    done_pairs.add(paired)
                    stack.pop()

        stack = [] # reset the stack count
        bracket_count += 1

        if bracket_count >= 30:
            warnings.warn(f"Warning: Too many bracket types needed in to write the structure. Stopping at 30, 'Z' to 'z'.", stacklevel=3)
            break

    return ''.join(dotbracket)

def dot_bracket_to_stacks(dot_bracket, only_opening=False):
    pair_map = dot_bracket_to_pair_map(dot_bracket)
    stacks = [] # a list of tuples containing starting index and the length of the stack
    reduced_dot_bracket = [] # a list of the reduced dot bracket symbols, that will be returned as a string

    stack_start = 0
    current_pair_sym = dot_bracket[0]
    last_stack_pair = pair_map[0] + 1 if pair_map[0] else None
    ### NOT NECESSARY # open_pair_sym = '.' + ''.join(db_pairs.keys())
    for i, symbol in enumerate(dot_bracket + '_'):
        ### NOT NECESSARY 
            # # ignore the symbols that are not '.' or not opening pairs, they should be already taken into accound
            # if symbol in open_pair_sym:
            #     continue

        stack_pair = pair_map.get(i) # get the paired index of the current symbol
        # if symbol change, or the pairing is not consecutive, add the last stack to the list
        if symbol != current_pair_sym or stack_pair is not None and stack_pair != last_stack_pair - 1:
            add_last_stack = True # add the last stack to the list by default
            if only_opening: # if only the opening stacks are needed, check if the last stack is an opening stack or unpaired
                add_last_stack = current_pair_sym in db_pairs.keys() or current_pair_sym == '.'

            if add_last_stack:
                reduced_dot_bracket.append(current_pair_sym)
                stacks.append((stack_start, i - 1))
            stack_start = i

        last_stack_pair = stack_pair
        current_pair_sym = symbol
    return ''.join(reduced_dot_bracket), stacks

def folding_barriers(structure: str,
                     kl_delay: int = 150) -> Tuple[str, int]:
    """
    Compute the folding barriers for a given RNA secondary structure.
    This function is based on ROAD: https://www.nature.com/articles/s41557-021-00679-1
    This function analyzes the dot-bracket representation of the RNA secondary structure 
    to determine folding barriers based on base pair topology and kinetic delay.

    Parameters
    ----------
    structure : str
        The dot-bracket notation representing the secondary structure of the RNA. 
        If folding barriers is called from a Motif or Origami, the structure is already
        provided by the object in the dot-bracket format.
    kl_delay : int, optional, default=150
        The number of nucleotides (nts) of delay before kissing loops (KLs) snap closed.
        A typical realistic value is around 350 (~1 second), while a more conservative 
        setting is 150 by default.

    Returns
    -------
    Tuple[str, int]
        A tuple containing:
        - A string representing the barrier map where:
        '─' = no barrier (0 penalty),
        '▂' = opening pair weak barrier (1 penalty),
        '▄' = cloding pair weak barrier (1 penalty),
        '█' = strong barrier (2 penalty).
        - An integer score indicating the total penalty based on barrier strengths.

    Notes
    -----
    The function assigns barrier strengths based on the topology of base pairs 
    and kinetic constraints, ensuring proper folding predictions.
    """
    structure = structure.replace('&', '')
    ss_len = len(structure)
    s_map = dot_bracket_to_pair_map(structure)
    barriers = [''] * len(structure)

    # don't count the nucleotides in the 3' terminal position
    terminal_3_bonus = 16

    current_barr_count = 0
    for ind in range(ss_len):
        # get the bracket at the position
        bra = structure[ind]
        # if bra == '&':
        #     barriers[ind] = '&'
        #     continue

        # set the default barrier to 'N'
        barriers[ind] = '─'

        if ind > ss_len - terminal_3_bonus:
            if s_map[ind] is not None:
                barriers[s_map[ind]] = '─'
            continue

        ### unpaired nts are not barriers and reset the current barrier count
        if bra == '.':
            barriers[ind] = '─'
            current_barr_count = 0

        ### opening pairs may be barriers or not, 
        # indicate them with 'w', reset the barrier count
        elif bra == '(':
            barriers[ind] = '▂'
            current_barr_count = 0

        ### closing pairs may be barriers or not, 
        # depending on the topology count
        elif bra == ')':
            # if the opening pair is not blocked, 
            # we close and reset the current barrier count
            if barriers[s_map[ind]] == '▂':
                barriers[ind] = '─'
                barriers[s_map[ind]] = '─'
                current_barr_count = 0

            # if the closing pair is already blocked, then 
            # we mark the barrier reached and start counting
            elif barriers[s_map[ind]] == '█':
                # if the current barrier count is greater than 5, wa mark it with 'S'
                if current_barr_count > 5:
                    barriers[ind] = '█'
                    barriers[s_map[ind]] = '▄'
                # if the current barrier count is less than 5, we mark it with 'W'
                else:
                    barriers[ind] = '▄'
                    barriers[s_map[ind]] = '▄'
                current_barr_count += 1

        # if the index is bigger than the kl_delay, 
        # we check if the internal KLs are closed
        if ind > kl_delay:
            
            # we check if the KLs are closed, 
            # if yes, we mark them with 'N'
            close_sym = structure[ind - kl_delay]
            if close_sym not in '(.)' and close_sym in db_pairs.values():
                barriers[ind - kl_delay] == '─'
                barriers[s_map[ind - kl_delay]] == '─'

                # for each nt in the delay, we check if there are any barriers, 
                # if yes, we mark them with 'S'
                for k in range(s_map[ind - kl_delay], ind-kl_delay):
                    if barriers[k] == '▂':
                        barriers[k] = '█'

    penalty = {'█': 2, '▄': 1, '▂': 1, '─': 0}#, '&': 0}
    repl = [penalty[i] for i in barriers]
    score = sum(repl)
    return ''.join(barriers), score

class Node():
    def __init__(self, index = None, label = "5'", parent=None, seq=None, **kwargs):
        self.index = index
        self.label = label
        self.paired_index = None
        self.parent = parent
        self.seq = seq
        self.children = []
        if parent is not None:
            parent.children.append(self)
        for key, value in kwargs.items():
            setattr(self, key, value)

    def __repr__(self, level=0):
        paired_str = f' -> {self.paired_index}' if self.paired_index is not None else ''
        if level == 0:
            ret = f"{self.label}\n"
        else:
            ret = f"╰{'──' * level} {self.label} - {self.seq} - {self.index} {paired_str} \n" 
        for child in self.children:
            ret += child.__repr__(level + 1)
        return ret

    def __str__(self):
        return self.__repr__()

    def search(self, target_index, target_label=None):
        # Check if the current node matches the target criteria
        if target_index == self.index and (target_label is None or target_label == self.label):
            return self
        elif self.index is None and self.parent is None and abs(target_index) == float('inf'):
            return self # special case for the root node

        # Recursively search in child nodes
        for child in self.children:
            result = child.search(target_index=target_index, target_label=target_label)
            if result:
                return result  # Return the node if found
    
        return None  # Return None if not found

def dot_bracket_to_tree(dot_bracket, sequence=None):
    # Remove pseudoknots from the dot bracket
    dot_bracket = dot_bracket.translate(pseudo_to_dot)
    root = Node()
    current_node = root
    for i, label in enumerate(dot_bracket):
        if sequence is not None:
            seq_const = sequence[i]
        else:
            seq_const = None
        if label == "(": # add a new node, move to the new node
            new_node = Node(i, label, current_node, seq = seq_const)
            current_node = new_node
        elif label == ")": # add the paired index to the parent, move to the parent node
            current_node.paired_index = i
            current_node = current_node.parent
        elif label == ".": # add a children to the current node, do not move
            new_node = Node(i, label, current_node, seq = seq_const)
        elif label == '&':
             new_node = Node(None, label, current_node)
    return root

def tree_to_dot_bracket(node, dot_bracket=None, seq_constraints=False):

    if isinstance(seq_constraints, bool) and seq_constraints:
        # initialize the sequence constraints only if it's a bool
        seq_constraints = ["N"] # assume the sequence is at least one nucleotide long

    if dot_bracket is None:
        dot_bracket = ["."] # assume the structure is at least one nucleotide long

    if node.parent is not None: # don't add the root node
        if node.index >= len(dot_bracket):
            add_length = node.index - len(dot_bracket) + 1
            dot_bracket += ["."] * add_length
            if seq_constraints:
                seq_constraints += ["N"] * add_length

        dot_bracket[node.index] = node.label
        if seq_constraints and node.seq:
            seq_constraints[node.index] = node.seq

        if node.paired_index is not None:

            if node.paired_index >= len(dot_bracket):
                add_length = node.paired_index - len(dot_bracket) + 1
                dot_bracket += "." * add_length
                if seq_constraints:
                    seq_constraints += ["N"] * add_length

            dot_bracket[node.paired_index] = ")"
            if seq_constraints and node.seq:
                seq_constraints[node.paired_index] = node.seq.translate(nucl_to_pair)

    for child in node.children: # recursively add the children
        tree_to_dot_bracket(child, dot_bracket, seq_constraints)
    
    if node.parent is None: # we reached back the root node
        if seq_constraints:
            return ''.join(dot_bracket), ''.join(seq_constraints)
        
        return ''.join(dot_bracket)

def hamming_distance(s1, s2, ignore_ind=(), **kwargs):
    return sum((1 for i, (x, y) in enumerate(zip(s1, s2)) if x != y and i not in ignore_ind))

def base_pair_difference(s1, s2, pair_map1=None, pair_map2=None, ignore_ind=(), accept_unpaired_ind=(), **kwargs):
    """ Accpet unpaired indexes is a list of indexes that are allowed to be unpaired in the second structure."""
    if pair_map1 is None:
        pair_map1 = dot_bracket_to_pair_map(s1)
    if pair_map2 is None:
        pair_map2 = dot_bracket_to_pair_map(s2) 

    # Determine if an index should be considered based on the unpaired indices
    def check_ind(main_pair_map, i, accept_unpaired=False):
        if i in ignore_ind or main_pair_map[i] is None:
            return False
        elif accept_unpaired and i in accept_unpaired_ind:
            return pair_map2[i] is not None and pair_map1[i] != pair_map2[i]  # Only count if paired in both but different
        return pair_map1[i] != pair_map2[i]

    # pairs in s1 that are not in s2
    diff_set = [i for i in pair_map1 if check_ind(pair_map1, i, True)]  # Pairs in s1 not matching in s2
    # pairs in s2 that are not in s1
    diff_set.extend(i for i in pair_map2 if check_ind(pair_map2, i))  # Pairs in s2 not matching in s1
    
    return diff_set

def base_pair_distance(s1, s2, pair_map1=None, pair_map2=None, ignore_ind=(), accept_unpaired_ind=(), **kwargs):
    """ Accpet unpaired indexes is a list of indexes that are allowed to be unpaired in the second structure."""
    return len(base_pair_difference(s1, s2, pair_map1, pair_map2, ignore_ind, accept_unpaired_ind, **kwargs))

def differently_paired_distance(s1, s2, pair_map1=None, pair_map2=None, ignore_ind=(), accept_unpaired_ind=(), **kwargs):
    """ Accpet unpaired indexes is a list of indexes that are allowed to be unpaired in the second structure."""
    if pair_map1 is None:
        pair_map1 = dot_bracket_to_pair_map(s1)
    if pair_map2 is None:
        pair_map2 = dot_bracket_to_pair_map(s2) 
    # check each symbol, if they have the same pair or different pairs
    diff_pair_dist = 0
    for i in range(len(s1)):
        if i in ignore_ind:
            continue
        paired1 = pair_map1[i]
        paired2 = pair_map2[i]

        if paired1 is not None and paired2 is None:
            if i not in accept_unpaired_ind:
                diff_pair_dist += 1
        elif paired1 != paired2:
            diff_pair_dist += 1

    return diff_pair_dist

###
### SEQUENCE CONSTRAINT FUNCTIONS
###

def check_sequence_constraints(sequence_constraints, pair_map):
    """ Check if the sequence constraints agree with the pair map and update the sequence constraints"""
    new_seq_const = list(sequence_constraints)
    for i, j in pair_map.items():
        if j is None:
            continue
        iupac_i = iupac_code[sequence_constraints[i]]
        iupac_j = iupac_code[sequence_constraints[j]]
        # give the priority to the set with the least number of elements, the other index is the slave index
        priority_set = iupac_i
        slave_set, slave_ind = iupac_j, j
        if len(iupac_i) > len(iupac_j):
            priority_set = iupac_j
            slave_set, slave_ind = iupac_i, i
        paired_set = {nucl.translate(nucl_to_pair) for nucl in priority_set}
        # add wobble pairing
        if sequence_constraints[i] == 'K' or sequence_constraints[j] == 'K':
            paired_set.update('U' if nucl == 'G' else 'G' for nucl in priority_set)
        elif sequence_constraints[i] == 'U' and sequence_constraints[j] == 'G':
            paired_set.add('G'); paired_set.add('U')
        elif sequence_constraints[i] == 'G' and sequence_constraints[j] == 'U':
            paired_set.add('G'); paired_set.add('U')
        slave_intersection = paired_set.intersection(slave_set)
        if not slave_intersection:
            raise ValueError(f"Pairing constraint violated at index {i} and {j}:  nucleotide {i} ({sequence_constraints[i]}) is supposed to pair nucleotide {j} ({sequence_constraints[j]}) but the sequence constraints do not allow it.")
        new_iupac = next(key for key, value in iupac_code.items() if value == slave_intersection)
        new_seq_const[slave_ind] = new_iupac
    return ''.join(new_seq_const)

###
### CUSTOM EXCEPTIONS
###

class MotifStructureError(Exception):
    pass

class AmbiguosStructure(Warning):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return repr(self.message)

def visualize_tree(structure):
    from igraph import Graph, EdgeSeq
    import plotly.graph_objects as go
    def make_annotations(pos, text, font_size=10, font_color='rgb(250,250,250)'):
        L = len(pos)
        if len(text) != L:
            raise ValueError('The lists pos and text must have the same length')
        annotations = []
        for k in range(L):
            annotations.append(
                dict(
                    text=text[k],  # Use text[k] instead of labels[k]
                    x=pos[k][0], y=2 * M - position[k][1],
                    xref='x1', yref='y1',
                    font=dict(color=font_color, size=font_size),
                    showarrow=False
                )
            )
        return annotations

    # Define your RNA structure
    if type(structure) == str:
        root = dot_bracket_to_tree(structure)
    elif type(structure) == Node:
        root = structure
    else:
        raise ValueError("The structure should be a string or a Node object.")

    # Helper function to recursively add nodes and edges to the igraph graph
    def add_nodes_edges(igraph_graph, node, parent_index=None):
        current_index = igraph_graph.vcount()  # Get the current index for the new node
        igraph_graph.add_vertex(name=f"{node.label}-{node.index}")

        # If the node has a parent, create an edge to the parent
        if parent_index is not None:
            igraph_graph.add_edge(parent_index, current_index)

        # Recursively add children nodes
        for child in node.children:
            add_nodes_edges(igraph_graph, child, current_index)

    # Create an igraph Graph
    G = Graph(directed=True)
    add_nodes_edges(G, root)

    # Layout using igraph's reingold-tilford algorithm for tree layout
    lay = G.layout("rt")

    # Extract coordinates for Plotly
    nr_vertices = G.vcount()
    position = {k: lay[k] for k in range(nr_vertices)}
    Y = [lay[k][1] for k in range(nr_vertices)]
    M = max(Y)

    # Get edges and adjust positions
    es = EdgeSeq(G)
    E = [e.tuple for e in G.es]

    L = len(position)
    Xn = [position[k][0] for k in range(L)]
    Yn = [2 * M - position[k][1] for k in range(L)]
    Xe = []
    Ye = []
    for edge in E:
        Xe += [position[edge[0]][0], position[edge[1]][0], None]
        Ye += [2 * M - position[edge[0]][1], 2 * M - position[edge[1]][1], None]

    # Helper function to traverse the tree and collect labels, separating leaf nodes
    def get_labels_and_leaf_status(node, labels=None, is_leaf=None):
        if labels is None:
            labels = []
        if is_leaf is None:
            is_leaf = []

        label_str  = f"{node.label} - {node.index}"
        if node.paired_index is not None:
            label_str += f" - ({node.paired_index})"
        if node.seq:
            label_str += f" - {node.seq}"

        # Determine if the node is a leaf (no children)
        labels.append(label_str)
        # labels.append("")
        is_leaf.append(len(node.children) == 0)  # True for leaf, False for non-leaf

        # Recursively collect labels and leaf status from children
        for child in node.children:
            get_labels_and_leaf_status(child, labels, is_leaf)

        return labels, is_leaf

    # Generate labels and leaf status for coloring
    v_label, leaf_status = get_labels_and_leaf_status(root)

    # Separate node colors based on leaf status
    node_colors = ['#696969' if is_leaf else '#00856A' for is_leaf in leaf_status]
    # node_colors = ['#696969' if is_leaf else '#696969' for is_leaf in leaf_status]

    # Plot using Plotly
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=Xe,
        y=Ye,
        mode='lines',
        line=dict(color='rgb(69,69,69)', width=1),
        # hoverinfo='none'
    ))
    fig.add_trace(go.Scatter(
        x=Xn,
        y=Yn,
        mode='markers',
        name='RNA Structure',
        marker=dict(
            symbol='circle-dot',
            size=15,
            color=node_colors,  # Apply color based on leaf status
            line=dict(color='rgb(50,50,50)', width=1)
        ),
        text=v_label,
        hoverinfo='text',
        opacity=1
    ))
    fig.update_layout(title="RNA Secondary Structure Tree", showlegend=False)

    axis = dict(showline=False,  # hide axis line, grid, ticklabels, and title
                zeroline=False,
                showgrid=False,
                showticklabels=False,
                )

    fig.update_layout(
        title='Tree with Reingold-Tilford Layout',
        annotations=make_annotations(position, v_label),
        font_size=12,
        showlegend=False,
        xaxis=axis,
        yaxis=axis,
        margin=dict(l=40, r=40, b=85, t=100),
        hovermode='closest',
        plot_bgcolor='rgb(255,255,255)'
    )
    fig.show()
