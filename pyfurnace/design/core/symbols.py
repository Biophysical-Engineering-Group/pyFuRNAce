"""
This module contains general utility functions and constants for the pyfurnace package.

It includes:
functions for manipulating RNA strands and sequences,
utilities for RNA structures (dot-bracket notation, pair maps, trees)
folding barriers calculation for canonical co-transcriptional RNA Origami
"""

import random
import warnings
from typing import Optional, Tuple, Dict, List, Union, Iterable
from .basepair import BasePair

###
### USEFUL SETS AND DICTIONARIES
###

#: Set of all accepted nucleotides, including standard bases and extended IUPAC codes.
nucleotides = {
    "A",
    "U",
    "C",
    "G",  # standard nucleotides
    "W",  # A or U
    "M",  # A or C
    "R",  # A or G
    "Y",  # C or U
    "K",  # G or U
    "S",  # G or C
    "D",  # A or G or U
    "H",  # A or C or U
    "V",  # A or C or G
    "B",  # C or G or U
    "N",  # any nucleotide
    "X",  # any nucleotide for external Kissing Loops in ROAD
    "&",  # separator
}

#: Dictionary mapping IUPAC codes to corresponding sets of nucleotides.
#: Source: https://www.bioinformatics.org/sms/iupac.html
iupac_code = {
    "A": {"A"},
    "U": {"U"},
    "C": {"C"},
    "G": {"G"},  # standard nucleotides
    "W": {"A", "U"},  # A or U
    "M": {"A", "C"},  # A or C
    "R": {"A", "G"},  # A or G
    "Y": {"C", "U"},  # C or U
    "K": {"G", "U"},  # G or U
    "S": {"G", "C"},  # G or C
    "D": {"A", "G", "U"},  # A or G or U
    "H": {"A", "C", "U"},  # A or C or U
    "V": {"A", "C", "G"},  # A or C or G
    "B": {"C", "G", "U"},  # C or G or U
    "N": {"A", "U", "C", "G"},  # any nucleotide
    "X": {"A", "U", "C", "G"},  # any nucleotide for external Kissing Loops in ROAD
    "&": {"&"},  # separator
}

#: Dictionary mapping IUPAC codes to the set of codes they can base pair with.
base_pairing = {
    "A": {"B", "D", "H", "K", "N", "U", "W", "Y"},
    "U": {"A", "B", "D", "G", "H", "K", "M", "N", "R", "S", "V", "W"},
    "C": {"B", "D", "G", "K", "N", "R", "S", "V"},
    "G": {"B", "C", "D", "H", "K", "M", "N", "S", "U", "V", "W", "Y"},
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
    # this symbol is used for external Kissing Loops in ROAD
    "X": {
        "A",
        "B",
        "C",
        "D",
        "G",
        "H",
        "K",
        "M",
        "N",
        "R",
        "S",
        "U",
        "V",
        "W",
        "Y",
        "X",
    },
    "&": {"&"},  # separator
}

#: Dictionary mapping dot-bracket symbols to their corresponding pairing symbols.
db_pairs = {
    "(": ")",
    "[": "]",
    "{": "}",
    "<": ">",
    "A": "a",
    "B": "b",
    "C": "c",
    "D": "d",
    "E": "e",
    "F": "f",
    "G": "g",
    "H": "h",
    "I": "i",
    "J": "j",
    "K": "k",
    "L": "l",
    "M": "m",
    "N": "n",
    "O": "o",
    "P": "p",
    "Q": "q",
    "R": "r",
    "S": "s",
    "T": "t",
    "U": "u",
    "V": "v",
    "W": "w",
    "X": "x",
    "Y": "y",
    "Z": "z",
}

#: Tuple containing all possible pseudoknot symbols used in dot-bracket notation.
all_pk_symbols = (
    "[",
    "]",
    "{",
    "}",
    "<",
    ">",
    "A",
    "B",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "J",
    "K",
    "L",
    "M",
    "N",
    "O",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "U",
    "V",
    "W",
    "X",
    "Y",
    "Z",
    "a",
    "b",
    "c",
    "d",
    "e",
    "f",
    "g",
    "h",
    "i",
    "j",
    "k",
    "l",
    "m",
    "n",
    "o",
    "p",
    "q",
    "r",
    "s",
    "t",
    "u",
    "v",
    "w",
    "x",
    "y",
    "z",
)

#: Set of all accepted symbols including nucleotides, dots, brackets, and special ROAD
#: symbols.
accept_symbol = (
    nucleotides
    | {".", "(", ")"}
    | set(all_pk_symbols)
    | {
        "─",
        "│",
        "╭",
        "╮",
        "╰",
        "╯",
        "^",
        "*",
        "┼",
        "┊",
        "~",
        "◦",
        "↑",
        "↓",
        "⊗",
        "⊙",
        "•",
        "▂",
        "▄",
        "█",
        "-",
        "|",
        "+",
        ":",
        "=",
        "!",
        " ",
        "/",
        "\\",
        "3",
        "5",
        "&",
    }
)

#: Set of symbols representing base pairs.
bp_symbols = {"┊", "=", ":", "!", "*"}

#: Sequence representing the T7 promoter.
T7_PROMOTER = "TAATACGACTCACTATA"

###
### STRING TRANSLATORS
###

#: String translator mapping all pseudoknot symbols to dots.
pseudo_to_dot = str.maketrans("".join(all_pk_symbols), "." * len(all_pk_symbols))

#: String translator to map each dot-bracket opening symbol to its closing counterpart
#: and vice versa.
pair_db_sym = str.maketrans(
    "".join(db_pairs.keys()) + "".join(db_pairs.values()),
    "".join(db_pairs.values()) + "".join(db_pairs.keys()),
)

#: String translator that removes all nucleotides from a string.
nucl_to_none = str.maketrans("", "", "".join(nucleotides))

#: String translator that removes all accepted symbols from a string.
symb_to_none = str.maketrans("", "", "".join(accept_symbol))

#: String translator that maps nucleotides (including IUPAC codes) to their
#: complementary base.
nucl_to_pair = str.maketrans("AUCGWMRYKSDHVBNX", "UAGCWKYRKSHDBVNX")

#: String translator that removes all non-nucleotide symbols from a string.
only_nucl = str.maketrans("", "", "".join(accept_symbol - nucleotides))

#: String translator for horizontal flip of strand turning symbols.
horiz_flip = str.maketrans("╭╮╰╯/\\", "╮╭╯╰\\/")

#: String translator for vertical flip of strand turning symbols.
verti_flip = str.maketrans("╭╮╰╯/\\", "╰╯╭╮\\/")

#: String translator for 90 degree rotation of strand symbols.
rotate_90 = str.maketrans("╭╮╰╯│|─-", "╮╯╭╰──││")

#: String translator for converting normal symbols to ROAD strand symbols.
symb_to_road = str.maketrans("-|+=:*!", "─│┼┊┊┊┊")

###
### SEQUENCE FUNCTIONS
###


def complement(sequence: str) -> str:
    """
    Return the complement of an RNA sequence using IUPAC base-pairing rules.

    Parameters
    ----------
    sequence : str
        RNA sequence composed of standard and extended IUPAC nucleotide symbols.

    Returns
    -------
    str
        Complementary sequence with each base replaced by its pair.
    """
    return sequence.translate(nucl_to_pair)


def reverse_complement(sequence: str) -> str:
    """
    Return the reverse complement of an RNA sequence.

    Parameters
    ----------
    sequence : str
        RNA sequence to reverse and complement.

    Returns
    -------
    str
        The reverse-complemented sequence.
    """
    return sequence.translate(nucl_to_pair)[::-1]


def pair_nucleotide(nucleotide: str, iupac_symbol: str = "N") -> str:
    """
    Return a nucleotide that can pair with the given nucleotide, optionally
    constrained by an IUPAC symbol.

    Parameters
    ----------
    nucleotide : str
        The base to find a pair for.
    iupac_symbol : str, optional
        Constraint based on allowed base-pairing (IUPAC code).

    Returns
    -------
    str
        A valid pairing nucleotide, or the complement of the input base.
    """
    if iupac_symbol in "AUCG":
        return iupac_symbol
    if iupac_symbol == "K":
        if nucleotide == "G":
            return "U"
        elif nucleotide == "U":
            return "G"
    return nucleotide.translate(nucl_to_pair)


def mutate_nucleotide(
    sequence: str,
    sequence_constraints: str,
    nucl_ind: int,
    pair_map: Dict[int, Optional[int]],
) -> Tuple[Optional[str], Optional[str]]:
    """
    Mutate a nucleotide in the sequence and update its paired counterpart accordingly.

    Parameters
    ----------
    sequence : str
        Original RNA sequence.
    sequence_constraints : str
        A string of IUPAC codes representing base constraints for each position.
    nucl_ind : int
        Index of the nucleotide to mutate.
    pair_map : dict of int to int or None
        A mapping of indices representing base pairs.

    Returns
    -------
    tuple of (str or None, str or None)
        A tuple with the new nucleotide and its paired base if applicable.
    """
    new_set = iupac_code[sequence_constraints[nucl_ind]] - {sequence[nucl_ind]}
    if not new_set:
        return None, None
    new_nucl = random.choice(list(new_set))

    # get the index of the paired nucleotide
    paired_nucl = pair_map[nucl_ind]

    # if the nucleotide is paired, get the paired nucleotide
    if paired_nucl is not None:
        # get the possible nucleotides that can pair with the new nucleotide
        paired_nucl = pair_nucleotide(new_nucl, sequence_constraints[paired_nucl])
    return new_nucl, paired_nucl


def gc_content(seq: str, extended_alphabet: bool = True) -> float:
    """
    Calculate GC content (optionally including ambiguous nucleotides).

    Parameters
    ----------
    seq : str
        The RNA sequence.
    extended_alphabet : bool, optional
        Whether to include ambiguous IUPAC symbols (e.g., S, K) as partial GC.

    Returns
    -------
    float
        Proportion of GC content in the sequence.
    """
    total_count = len(seq)
    if not total_count:
        return 0
    gc_count = seq.count("G") + seq.count("C")
    if extended_alphabet:
        gc_count += (
            seq.count("S")
            + sum(map(seq.count, ["M", "R", "Y", "K", "N", "X"])) / 2
            + sum(map(seq.count, ["D", "H"])) / 3
            + sum(map(seq.count, ["V", "B"])) * 2 / 3
        )
    percentage = gc_count / total_count
    return percentage


###
### STRUCTURE FUNCTIONS
###


def rotate_dot_bracket(dot_bracket: str, shift_left: int) -> str:
    """
    Rotate a dot-bracket structure leftward by a given number of positions.

    Parameters
    ----------
    dot_bracket : str
        RNA secondary structure in dot-bracket notation.
    shift_left : int
        Number of positions to rotate to the left.

    Returns
    -------
    str
        Rotated dot-bracket structure.
    """
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


def dot_bracket_to_pair_map(dot_bracket: str) -> BasePair:
    """
    Convert dot-bracket notation into a pair map.

    Parameters
    ----------
    dot_bracket : str
        Secondary structure in dot-bracket notation.

    Returns
    -------
    BasePair
        A bidirectional dictionary mapping nucleotide indices to their paired
        partner (or None).
    """
    pair_map = BasePair()
    stacks = [[] for _ in range(len(db_pairs))]
    for i, sym in enumerate(dot_bracket):
        if sym in db_pairs.keys():
            stack_num = list(db_pairs.keys()).index(sym)
            stacks[stack_num].append(i)
        elif sym in db_pairs.values():
            stack_num = list(db_pairs.values()).index(sym)
            if stacks[stack_num]:  # pop the last element if the stack is not empty
                j = stacks[stack_num].pop()
                pair_map[j] = i
        else:  # unpaired nucleotide: No pair
            pair_map[i] = None

    # If there are still elements that are not in the dictionary, add them as unpaired
    for i in range(len(dot_bracket)):
        if i not in pair_map:
            pair_map[i] = None

    return pair_map


def pair_map_to_dot_bracket(
    pair_map: Dict[int, Optional[int]], structure_length: Optional[int] = None
) -> str:
    """
    Convert a base pair map into dot-bracket notation.

    Parameters
    ----------
    pair_map : dict
        Mapping from nucleotide index to paired index or None.
    structure_length : int, optional
        Total length of the structure. If None, computed from the map.

    Returns
    -------
    str
        Dot-bracket notation representing the structure.
    """
    if structure_length is None:
        keys_max = max(pair_map.keys(), default=-1)
        values_max = max((k for k in pair_map.values() if k is not None), default=-1)
        structure_length = max(keys_max, values_max) + 1

    ### CREATE THE DOT BRACKET NOTATION ###
    done_pairs = set()

    # Prepare the variables
    dotbracket = ["."] * structure_length
    bracket_count = 0
    recheck = True
    stack = []

    while recheck:  # recheck the structure, for pseudoknots
        recheck = False

        for i in range(structure_length):

            paired = pair_map[i]
            if paired is None or i in done_pairs or paired in done_pairs:
                continue  # ignore unpaired positions or already paired positions

            if paired > i:  # the current position will close the pair later

                # may be a pseudknot:  Check for clash
                if stack and paired > stack[-1]:
                    recheck = True  # clash detected: recheck the structure later
                    continue
                else:  # no clash: add the pair to the stack
                    stack.append(paired)
                    dotbracket[i] = list(db_pairs.keys())[bracket_count]
                    dotbracket[paired] = list(db_pairs.values())[bracket_count]

            else:  # paired already analyzed
                if not stack:  # no pair in the stack
                    continue
                if i == stack[-1]:  # the current position closes the last pair
                    # Remove the pair from the map
                    done_pairs.add(i)
                    done_pairs.add(paired)
                    stack.pop()

        stack = []  # reset the stack count
        bracket_count += 1

        if bracket_count >= 30:
            warnings.warn(
                "Warning: Too many bracket types needed in to write the "
                "structure. Stopping at 30, 'Z' to 'z'.",
                stacklevel=3,
            )
            break

    return "".join(dotbracket)


def dot_bracket_to_stacks(
    dot_bracket: str, only_opening: bool = False
) -> Tuple[str, List[Tuple[int, int]]]:
    """
    Identify contiguous structural elements (stacks) in a dot-bracket string.

    Parameters
    ----------
    dot_bracket : str
        RNA structure in dot-bracket format.
    only_opening : bool, optional
        Whether to return only stacks that start with an opening symbol.

    Returns
    -------
    tuple of str and list of tuple
        - A reduced dot-bracket string marking distinct stacks.
        - A list of (start_index, end_index) tuples for each stack.
    """
    pair_map = dot_bracket_to_pair_map(dot_bracket)
    # a list of tuples containing starting index and the length of the stack
    stacks = []
    # a list of the reduced dot bracket symbols, that will be returned as a string
    reduced_dot_bracket = []

    stack_start = 0
    current_pair_sym = dot_bracket[0]
    last_stack_pair = pair_map[0] + 1 if pair_map[0] else None

    ### NOT NECESSARY
    # open_pair_sym = '.' + ''.join(db_pairs.keys())

    for i, symbol in enumerate(dot_bracket + "_"):

        ### NOT NECESSARY
        # # ignore the symbols that are not '.' or not opening pairs,
        # # they should be already taken into accound
        # if symbol in open_pair_sym:
        #     continue

        # get the paired index of the current symbol
        stack_pair = pair_map.get(i)

        # if symbol change, or the pairing is not consecutive,
        # add the last stack to the list

        if (
            symbol != current_pair_sym
            or stack_pair is not None
            and stack_pair != last_stack_pair - 1
        ):

            add_last_stack = True  # add the last stack to the list by default
            # if only the opening stacks are needed,
            # check if the last stack is an opening stack or unpaired

            if only_opening:
                add_last_stack = (
                    current_pair_sym in db_pairs.keys() or current_pair_sym == "."
                )

            if add_last_stack:
                reduced_dot_bracket.append(current_pair_sym)
                stacks.append((stack_start, i - 1))
            stack_start = i

        last_stack_pair = stack_pair
        current_pair_sym = symbol

    return "".join(reduced_dot_bracket), stacks


class Node:
    """
    A tree node representing a region of an RNA secondary structure.

    Attributes
    ----------
    index : int or None
        Position of the nucleotide in the sequence.
    label : str
        Structure symbol (e.g., '(', ')', '.', '&').
    paired_index : int or None
        Index of the nucleotide this node pairs with, if any.
    parent : Node or None
        Parent node in the tree.
    seq : str or None
        Sequence constraint at this node.
    children : list of Node
        Child nodes representing nested structures.
    """

    def __init__(
        self,
        index: Optional[int] = None,
        label: str = "5'",
        parent: Optional["Node"] = None,
        seq: Optional[str] = None,
        **kwargs,
    ) -> None:
        """
        Initialize a Node object.
        Parameters
        ----------
        index : int or None
            Index of the nucleotide in the sequence.
        label : str
            Structure symbol (e.g., '(', ')', '.', '&').
        parent : Node or None
            Parent node in the tree.
        seq : str or None
            Sequence constraint at this node.
        kwargs : dict
            Additional attributes to set on the node.
        """
        # Initialize the node with the given parameters
        self.index = index
        self.label = label
        self.paired_index = None
        self.parent = parent
        self.seq = seq
        self.children: List[Node] = []
        if parent is not None:
            parent.children.append(self)
        for key, value in kwargs.items():
            setattr(self, key, value)

    def __repr__(self, level=0):
        """
        Return a string representation of the node and its children.
        """
        paired_str = ""
        if self.paired_index is not None:
            paired_str = f" -> {self.paired_index}"

        if level == 0:
            ret = f"{self.label}\n"
        else:
            ret = (
                f"╰{'──' * level} {self.label} - {self.seq} "
                f"- {self.index} {paired_str} \n"
            )
        for child in self.children:
            ret += child.__repr__(level + 1)
        return ret

    def __str__(self):
        """
        Return a string representation of the node.
        """
        return self.__repr__()

    def search(
        self, target_index: Optional[int], target_label: Optional[str] = None
    ) -> Optional["Node"]:
        """
        Recursively search for a node with a matching index and label.

        Parameters
        ----------
        target_index : int or None
            Index to match.
        target_label : str or None
            Optional label to also match.

        Returns
        -------
        Node or None
            The first matching node found, or None if not found.
        """
        # Check if the current node matches the target criteria
        if target_index == self.index and (
            target_label is None or target_label == self.label
        ):
            return self
        elif (
            self.index is None
            and self.parent is None
            and abs(target_index) == float("inf")
        ):
            return self  # special case for the root node

        # Recursively search in child nodes
        for child in self.children:
            result = child.search(target_index=target_index, target_label=target_label)
            if result:
                return result  # Return the node if found

        return None  # Return None if not found


def dot_bracket_to_tree(dot_bracket: str, sequence: Optional[str] = None) -> Node:
    """
    Convert a dot-bracket RNA structure into a hierarchical tree of nodes.

    Parameters
    ----------
    dot_bracket : str
        RNA structure in dot-bracket format. Pseudoknots are ignored.
    sequence : str, optional
        Sequence constraints associated with each index.

    Returns
    -------
    Node
        Root node of the tree representing the full structure.
    """
    # Remove pseudoknots from the dot bracket
    dot_bracket = dot_bracket.translate(pseudo_to_dot)
    if not isinstance(sequence, str) and sequence is not None:
        sequence = str(sequence)
    root = Node()
    current_node = root
    for i, label in enumerate(dot_bracket):
        if sequence is not None:
            seq_const = sequence[i]
        else:
            seq_const = None
        if label == "(":  # add a new node, move to the new node
            new_node = Node(i, label, current_node, seq=seq_const)
            current_node = new_node
        # add the paired index to the parent, move to the parent node
        elif label == ")":
            current_node.paired_index = i
            current_node = current_node.parent
        elif label == ".":  # add a children to the current node, do not move
            new_node = Node(i, label, current_node, seq=seq_const)
        elif label == "&":
            new_node = Node(None, label, current_node)
    return root


def tree_to_dot_bracket(
    node: Node,
    dot_bracket: Optional[List[str]] = None,
    seq_constraints: Union[bool, List[str]] = False,
) -> Union[str, Tuple[str, str]]:
    """
    Convert a tree of Nodes into a dot-bracket structure and optional sequence
    constraints.

    Parameters
    ----------
    node : Node
        Root node of the RNA structure tree.
    dot_bracket : list of str, optional
        Accumulator for the dot-bracket string.
    seq_constraints : bool or list of str
        If True, collect sequence constraints using the `seq` attribute.

    Returns
    -------
    str or tuple of (str, str)
        Dot-bracket string, and sequence constraints if requested.
    """

    if isinstance(seq_constraints, bool) and seq_constraints:
        # initialize the sequence constraints only if it's a bool
        seq_constraints = ["N"]  # assume the sequence is at least one nucleotide long

    if dot_bracket is None:
        dot_bracket = ["."]  # assume the structure is at least one nucleotide long

    if node.parent is not None:  # don't add the root node
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

    for child in node.children:  # recursively add the children
        tree_to_dot_bracket(child, dot_bracket, seq_constraints)

    if node.parent is None:  # we reached back the root node
        if seq_constraints:
            return "".join(dot_bracket), "".join(seq_constraints)

        return "".join(dot_bracket)


###
### DISTANCE FUNCTIONS
###


def hamming_distance(s1: str, s2: str, ignore_ind: Iterable[int] = (), **kwargs) -> int:
    """
    Compute the Hamming distance between two strings, ignoring specified indices.

    Parameters
    ----------
    s1 : str
        First sequence or structure.
    s2 : str
        Second sequence or structure.
    ignore_ind : iterable of int, optional
        Indices to exclude from comparison.

    Returns
    -------
    int
        Number of differing characters at unignored positions.
    """
    return sum(
        (1 for i, (x, y) in enumerate(zip(s1, s2)) if x != y and i not in ignore_ind)
    )


def base_pair_difference(
    s1: str,
    s2: str,
    pair_map1: Optional[Dict[int, Optional[int]]] = None,
    pair_map2: Optional[Dict[int, Optional[int]]] = None,
    ignore_ind: Iterable[int] = (),
    accept_unpaired_ind: Iterable[int] = (),
    **kwargs,
) -> List[int]:
    """
    Find indices where base pairing differs between two structures.

    Parameters
    ----------
    s1 : str
        First structure (dot-bracket notation).
    s2 : str
        Second structure.
    pair_map1 : dict, optional
        Precomputed pair map for s1.
    pair_map2 : dict, optional
        Precomputed pair map for s2.
    ignore_ind : iterable of int, optional
        Indices to ignore in the comparison.
    accept_unpaired_ind : iterable of int, optional
        Indices allowed to be unpaired in s2 without penalty.

    Returns
    -------
    list of int
        Indices of base pairs that differ.
    """
    if pair_map1 is None:
        pair_map1 = dot_bracket_to_pair_map(s1)
    if pair_map2 is None:
        pair_map2 = dot_bracket_to_pair_map(s2)

    # Determine if an index should be considered based on the unpaired indices
    def check_ind(main_pair_map, i, accept_unpaired=False):
        if i in ignore_ind or main_pair_map[i] is None:
            return False
        elif accept_unpaired and i in accept_unpaired_ind:
            # Only count if paired in both but different
            return pair_map2[i] is not None and pair_map1[i] != pair_map2[i]
        return pair_map1[i] != pair_map2[i]

    # pairs in s1 that are not in s2
    diff_set = [i for i in pair_map1 if check_ind(pair_map1, i, True)]
    # pairs in s2 that are not in s1
    diff_set.extend(i for i in pair_map2 if check_ind(pair_map2, i))

    return diff_set


def base_pair_distance(
    s1: str,
    s2: str,
    pair_map1: Optional[Dict[int, Optional[int]]] = None,
    pair_map2: Optional[Dict[int, Optional[int]]] = None,
    ignore_ind: Iterable[int] = (),
    accept_unpaired_ind: Iterable[int] = (),
    **kwargs,
) -> int:
    """
    Compute the number of base pair differences between two structures.

    Parameters
    ----------
    s1 : str
        First dot-bracket structure.
    s2 : str
        Second dot-bracket structure.
    pair_map1 : dict, optional
        Pair map for s1.
    pair_map2 : dict, optional
        Pair map for s2.
    ignore_ind : iterable of int, optional
        Indices to ignore.
    accept_unpaired_ind : iterable of int, optional
        Positions where being unpaired in s2 is allowed.

    Returns
    -------
    int
        Total number of mismatched base pairs.
    """
    return len(
        base_pair_difference(
            s1, s2, pair_map1, pair_map2, ignore_ind, accept_unpaired_ind, **kwargs
        )
    )


###
### FOLDING BARRIERS
###


def folding_barriers(structure: str, kl_delay: int = 150) -> Tuple[str, int]:
    """
    Compute the folding barriers for a given RNA secondary structure.
    This function is based on
    ROAD: https://www.nature.com/articles/s41557-021-00679-1
    This function analyzes the dot-bracket representation of the RNA secondary
    structure to determine folding barriers based on base pair topology and
    kinetic delay.

    Parameters
    ----------
    structure : str
        The dot-bracket notation representing the secondary structure of the RNA.
        If folding barriers is called from a Motif or Origami, the structure is
        already provided by the object in the dot-bracket format.
    kl_delay : int, optional, default=150
        The number of nucleotides (nts) of delay before kissing loops (KLs) snap
        closed. A typical realistic value is around 350 (~1 second), while a
        more conservative  setting is 150 by default.

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
    structure = structure.replace("&", "")
    ss_len = len(structure)
    s_map = dot_bracket_to_pair_map(structure)
    barriers = [""] * len(structure)

    # don't count the nucleotides in the 3' terminal position
    terminal_3_bonus = 16

    current_barr_count = 0
    for ind in range(ss_len):
        # get the bracket at the position
        bra = structure[ind]
        # if bra == '&':
        #     barriers[ind] = '&'
        #     continue

        # set the default barrier to '─', no barrier
        barriers[ind] = "─"

        if ind > ss_len - terminal_3_bonus:
            if s_map[ind] is not None:
                barriers[s_map[ind]] = "─"
            continue

        ### unpaired nts are not barriers and reset the current barrier count
        if bra == ".":
            barriers[ind] = "─"
            current_barr_count = 0

        ### opening pairs may be barriers or not,
        # indicate them with '▂', reset the barrier count
        elif bra == "(":
            barriers[ind] = "▂"
            current_barr_count = 0

        ### closing pairs may be barriers or not,
        # depending on the topology count
        elif bra == ")":
            # if the opening pair is not blocked,
            # we close and reset the current barrier count
            if barriers[s_map[ind]] == "▂":
                barriers[ind] = "─"
                barriers[s_map[ind]] = "─"
                current_barr_count = 0

            # if the closing pair is already blocked, then
            # we mark the barrier reached and start counting
            elif barriers[s_map[ind]] == "█":
                # if the current barrier count is greater than 5, wa mark it with '█'
                if current_barr_count > 5:
                    barriers[ind] = "█"
                    barriers[s_map[ind]] = "▄"
                # if the current barrier count is less than 5, we mark it with '▄'
                else:
                    barriers[ind] = "▄"
                    barriers[s_map[ind]] = "▄"
                current_barr_count += 1

        # if the index is bigger than the kl_delay,
        # we check if the internal KLs are closed
        if ind > kl_delay:

            # we check if the KLs are closed,
            # if yes, we mark them with 'N'
            close_sym = structure[ind - kl_delay]
            if close_sym not in "(.)" and close_sym in db_pairs.values():
                barriers[ind - kl_delay] == "─"
                barriers[s_map[ind - kl_delay]] == "─"

                # for each nt in the delay, we check if there are any barriers,
                # if yes, we mark them with '█'
                for k in range(s_map[ind - kl_delay], ind - kl_delay):
                    if barriers[k] == "▂":
                        barriers[k] = "█"

    penalty = {"█": 2, "▄": 1, "▂": 1, "─": 0}  # , '&': 0}
    repl = [penalty[i] for i in barriers]
    score = sum(repl)
    return "".join(barriers), score


###
### CUSTOM EXCEPTIONS
###


class MotifStructureError(Exception):
    """
    Exception raised when a motif's secondary structure is invalid or inconsistent.

    This error typically indicates issues such as:
    - Overlapping or conflicting strands
    - Joining strands with different directionalities
    - Logical inconsistencies in structural motifs
    """

    pass


class AmbiguosStructure(Warning):
    """
    Warning raised when a structure is ambiguous or potentially problematic.

    This can occur in cases where:
    - The strand routing is not clear

    Parameters
    ----------
    message : str
        Description of the ambiguity or issue.
    """

    def __init__(self, message: str) -> None:
        self.message: str = message

    def __str__(self) -> str:
        return repr(self.message)
