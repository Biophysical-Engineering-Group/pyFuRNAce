import copy
from pathlib import Path
from functools import wraps
from typing import List, Tuple, Union, Literal, Callable, Optional
from .symbols import (
    iupac_code,
    rotate_dot_bracket,
    pair_map_to_dot_bracket,
    dot_bracket_to_pair_map,
    Node,
    tree_to_dot_bracket,
    dot_bracket_to_tree,
    dot_bracket_to_stacks,
    folding_barriers,
)
from .position import Position, Direction
from .callback import Callback
from .sequence import Sequence
from .strand import Strand
from .basepair import BasePair
from .motif import Motif


class Origami(Callback):
    """
    A class for building and manipulating RNA Origami structures.

    The Origami class organizes RNA `Motif` objects into a 2D matrix to
    represent spatial arrangements of strands. It supports stacking motifs
    horizontally and vertically, calculating vertical junctions and connections,
    assembling the full origami structure in a motif object, and exporting
    the full structure and sequence.

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
    align : str
        Current vertical alignment mode ('left', 'first', 'center').
    assembled : Motif
        Combined Motif representing the full assembled Origami.
    num_char : List[int]
        Number of characters per line, used for alignment.
    num_lines : int
        Number of horizontal lines (rows) in the Origami.
    num_motifs : int
        Total number of motifs in the Origami.
    pair_map : dict
        Dictionary of paired nucleotide indices (alternative to dot-bracket notation).
    pos_index_map : Dict[Position, Tuple[int, int]]
        Map from character position (x, y) in the assembled motif to the
        corresponding index in the original matrix (y, x).
    index_shift_map : Dict[Position, Tuple[int, int]]
        Map from motif matrix indexes (y, x) to spatial shifts (x, y).
    pseudoknots : dict
        Pseudoknot metadata including indices and energies.
    sequence : Sequence
        Full nucleotide sequence of the Origami.
    seq_positions : Tuple[Position]
        The positions of each nucleotide in the origami sequence (x,y coordinates).
        Same as calling origami.assembled.seq_positions. Always 5' to 3'.
    ss_assembly : bool
        Whether to assemble the 3D structure of the origami without locking
        the coordinates of the motifs.
    strands : List[Strand]
        List of individual strands in the Origami.
    structure : str
        Dot-bracket notation of the RNA secondary structure.

    See Also
    --------
    Motif, Strand, Sequence
    """

    @classmethod
    def from_structure(
        cls,
        structure: Optional[Union[str, dict, BasePair, Node]] = None,
        sequence: Optional[str] = None,
        pk_energy=-8.5,
        pk_denergy=0.5,
        motif_list: Optional[List[Motif]] = None,
        **kwargs,
    ) -> "Origami":
        """
        Parse a structure or sequence representation to an Origami object.
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
        motif_list : List[Motif], optional
            A list of specific motifs to parse the structure.
            By default, the motifs are stems and aptamers.
        **kwargs : dict
            Additional arguments to pass to the Motif constructor.


        Returns
        -------
        Origami
            The Origami object created from the structure representation.
        """
        from RNA import fold
        from ..motifs import Stem, aptamers, aptamers_list, Loop
        from ..utils import vertical_double_link, stem_cap_link

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

        if motif_list and not all(isinstance(m, Motif) for m in motif_list):
            raise ValueError("The motif_list must contain only Motif objects.")
        else:
            motif_list = []

        for name in aptamers_list:
            motif = aptamers.__dict__[name]()
            # the motif is made of multiple strands,
            #  the flipped version has a different tree and should be checked
            if len(motif.strands) > 1:
                motif_list.append(motif)
                motif_list.append(motif.copy().flip(reorder=True))

            # if the motif is a loop, so for the loop convention is opened
            # on the right side, but the origami is built with the left side open
            elif isinstance(motif, Loop):
                motif_list.append(motif.flip())

            # just add the motif
            else:
                motif_list.append(motif)

        ### Idea in principle: if the structure is folded with ViennaRNA,
        #  fold the aptamers with ViennaRNA too, so you can find them in the tree.
        #  But single aptamers don't fold in the same way when when they are in a
        # structure, so we cannot use this approach.
        # CODE:
        # if folded:
        #     tree_to_mot = {dot_bracket_to_tree(fold(str(m.sequence))[0],
        #                                        sequence=str(m.sequence)): m
        #                         for m in motif_list}
        # else:

        mot_trees = [
            dot_bracket_to_tree(m.structure, sequence=str(m.sequence))
            for m in motif_list
        ]

        # initialize the origami object
        origami = Origami([[]], align="first", ss_assembly=True)
        current_index = [0, 0]
        m_seq = ["", ""]

        def match_subtree(node: Node, motif_node: Node, depth: int = 0) -> bool:
            """
            Recursively checks if motif_root matches the subtree rooted at node.

            Parameters
            ----------
            node : Node
                The root of the subtree to compare.
            motif_node : Node
                The root of the motif tree.
            depth : int, optional
                The current depth in the tree, used to track the last matched node.

            Returns
            -------
            bool
                True if motif_root matches the subtree rooted at node.
            """
            if (
                motif_node.parent  # not valid for the root node
                and node.label != motif_node.label
            ):
                return None

            if (
                motif_node.parent  # not valid for the root node
                and motif_node.seq is not None
                and node.seq is not None
                and not (iupac_code[node.seq] & iupac_code[motif_node.seq])
            ):
                return None

            # sanitize the children comparison
            motif_child = [c for c in motif_node.children if c.label != "&"]
            node_child = [c for c in node.children if c.label != "&"]

            if not motif_child:
                return node, depth  # Leaf node matched

            if len(motif_child) != len(node_child):
                return None

            node_depths = [
                match_subtree(nc, mc, depth=depth + 1)
                for nc, mc in zip(node_child, motif_child)
            ]
            if not all(node_depths):
                return None

            # save the node at the maximum depth
            # exclude invalid nodes and nodes that are unpaired
            # (the successive tree can grow only from paired nodes)
            node, max_depth = max(
                node_depths,
                key=lambda x: x[1] if (x is not None and x[0].label != ".") else 0,
            )

            if depth == 0:
                # if we are at the root node, we return the node
                return node

            return node, max_depth

        def recursive_build_origami(node, insert_at=None, flip=False):
            """
            Recursively build the origami from the tree representation.

            Parameters
            ----------
            node : Node
                The current node in the tree representation.
            insert_at : Tuple[int, int], optional
                The position in the origami to insert the motif.
            flip : bool, optional
                Whether to flip the motif horizontally and vertically.
            current_index : List[int]
                The current index in the origami matrix, used to track the position
                where the next motif should be inserted.
            """

            nonlocal m_seq  # , current_index

            # initialize the variables
            if insert_at is None:
                insert_at = current_index.copy()
            motif = None

            ### BASE CASES: sequence break, stem, unpaired nucleotide
            if node.label == "&":
                return

            ### Check if the current node matches any motif in the motif list
            if node.parent:
                for tree_mot, mot in zip(mot_trees, motif_list):
                    found_node = match_subtree(node.parent, tree_mot)

                    if found_node:
                        # First, flush stems accumulated
                        # this cause problems with bulges before aptamers
                        # but too many edge cases to handle
                        if m_seq[0]:
                            if m_seq[0] == "N" * len(m_seq[0]):
                                stem = Stem(len(m_seq[0]))
                            else:
                                stem = Stem(sequence=m_seq[0])
                                stem[1].sequence = m_seq[1][::-1]
                            origami.insert(insert_at, stem.flip(flip, flip))
                            insert_at[1] += 1
                            current_index[1] += 1
                            m_seq = ["", ""]

                        motif = mot.copy()
                        node = found_node
                        break

            # find the index of the current node in the parent children
            if node.parent is not None:
                child_seq_ind = [c.index for c in node.parent.children]
                n_c_ind = child_seq_ind.index(node.index)

            if motif is not None:  # if a motif was found, insert it
                pass

            elif node.label == "(":
                m_seq[0] += node.seq if node.seq else "N"
                m_seq[1] += sequence[pair_map[node.index]]

                # if the next node is not a stem, create a stem motif
                if (
                    not node.children
                    or len(node.children) > 1
                    or any(c.label != "(" for c in node.children)
                ):
                    # add wobbles in the stem
                    if m_seq[0] == "N" * len(m_seq[0]):
                        motif = Stem(len(m_seq[0]))
                    else:
                        motif = Stem(sequence=m_seq[0])
                        motif[1].sequence = m_seq[1][::-1]

            elif node.label == ".":
                m_seq[0] += node.seq if node.seq else "N"

                # if the next node adjacent node is not unpaired, create a motif
                if (
                    n_c_ind == len(child_seq_ind) - 1
                    or node.parent.children[n_c_ind + 1].label != "."
                ):
                    motif = Motif(
                        Strand(m_seq[0]), Strand("-" * len(m_seq[0]), start=(0, 2))
                    )

            # add the motif and update the current index
            if motif:
                origami.insert(insert_at, motif.flip(flip, flip))
                current_index[1] += 1  # increment the x index
                m_seq = ["", ""]  # reset the motif sequence

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
                        connect_up = stem_cap_link(vflip=True)
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
                                (i, 0), vertical_double_link().shift((shift_x, 0))
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

        ### ADD THE PSEUDOKNOTS ###
        # dictionary with index as key and pseudoknot id as value
        full_map = dict()
        struct = structure.replace("&", "")
        pair_map = dot_bracket_to_pair_map(struct)

        # map the sequence index to the slice
        pos_to_slice = origami.pos_index_map
        seq_positions = origami.seq_positions
        motif_shifts = origami.index_shift_map

        # iterate over the structure
        i = 0
        while i < len(struct):
            new_pk_info = {"id": [], "ind_fwd": [], "E": [], "dE": []}

            # iterate over the subsequences structure
            sym = struct[i]

            # found pseudoknot
            if sym not in ".()" and i not in full_map:
                # get the length of the pseudoknot
                length = 1
                while struct[i + length] == sym:
                    length += 1

                # get the pseudoknot id of get a new one
                if pair_map[i] in full_map:
                    pk_id = full_map[pair_map[i]] + "'"
                else:
                    inds = [k.split("_")[1].strip("'") for k in full_map.values()]
                    pk_id = "1_" + str(int(max(inds, default="-1")) + 1)

                # get the pseudoknot motif, stand an insert offset
                pos = seq_positions[i]
                motif_yx = pos_to_slice[pos]
                shift_yx = motif_shifts[motif_yx]
                motif = origami._matrix[motif_yx[0]][motif_yx[1]]
                original_pos = (pos[0] - shift_yx[0], pos[1] - shift_yx[1])
                strand_ind = next(
                    i for i, s in enumerate(motif) if original_pos in s.seq_positions
                )
                seq_offset = motif[strand_ind].seq_positions.index(original_pos)

                # add the pseudoknot to the motif
                new_pk_info["id"].append(pk_id)
                new_pk_info["ind_fwd"].append((seq_offset, seq_offset + length - 1))
                indices = range(seq_offset + i, seq_offset + i + length)
                # update the full map
                full_map.update({k: pk_id for k in indices})

                new_pk_info["E"].append(pk_energy)
                new_pk_info["dE"].append(pk_denergy)

                # add the pseudoknots info to the strand
                motif[strand_ind].pk_info = new_pk_info
                i += length
            i += 1

        return origami

    def __init__(
        self,
        matrix: Union[Motif, List[Motif], List[List[Motif]]] = None,
        *args: Union[Motif, List[Motif]],
        align: Literal["left", "first", "center"] = "left",
        copy: bool = False,
        ss_assembly: bool = False,
        **kwargs,
    ) -> None:
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
        self._pos_index_map = dict()
        self._index_shift_map = dict()
        self._assembled = None
        self._ss_assembly = bool(ss_assembly)
        self._pseudoknots = None

        # initialize the matrix
        if not matrix:
            matrix = []

        ### CHECK THE MATRXI
        # the matrix is a proper 2D list
        if (
            isinstance(matrix, (list, tuple))
            and all(isinstance(row, (list, tuple)) for row in matrix)
            and all(isinstance(m, Motif) for row in matrix for m in row)
        ):
            pass

        # the matrix is a list of motif
        elif isinstance(matrix, (list, tuple)) and any(
            isinstance(row, Motif) for row in matrix
        ):
            matrix = [matrix]  # make it a 1D list

        # the matrix is a motif
        elif isinstance(matrix, Motif):
            matrix = [[matrix]]

        # unsupported type for matrix
        else:
            raise ValueError(
                f"The matrix variable may only contain lists of motifs"
                f" or motifs, but it contains {type(matrix)}."
            )

        ### check the type of the args variable
        if args:
            # the args contanes lines (lists of motifs), so args it's a matrix
            if all(isinstance(item, (list, tuple)) for item in args):
                # add the lines to the matrix
                matrix.extend(args)
                # the args contains only motifs, so it's a line
            elif all(isinstance(item, Motif) for item in args):
                matrix[-1].extend(args)  # add the motifs to the last line
            else:
                raise ValueError(
                    f"The args variable may only contain lists of"
                    f" motifs or motifs, but it contains {type(args)}."
                )

        ### add the matrix to the object
        if copy:
            # make a copy of the matrix
            self._matrix = [[m.copy() for m in row] for row in matrix]
        else:
            self._matrix = matrix

        # registest the callbacks
        for row in self._matrix:
            for m in row:
                if self._updated_motif not in m._callbacks:
                    m.register_callback(self._updated_motif)

        ### set the alignment type
        if align not in ("left", "first", "center"):
            raise ValueError(
                f'"{align}" is not an accepted value for the'
                "align_type variable. The align_type variable"
                ' can only be: "left", "first" or "center".'
            )
        else:
            self._align = align

    def __str__(self) -> str:
        """Return a string representation of the assmebled origami
        (the origami motif)."""
        return str(self.assembled)

    def __repr__(self):
        """Return a string representation of the origami object, by
        iterating through the matrix and calling the repr method of
        each motif."""
        reprs = ""
        for line in self._matrix:
            for item in line:
                reprs += repr(item) + ", "
            reprs += ";\n"
        return reprs

    def __getitem__(
        self,
        key: Union[
            int, slice, Tuple[int, int], Tuple[slice, slice], Callable[[Motif], bool]
        ],
    ) -> Union[Motif, List[Motif], List[List[Motif]]]:
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
                x = slice(x, x + 1)
            return [m for row in self._matrix[y] for m in row[x]]

        ### 1D slice, return the row
        elif isinstance(key, slice) or isinstance(key, int):
            return self._matrix[key]

        ### Function, return a 2D list of motifs that satisfy the function
        elif hasattr(key, "__call__"):
            return [[m for m in row if key(m)] for row in self._matrix]

        else:
            raise TypeError(
                "Index must be a single int/slice, "
                "a tuple of (row, col) of int/slice, "
                "or a function to screen the motifs, got: "
                f"{key}, of type: {type(key)}"
            )

    def __setitem__(
        self,
        key: Union[
            int, slice, Tuple[int, int], Tuple[slice, slice], Callable[[Motif], bool]
        ],
        value: Union[Motif, List[Motif], List[List[Motif]]],
    ) -> None:
        """
        Set motif(s) at the specified position in the matrix, trying to
        match the shape of the value to the shape of the key. The value
        is always copied when added to the matrix to avoid references
        problem when repeating motifs in the matrix.
        If the value is a single of motifs or a list of motifs, it will be
        set in the selected row(s). If the value is a 2D list of motifs, it will be
        set in the selected region of the matrix only if the selected region is a
        2D region.

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
        elif isinstance(value, list) and all(isinstance(item, Motif) for item in value):
            value_dimension = 1

        # value is a 2D list of motifs
        elif (
            isinstance(value, list)
            and all(isinstance(item, (list, tuple)) for item in value)
            and all(isinstance(m, Motif) for sublist in value for m in sublist)
        ):
            value_dimension = 2

        else:
            raise ValueError(
                f"Only motifs, lists of motifs, or 2D lists "
                f"of motifs can be added to the Origami, but "
                f"the object {value} was added."
            )

        ### CHECK THE DIMENSIONALITY OF THE KEY

        mask = None
        # if the key is a function, we need to create a mask wich is a submatrix
        # with the slices of the motifs that satisfy the function
        if hasattr(key, "__call__"):
            mask = [
                [(i, slice(j, j + 1)) for j, m in enumerate(row) if key(m)]
                for i, row in enumerate(self._matrix)
            ]

        # the key is a single int, we select a row
        elif isinstance(key, int):
            key_dimension = 1
            y_int = key % len(self._matrix)

            y_slice = slice(y_int, y_int + 1)
            x_slice = slice(0, len(self._matrix[key]))

        # the key is a tuples of ints, we select a single motif
        elif isinstance(key, (tuple, list)) and all(isinstance(i, int) for i in key):
            key_dimension = 0
            # Convert the keys to a positive integer
            y_int = key[0] % len(self._matrix)
            x_int = key[1] % len(self._matrix[y_int])

            # Convert the keys to slices
            y_slice = slice(y_int, y_int + 1)
            x_slice = slice(x_int, x_int + 1)

        # the key is a slice of a row and slice of a column
        # so this is still a 1D region
        elif (
            isinstance(key, (tuple, list))
            and isinstance(key[0], int)
            and isinstance(key[1], slice)
        ):
            key_dimension = 1  # select a row
            y_int = key[0] % len(self._matrix)

            y_slice = slice(y_int, y_int + 1)
            x_slice = key[1]  # get the slice

        # key selects a submatrix, so this is a 2D region
        elif isinstance(key, (tuple, list)) and all(isinstance(i, slice) for i in key):
            key_dimension = 2  # select a 2D region
            y_slice = key[0]
            x_slice = key[1]

        # special case: vertical selection (theoretically a 1D region,
        # but for code purposes we need to treat it as a 2D region)
        elif (
            isinstance(key, (tuple, list))
            and isinstance(key[0], slice)
            and isinstance(key[1], int)
        ):
            key_dimension = 2  # still select a 2D region, but VERTICAL
            x_int = key[1] % len(self._matrix[key[0]])

            y_slice = key[0]
            x_slice = slice(x_int, x_int + 1)

        else:
            raise TypeError(
                "Origami indexes can be: \n"
                "\t - a function to screen the motifs, \n"
                "\t - an int/slice to select a row, \n"
                "\t - a tuple of two int/slice to select a region. \n"
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
                                f"Error while setting the value to the Origami. "
                                f"The lists do not match. Origami indexes: y: {i}, "
                                f"x: {j}."
                            ) from e

        ### REDUCED ALL CASES TO A 2D SLICING
        else:
            for i, line in enumerate(self._matrix[y_slice]):

                # dimensionality 0 or 1: just set the value
                if value_dimension in (0, 1):
                    line[x_slice] = [m.copy() for m in value]

                # dimensionality 2 math the two submatrices
                elif key_dimension == 2 and value_dimension == 2:
                    # try to set the value matching the indices
                    try:
                        line[x_slice] = [m.copy() for m in value[i]]
                    except IndexError as e:
                        raise IndexError(
                            f"Error while setting the value to the Origami. "
                            f"The lists do not match. Origami indexes: "
                            f"y index: {i}, x slice: {x_slice}."
                        ) from e

        ### update the motif
        self._updated_motif()

    def __len__(self):
        """Get the number of rows in the origami."""
        return len(self._matrix)

    def __add__(self, other: "Origami") -> "Origami":
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
            raise TypeError(
                "Unsupported operand type(s) for +: "
                "'Origami' and '{type(other).__name__}'"
            )

        # add extra rows to the new matrix
        diff_len = len(other._matrix) - len(new_matrix)
        if diff_len > 0:
            new_matrix.extend([[] for _ in range(diff_len)])

        for i, row in enumerate(other._matrix):
            new_matrix[i].extend([m.copy() for m in row])

        return Origami(
            new_matrix, align=self.align, ss_assembly=self.ss_assembly, copy=False
        )

    def __bool__(self) -> bool:
        """Return False there are no motifs or all motifs are empty."""
        if not self._matrix:
            return False
        for row in self:
            for motif in row:
                if motif:
                    return True
        return False

    ###
    ### PROPERTIES
    ###

    @property
    def align(self) -> Literal["left", "first", "center"]:
        """
        The alignment type of the rows of the origami.
        """
        return self._align

    @align.setter
    def align(self, new_align):
        """
        Set the alignment type of the rows of the origami.

        Parameters
        ----------
        new_align : {'left', 'first', 'center'}
            The new alignment type for the origami.
            When set to 'left', the motif rows are aligned to the left.
            When set to 'first', the motifs rows are aligned to match the first
            vertical junction.
            When set to 'center', the motifs rows are aligned to the center.

        """

        if new_align not in ("left", "first", "center"):
            raise ValueError(
                f'"{new_align}" is not an accepted value for '
                "the align_type variable. The align_type variable "
                'may only a string reading "left", "first" or "center".'
            )
        self._align = new_align
        self._updated_motif()

    @property
    def assembled(self):
        """
        The matrix of the origami with the motif shifted in the correct position
        for the assembly. The assembled matrix contains rows with the vertical
        connection motifs.
        """
        if self._assembled is None:
            self._assemble()
        return self._assembled

    @property
    def num_char(self) -> List[int]:
        """
        The number of characters in each line of the origami.
        """
        if not self._matrix:
            return 0
        return [sum(m.num_char for m in line) for line in self._matrix]

    @property
    def num_lines(self) -> int:
        """
        The number of lines in the origami.
        """
        if not self._matrix:
            return 0
        return len(self._matrix)

    @property
    def num_motifs(self) -> int:
        """
        The number of motifs in the origami.
        """
        return sum(1 for line in self._matrix for item in line if item is not None)

    @property
    def pair_map(self) -> dict:
        """
        The dictionary of the paired indexes (alternative to the dot bracket
        notation).
        """
        return self.assembled.pair_map

    @property
    def pos_index_map(self) -> dict:
        """
        A dictionary with the symbols position (x, y) as keys and the
        matrix index (y, x) of the motif that contains it as values.
        """
        if self._assembled is None:
            self._assemble()
        return self._pos_index_map

    @property
    def index_shift_map(self) -> dict:
        """
        A dictionary with the slice of the motif in the matrix as key (y, x)
        and positional shift of the motif as values (y, x). The shift is the
        difference between the position of the motif in the matrix and the position
        of the motif in the assembled origami.
        """
        if self._assembled is None:
            self._assemble()
        return self._index_shift_map

    @property
    def positions(self) -> List[Position]:
        """
        The positions of the characters in the origami. Is the same as
        calling the assembled.motif.positions.
        """
        if self._assembled is None:
            self._assemble()
        return self.assembled.positions

    @property
    def pseudoknots(self) -> dict:
        """
        A dictionary with the pseudoknot information.
        The dictionary has pseudoknot IDs as keys and the pseudoknot information as
        values.
        The pseudoknot information is a dictionary with the following keys:
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
        pos_to_ind = {pos: ind for ind, pos in enumerate(self.assembled.seq_positions)}

        def add_pk(strand, pk_index, info_nr, shift, forward=True):
            """Add the pseudoknot information to the pk_dict."""
            # get the pseudoknot information
            pk_info = strand.pk_info

            # add the pseudoknot information to the pk_dict
            pk_dict.setdefault(
                pk_index, {"ind_fwd": [], "ind_rev": [], "E": [], "dE": []}
            )
            pk_dict[pk_index]["E"].append(pk_info["E"][info_nr])
            pk_dict[pk_index]["dE"].append(pk_info["dE"][info_nr])

            # indicate the index of the pseudoknot in the sequence
            start_pos_ind = 0
            if strand.directionality == "35":
                start_pos_ind = -1
            pos = strand.seq_positions[start_pos_ind]

            # get the index of the sequence in the strand
            offset_ind = pos_to_ind[(shift[0] + pos[0], shift[1] + pos[1])]

            # get the start and end positions of the pseudoknot
            pk_start, pk_end = pk_info["ind_fwd"][info_nr]
            start_end_tuple = (offset_ind + pk_start, offset_ind + pk_end)
            # add the start and end positions to the pk_dict
            if forward:
                pk_dict[pk_index]["ind_fwd"].append(start_end_tuple)
            else:
                pk_dict[pk_index]["ind_rev"].append(start_end_tuple)

        pk_motifs = []

        ### collect all the motifs with pseudoknot information
        for i, line in enumerate(self._matrix):
            for j, m in enumerate(line):
                if any(hasattr(s, "pk_info") for s in m):
                    pk_motifs.append((i, j))

        ### Iterate through the strands of the motifs with pseudoknot information
        for i, j in pk_motifs:
            m = self._matrix[i][j]
            shift = self.index_shift_map[(i, j)]

            # get pseudoknot IDs from the strands
            pk_strands = [s for s in m if s.pk_info]
            pk_indexes = [pk_id for s in pk_strands for pk_id in s.pk_info["id"]]

            ### Adjust the pk_index for unique pseudoknots
            if any(ind[0] == "0" for ind in pk_indexes):  # new 0 pseudoknot
                current_n_zero = sum(1 if key[0] == "0" else 0 for key in pk_dict)
                pk_index_0 = "0_" + str(current_n_zero + 1)

            # add the pseudoknots
            for strand in pk_strands:
                for info_nr, pk_index in enumerate(strand.pk_info["id"]):
                    reverse = pk_index[-1] == "'"
                    if pk_index[0] == "0":
                        pk_index = pk_index_0
                    elif reverse:
                        pk_index = pk_index[:-1]
                    add_pk(strand, pk_index, info_nr, shift, forward=not reverse)

        # make the average energy and average tolerance
        for pk in pk_dict.values():
            pk["E"] = sum(pk["E"]) / len(pk["E"])
            pk["dE"] = sum(pk["dE"]) / len(pk["dE"])

        self._pseudoknots = pk_dict
        return self._pseudoknots

    @property
    def sequence(self) -> "Sequence":
        """
        The sequence of the origami, as a Sequence.
        """
        return self.assembled.sequence

    @sequence.setter
    def sequence(self, new_seq):
        """
        Set the sequence of the origami.
        """

        # remove the '&' symbol
        new_seq = new_seq.replace("&", "")
        current_seq = self.sequence.replace("&", "")

        if not isinstance(new_seq, (str, Sequence)) or len(new_seq) != len(current_seq):
            raise ValueError(
                f"The new sequence must be a string or a Sequence object"
                f" with the same lenght of the current sequence "
                f"({len(current_seq)}). Got type: {type(new_seq)}; with "
                f"length: {len(new_seq)}, excluding the '&' symbols."
            )

        # adjust the offset if there are multiple strands
        offset = 0

        # read the maps once to avoid triggering the callback and origami assembly
        pos_to_slice = self.pos_index_map
        origami_motif = self.assembled
        motif_shifts = self.index_shift_map

        # iterate over the strands in the origami motif
        for s in origami_motif:

            # a tuple to identify a specific strand in a motif in the
            # origami matrix
            strand_ID = None
            # initialize/reset the current base map
            new_strand_seq = ""

            # iterate over the nucleotides in the strand
            for ind, pos in enumerate(s.seq_positions):

                ### GET THE STRAND ID FOR THIS NUCLEOTIDE
                # get the y, x cooridnates of the motif in the matrix
                motif_yx = pos_to_slice[pos]
                # get the x, y shift of the motif in the origami positions
                shift_yx = motif_shifts[motif_yx]
                # remove the shifts from the position of the base
                original_pos = (pos[0] - shift_yx[0], pos[1] - shift_yx[1])
                # get the motif at the position
                motif = self._matrix[motif_yx[0]][motif_yx[1]]
                # get the strand index of the motif at the base position
                strand_ind = next(
                    i for i, s in enumerate(motif) if original_pos in s.seq_positions
                )
                if strand_ID is None:
                    strand_ID = (motif_yx[0], motif_yx[1], strand_ind)

                ### NEW STRAND ID?
                # update the sequence of the previous strand before
                # moving to the next one

                if strand_ID != (motif_yx[0], motif_yx[1], strand_ind):

                    # get the strand and set the curent base maps
                    strand = self._matrix[strand_ID[0]][strand_ID[1]][strand_ID[2]]
                    strand_dir = 1
                    if strand.directionality == "35":
                        strand_dir = -1
                    strand.sequence = new_strand_seq[::strand_dir]

                    # reset the current base map with the new strand
                    new_strand_seq = ""

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
    def seq_positions(self) -> Tuple[Position]:
        """
        The positions of each nucleotide in the motif sequence (x,y coordinates).
        The sequence has always the directionality 5' to 3'
        """
        return self.assembled.seq_positions

    @property
    def ss_assembly(self) -> bool:
        """
        Boolean indicating if the origami 3d structure
        is assembled without locking the coordinates of the motifs.
        """
        return bool(self._ss_assembly)

    @ss_assembly.setter
    def ss_assembly(self, new_ss_assembly):
        """
        Set the ss_assembly attribute to True or False.
        """
        self._ss_assembly = bool(new_ss_assembly)
        self._updated_motif()

    @property
    def strands(self) -> List[Strand]:
        """
        The strands of the origami.
        """
        return self.assembled.strands

    @property
    def structure(self) -> str:
        """
        The dot-bracket structure of the origami.
        """
        return self.assembled.structure

    ###
    ###  STATIC METHODS
    ###

    @staticmethod
    def _calculate_connections(
        junctions1: dict,
        junctions2: dict,
        directionalities: List[str],
        x_shift: Union[tuple[int, int], Position] = (0, 0),
        start_y: int = 0,
    ) -> Tuple[Motif, Position]:
        """
        Creates the connection between the rows of the origami.

        Parameters
        ----------
        junctions1: dict
            junctions of the first line
        junctions2: dict
            junctions of the second line
        directionalities: list
            the directionalities of the top junctions
        x_shift: tuple
            The x shift of the junctions of the first and second line
        start_y: int
            The y position of the first line

        Returns
        -------
        Tuple[Motif, Position]
            The connection motifs and height of the vertical connections
        """
        ### take the junctions of the two lines
        j1 = [pos[0] + x_shift[0] for pos in junctions1[Direction.DOWN]]
        j2 = [pos[0] + x_shift[1] for pos in junctions2[Direction.UP]]

        # a junction is missing, then no connection
        if not j2 or not j1:
            return Motif(), 0

        # the number of connections is the minimum of the two junctions
        n_connect = min((len(j1), len(j2)))
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
            crossed = {
                i
                for i, x12 in enumerate(pairs)
                if (
                    i not in closed_crossings
                    and (x_min <= x12[0] <= x_max or x_min <= x12[1] <= x_max)
                )
            }

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

            if x1 < x2:  # the first motif is on the left
                strand = (
                    "│" * n_crossings
                    + "╰"
                    + "─" * (x2 - x1 - 1)
                    + "╮"
                    + "│" * (max_crossing - n_crossings)
                )

            elif x1 > x2:  # the first motif is on the right
                strand = (
                    "│" * (max_crossing - n_crossings)
                    + "╯"
                    + "─" * (x1 - x2 - 1)
                    + "╭"
                    + "│" * n_crossings
                )

            else:  # the motifs are on the same position vertically
                strand = "│" * (max_crossing + 1)

            # can add the symbol "^" for retrocompatibility with ROAD
            # instead add arrows for the directionality of the strand
            # if you  do this, increase the max_crossing by 1
            strand += "↑"

            strands.append(
                Strand(
                    strand,
                    directionality=directionalities[ind],
                    start=(x1, start_y),
                    direction=Direction.DOWN,
                )
            )

        # Extra +1 to the max_crossing to add the symbol "^" or "↑"
        connection_height = Position((0, max_crossing + 1 + 1))
        return Motif(strands, join=False), connection_height

    ###
    ### PROTECTED METHODS
    ###

    def _assemble(self) -> List[List[Motif]]:
        """
        Assemble the origami by shifting the motifs in the correct position,
        concatenating the motifs in the lines, and creating the connection motifs.
        """
        ### Screen the matrix to remove the empty motifs
        self._matrix = [[m for m in line if m] for line in self._matrix]

        ### initialize the variables
        motif_lines = []
        shifts = [[Position.zero() for _ in line] for line in self._matrix]
        align_shifts = [Position.zero() for _ in range(self.num_lines)]

        ### Center the rows, can precompute the shift
        if self.align == "center":
            # take the maximum center position of the motifs in all lines
            max_center = max([num_char // 2 for num_char in self.num_char], default=0)
            # shift to the right to align the center of the motifs
            align_shifts = [
                Position((max_center - num_char // 2, 0)) for num_char in self.num_char
            ]

        for ind, line in enumerate(self._matrix):

            # create the line
            mot_line, vh_shifts = Motif.concat(
                line,
                copy=True,
                align=True,
                extend=True,
                return_shifts=True,
                unlock_strands=self._ss_assembly,
                lock_coords=False,
            )

            mot_line.shift(align_shifts[ind])

            shifts[ind] = [
                h + vh_shifts[i] + align_shifts[ind] for i, h in enumerate(shifts[ind])
            ]

            motif_lines.append(mot_line)

        # shift the motif horizontally to align the first junction
        if self._align == "first":
            _, h_shifts = Motif.align(
                motif_lines,
                axis=0,
                return_shifts=True,
            )
            shifts = [[h + h_shifts[i] for h in line] for i, line in enumerate(shifts)]

        ### calculate the junctions
        ind1 = 0
        while ind1 < len(motif_lines) - 1:
            # get the motifs
            top_motif = motif_lines[ind1]
            bot_motif = motif_lines[ind1 + 1]

            # get the junctions of the motifs
            j1 = top_motif.junctions
            j2 = bot_motif.junctions

            # get the directionalities of the top junction
            mot_to_strand = motif_lines[ind1].get_strand_index_map()
            directs = []
            for pos in j1[Direction.DOWN]:
                strand = top_motif[mot_to_strand[pos]]
                if pos == strand.end and strand.end_direction == Direction.DOWN:
                    directs.append(strand.directionality)
                else:
                    directs.append(strand.directionality[::-1])

            # create the connection motifs
            m_connect, _ = self._calculate_connections(j1, j2, directs)
            # intercalate the connections into the motif lines
            motif_lines.insert(ind1 + 1, m_connect)
            ind1 += 2

        # assemble the origami, piece by piece
        mot = Motif()
        for ind, line in enumerate(motif_lines):
            # add the line to the motif
            mot, v_shifts = Motif.concat(
                mot,
                line,
                axis=0,
                copy=False,
                align=False,
                position_based=True,
                return_shifts=True,
                unlock_strands=self._ss_assembly,
                lock_coords=False,
            )
            if len(v_shifts) > 1 and ind % 2 == 0:
                # add the vertical shift to the horizontal shift
                shifts[ind // 2] = [h + v_shifts[1] for h in shifts[ind // 2]]

        self._assembled = mot
        self._index_shift_map = {
            (i, j): shift
            for i, line in enumerate(shifts)
            for j, shift in enumerate(line)
        }

        self._pos_index_map = {
            (pos + shifts[i][j]): (i, j)
            for i, line in enumerate(self._matrix)
            for j, m in enumerate(line)
            for pos in m.positions
        }

        for s in self._assembled:
            if s.directionality == "35":
                s.invert()

    def _updated_motif(self, **kwargs) -> None:
        """
        Reset cached motif-derived properties and trigger callbacks.

        Parameters
        ----------
        **kwargs : dict
            Optional keyword arguments passed to callbacks.
        """
        self._assembled = None
        self._pseudoknots = None
        self._trigger_callbacks(**kwargs)

    ###
    ### METHODS
    ###

    def append(self, item: Union[Motif, List[Motif]], copy: bool = True) -> None:
        """
        Append a Motif or a list of Motifs to the end of the matrix.
        If the item is a single Motif, it is appended to the last line of the matrix.
        If the item is a list of Motifs, it is appended as a new line in the matrix.

        Parameters
        ----------
        item : Motif or list of Motif
            The motif(s) to append.
        copy : bool, default=True
            Whether to copy motifs before appending.

        Raises
        ------
        TypeError
            If `item` is not a Motif or a list of Motifs.
        """
        if isinstance(item, Motif):
            if not self._matrix:
                self._matrix.append([])
            if copy:
                item = item.copy()

            # update the callbacks:
            if self._updated_motif not in item._callbacks:
                item.register_callback(self._updated_motif)

            self._matrix[-1].append(item)

        elif isinstance(item, (list, tuple)) and all(
            isinstance(m, Motif) for m in item
        ):
            if copy:
                item = [m.copy() for m in item]

            # update the callbacks:
            for m in item:
                if self._updated_motif not in m._callbacks:
                    m.register_callback(self._updated_motif)

            self._matrix.append(item)
        else:
            raise TypeError(
                f"Only motifs or lists of motifs can be added to the "
                f"Origami, but the object {item} was added."
            )

        self._updated_motif()

    def barrier_repr(
        self,
        kl_delay: int = 150,
        barriers: Optional[str] = None,
        return_list: bool = False,
    ) -> Union[str, List[str]]:
        """
        Overlay folding barrier characters onto the structure visualization.

        Parameters
        ----------
        kl_delay : int, default=150
            Delay parameter for computing folding barriers.
        barriers : str, optional
            Precomputed folding barrier string. If None, it will be recomputed.
        return_list : bool, default=False
            Whether to return the result as a list of lines instead of a string.

        Returns
        -------
        str or list of str
            The annotated structure as a single string or a list of lines.
        """
        motif = self.assembled
        origami_lines = str(self).split("\n")
        if barriers is None:
            barriers = motif.folding_barriers(kl_delay=kl_delay)[0]
        for i, (x, y) in enumerate(motif.seq_positions):
            origami_lines[y] = (
                origami_lines[y][:x] + barriers[i] + origami_lines[y][x + 1 :]
            )
        if return_list:
            return origami_lines
        return "\n".join(origami_lines)

    def copy(self) -> "Origami":
        """
        Create a deep copy of the Origami.

        Returns
        -------
        Origami
            A new instance identical to the current one.
        """
        new = Origami.__new__(Origami)
        new._matrix = [
            [m.copy(callback=new._updated_motif) for m in line] for line in self._matrix
        ]
        new._align = self._align
        new._ss_assembly = self._ss_assembly
        new._assembled = self._assembled.copy() if self._assembled else None
        new._pos_index_map = {k: val for k, val in self._pos_index_map.items()}
        new._index_shift_map = {k: val for k, val in self._index_shift_map.items()}
        new._pseudoknots = copy.deepcopy(self._pseudoknots)

        return new

    def duplicate_line(self, idx: int, insert_idx: Optional[int] = None) -> None:
        """
        Duplicate a line of motifs and optionally insert it elsewhere.

        Parameters
        ----------
        idx : int
            Index of the line to duplicate.
        insert_idx : int, optional
            Line index at which to insert the duplicated line. If None, it will be
            added at the end.

        Raises
        ------
        ValueError
            If the given `idx` is not an integer.
        """
        if not isinstance(idx, int):
            raise ValueError(f"The index must be an integer, but {idx} was given.")
        line = self._matrix[idx]
        new_line = [m.copy(callback=self._updated_motif) for m in line]
        if insert_idx is None:
            insert_idx = len(self._matrix)
        self._matrix.insert(insert_idx, new_line)
        self._updated_motif()

    # inherit the documentation from the function
    @wraps(Motif.folding_barriers)
    def folding_barriers(self, kl_delay: int = 150) -> Tuple[str, int]:
        return self.assembled.folding_barriers(kl_delay=kl_delay)

    def get_motif_at_position(self, position: Tuple[int, int]) -> Motif:
        """
        Get a motif from its position in 2D coordinates.

        Parameters
        ----------
        position : tuple of int
            Global coordinate (x, y) in the assembled structure.

        Returns
        -------
        Motif
            The motif located at the position.

        Raises
        ------
        ValueError
            If the position is not valid or not found in the map.
        """
        if (
            not isinstance(position, (tuple, list))
            or len(position) < 2
            or not all(isinstance(i, int) for i in position)
        ):
            raise ValueError(
                f"The position must be a tuple of two integers,"
                f" but {position} was given."
            )

        position = tuple(position)
        if position not in self.pos_index_map:
            raise ValueError(
                f"The position {position} is not in the map" f" of the origami."
            )

        return self[self.pos_index_map[position]]

    def get_motif_at_seq_index(self, index: int) -> Motif:
        """
        Get the motif that contains the given sequence index.

        Parameters
        ----------
        index : int
            Sequence index in the full assembled sequence.

        Returns
        -------
        Motif
            Motif containing the base at the given index.
        """
        return self[self.get_slice_at_seq_index(index)]

    def get_motif_type(self, motif_type: type) -> List[Motif]:
        """
        Get all motifs in the Origami that match the given type.

        Parameters
        ----------
        motif_type : type
            Motif subclass/type to filter.

        Returns
        -------
        List[Motif]
            All motifs of the specified type.
        """
        return [m for line in self._matrix for m in line if isinstance(m, motif_type)]

    def get_slice_at_seq_index(self, index: int) -> Tuple[int, int]:
        """
        Get matrix coordinates of the motif containing the given sequence index.

        Parameters
        ----------
        index : int
            Index in the full sequence.

        Returns
        -------
        Tuple[int, int]
            (row, column) coordinates in the matrix.

        Raises
        ------
        ValueError
            If the index is not found.
        """
        if not isinstance(index, int) or index >= len(self.sequence):
            raise ValueError(
                f"The sequence index must be an integer lower than "
                f"the length of the sequence ({len(self.sequence)}), "
                f"but {index} (type: {type(index)}) was given."
            )

        # map the sequence index to the slice
        ind_to_slice = {
            ind: self.pos_index_map[seq_pos]
            for ind, seq_pos in enumerate(self.assembled.seq_positions)
        }

        return ind_to_slice.get(index)

    def improve_folding_pathway(self, kl_delay: int = 150) -> "Origami":
        """
        Suggest a better folding pathway by circularly shifting the structure.
        This method attempts to find a better folding pathway by shifting the
        structure to a new position. IMPORTANT: this method is designed for
        simple origami blueprints based on DAE crossovers and may not work
        correctly for different structures.

        Parameters
        ----------
        kl_delay : int, default=150
            Delay parameter for kinetic loop folding.

        Returns
        -------
        Origami
            A new Origami object with an optimized folding pathway.
        """
        # import here to avoid circular imports
        from ..motifs import Stem
        from ..utils import start_end_stem

        # remove the motif that start the Origami
        ori = self.copy()
        start_ind = ori.index(lambda m: "5" in m)
        if start_ind:
            ori.pop(start_ind[0])

        # calculate the folding barriers
        start_barrier = ori.folding_barriers(kl_delay=kl_delay)[1]

        ### Check the folding barriers of starting in each possible stem
        ### of at least 5 bases, assuming the motif is has a length property

        # initialize the structures
        db, stacks = dot_bracket_to_stacks(ori.structure)
        min_bar = start_barrier
        best_middle = 0

        # map the sequence index to the slice
        ind_to_slice = {
            ind: ori.pos_index_map[pos]
            for ind, pos in enumerate(ori.assembled.seq_positions)
        }

        # rotate the dot-bracket structure
        for db, (start, end) in zip(db, stacks):
            if db in "()" and (end - start) > 4:
                middle = (start + end) // 2
                new_strucutre = rotate_dot_bracket(ori.structure, middle)
                _, new_bar = folding_barriers(
                    kl_delay=kl_delay, structure=new_strucutre
                )
                # save the best folding barrier
                if new_bar < min_bar:
                    min_bar = new_bar
                    best_middle = middle

        # replace the starting motif, check the two possible orientations
        for flip in range(2):

            # create a copy of the origami, get the slice and the motif
            ori_copy = ori.copy()
            start_slice = ind_to_slice[best_middle]
            m = ori[start_slice]

            ### IMPORTANT REMARK:
            ### THIS ASSUMENTS THE MOTIF HAS A LENGTH PROPERTY
            stem_1 = Stem(m.length // 2)
            start_end = start_end_stem()
            if flip:
                start_end.flip()
            stem2 = Stem(m.length - stem_1.length)
            ori_copy[start_slice] = [stem_1, start_end, stem2]

            # This is the good origami, save it
            if ori_copy.folding_barriers(kl_delay=kl_delay)[1] == min_bar:
                ori = ori_copy
                break

        return ori

    def index(
        self, condition: Union[Callable[[Motif], bool], Motif]
    ) -> List[Tuple[int, int]]:
        """
        Find the matrix coordinates of motifs that satisfy a given condition.

        Parameters
        ----------
        condition : Callable[[Motif], bool] or Motif
            A function that takes a Motif and returns True if it matches,
            or a Motif instance to match directly.

        Returns
        -------
        List[Tuple[int, int]]
            List of (row, column) indices of matching motifs.

        Raises
        ------
        ValueError
            If `condition` is neither a callable nor a Motif.
        """
        if isinstance(condition, Motif):
            motif = condition

            def condition(m: Motif) -> bool:
                return m == motif

        # if we submit a function, return the matrix filtered by the function
        elif not hasattr(condition, "__call__"):
            raise ValueError(
                f"The condition must be a function or a Motif "
                f"object, but {condition} was given."
            )

        return [
            (y, x)
            for y, line in enumerate(self._matrix)
            for x, m in enumerate(line)
            if condition(m)
        ]

    def insert(
        self,
        idx: Union[int, slice, Tuple[int, int]],
        item: Union[Motif, List[Motif]],
        copy: bool = True,
    ) -> None:
        """
        Insert a Motif or list of Motifs at a specific position.

        Parameters
        ----------
        idx : int, slice, tuple of int
            The index or coordinate at which to insert the item(s).
        item : Motif or list of Motif
            The motif(s) to insert.
        copy : bool, default=True
            Whether to copy motifs before inserting.

        Raises
        ------
        ValueError
            If the index or item is invalid.
        """
        ### check the item variable
        if isinstance(item, (list, tuple)) and all(isinstance(m, Motif) for m in item):
            dimension = 2
            if copy:
                item = [m.copy() for m in item]
            # update the callbacks:
            for m in item:
                if self._updated_motif not in m._callbacks:
                    m.register_callback(self._updated_motif)
        elif isinstance(item, Motif):
            dimension = 1
            if copy:
                item = item.copy()
            # update the callbacks:
            if self._updated_motif not in item._callbacks:
                item.register_callback(self._updated_motif)
        else:
            raise ValueError(
                f"Only motifs or lists of motifs can be added to the "
                f"Origami, but the object type {type(item)} was added."
            )

        ### add the item to the matrix according to the index
        if isinstance(idx, (int, slice)):
            if dimension == 1:
                self._matrix.insert(idx, [item])
            elif dimension == 2:
                self._matrix.insert(idx, item)

        elif isinstance(idx, (tuple, list)) and len(idx) == 2:
            if dimension == 2:
                for i, m in enumerate(item):
                    self._matrix[idx[0]].insert(idx[1] + i, m)
            if dimension == 1:
                self._matrix[idx[0]].insert(idx[1], item)

        else:
            raise ValueError(
                f"Index must be a single index or a tuple of "
                f"(row, col) or a list of two indices. "
                f" Got {idx} instead."
            )

        self._updated_motif()

    def pop(self, idx: Union[int, slice, Tuple[int, int]]) -> Union[Motif, List[Motif]]:
        """
        Remove and return a motif or line of motifs at the given index.

        Parameters
        ----------
        idx : int, slice, tuple of int
            The index or coordinates to remove from.

        Returns
        -------
        Motif or list of Motif
            The removed motif(s).

        Raises
        ------
        ValueError
            If the index is not valid.
        """
        if isinstance(idx, (int, slice)):
            if self._matrix:
                popped = self._matrix.pop(idx)
            else:
                return
        elif isinstance(idx, (tuple, list)) and len(idx) == 2:
            if self._matrix[idx[0]]:
                popped = self._matrix[idx[0]].pop(idx[1])
            else:
                return
        else:
            raise ValueError(
                f"Index must be a single index or a tuple of "
                f"(row, col) or a list of two indices. "
                f" Got {idx} instead."
            )

        # remove the callbacks
        if isinstance(popped, list):
            for m in popped:
                m._clear_callbacks()
        else:
            popped._clear_callbacks()

        self._updated_motif()
        return popped

    def reload(self) -> None:
        """
        Recompute the internal structure and regenerate the assembled motif.
        """
        self._updated_motif()
        self._assemble()

    def remove(self, motif: Motif) -> None:
        """
        Remove a specific motif from the matrix.

        Parameters
        ----------
        motif : Motif
            The motif to remove.
        """
        for line in self._matrix:
            if motif in line:
                # remove the callback
                motif._clear_callbacks()
                line.remove(motif)
                break
        self._updated_motif()

    # inherit the documentation from the function
    @wraps(Motif.save_3d_model)
    def save_3d_model(self, *args, **kwargs) -> Optional[Tuple[str, str]]:
        return self.assembled.save_3d_model(*args, **kwargs)

    def save_fasta(
        self, filename_path: str, return_text: bool = False
    ) -> Optional[str]:
        """
        Save the sequence of the Origami to a FASTA file.

        Parameters
        ----------
        filename_path : str
            Path to the output file.
        return_text : bool, default=False
            If True, return the text instead of saving it to a file.

        Returns
        -------
        Optional[str]
            The FASTA text if return_text is True.
        """
        path = Path(filename_path).with_suffix(".fasta")
        name = path.stem
        text = f">{name}\n" f"{self.sequence}\n" f"{self.structure}\n"

        if return_text:
            return text

        with open(str(path), "w", encoding="utf-8") as f:
            f.write(text)

    def save_text(
        self, filename_path: str, to_road: bool = False, return_text: bool = False
    ) -> Optional[str]:
        """
        Save only the structure part of the Origami to a text file.

        Parameters
        ----------
        filename_path : str
            Path to the output file.
        to_road : bool, default=False
            If True, convert to ROAD-compatible format.
        return_text : bool, default=False
            If True, return the text instead of saving it to a file.

        Returns
        -------
        Optional[str]
            The text if return_text is True.
        """
        path = Path(filename_path).with_suffix(".txt")
        name = path.stem
        text = (
            f">{name}\n"
            f"Sequence:\n{self.sequence}\n"
            f"Structure:\n{self.structure}\n"
            f"Pseudoknots info:\n{self.pseudoknots}\n\n"
            f"Blueprint:\n\n"
            f"{self.to_road() if to_road else str(self)}\n\n"
            f"Folding Barriers:\n\n"
            f"{self.barrier_repr()}\n"
        )

        if return_text:
            return text
        with open(str(path), "w", encoding="utf-8") as f:
            f.write(text)

    def to_road(self) -> str:
        """
        Try to convert the Origami's text representation into ROAD-compatible format.

        Returns
        -------
        str
            ROAD-compatible structure representation.
        """
        ori_str = str(self)
        ori_str = ori_str.replace("↑", "^")
        ori_str = ori_str.replace("↓", "^")
        ori_str = ori_str.replace("│ ┊┊┊┊┊┊ │", "│ ****** │")
        ori_str = ori_str.replace(" ┊┊ ", " !! ")
        ori_str = ori_str.replace(" ┊ ", " ! ")
        return ori_str
