from typing import List, Optional, Tuple, Union
from ..design.core.symbols import Node, dot_bracket_to_tree


def find_stems_in_multiloop(
    node: Union[str, Node],
    stems: Optional[List[Tuple[int, int]]] = None,
    parent_mloop: Optional[List[int]] = None,
) -> List[Tuple[int, int]]:
    """
    Recursive function to find stems that connect through multiloops in a given
    RNA secondary structure.

    A stem is defined as a pair of base-paired regions connected via a multiloop.
    This function does not detect 0-distance crossovers (e.g., m6),
    only canonical m4 crossover patterns.

    Parameters
    ----------
    node : str or Node
        Either a dot-bracket string or a pre-parsed structure tree (`Node`).
    stems : list of tuple of int, optional
        List to accumulate detected stems. Initialized automatically if None.
    parent_mloop : list of int, optional
        List of indices involved in the current multiloop context.

    Returns
    -------
    list of tuple of int
        A list of (i, j) pairs representing stems found across multiloops.
    """
    if not isinstance(node, Node):
        node = dot_bracket_to_tree(node)

    # Eventually initialize variables
    if stems is None:
        stems = []
    if parent_mloop is None:
        parent_mloop = []

    # A multiloop is a node with at least two paired children
    if len([1 for child in node.children if child.label == "("]) > 1:
        # if the parent mloop is not empty,
        # append the last child index to the dovetails
        if parent_mloop:
            stems.append((parent_mloop[-1], node.index))

        # append the children to the parent mloop
        # and recursively search in the child nodes
        for child in node.children:
            parent_mloop.append(child.index)

            find_stems_in_multiloop(child, stems, parent_mloop)
            # remove the last child index from the parent mloop
            parent_mloop.pop()

    else:  # Not a mloop recursively search in the child nodes
        for child in node.children:
            find_stems_in_multiloop(child, stems, parent_mloop)

    return stems
