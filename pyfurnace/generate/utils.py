from ..design.core.symbols import Node, dot_bracket_to_tree

def find_stems_in_multiloop(node, stems = None, parent_mloop = None):
    """ Find the stems in a multiloop. A stem is a pair of paired bases that are connected by a loop.
    This function doesn't detect the 0 DT crossovers, because that's 
    doesn't form canonical two m4 crossovers, but one m6 crossovers"""
    if not isinstance(node, Node):
        node = dot_bracket_to_tree(node)

    # Eventually initialize variables
    if stems is None:
        stems = []
    if parent_mloop is None:
        parent_mloop = []

    # A multiloop is a node with at least two paired children
    if node.label == '(' and len([1 for child in node.children if child.label == '(']) > 1:
        # if the parent mloop is not empty, append the last child index to the dovetails
        if parent_mloop:
            stems.append((parent_mloop[-1], node.index))
        # append the children to the parent mloop and recursively search in the child nodes
        for child in node.children:
            parent_mloop.append(child.index)
            find_stems_in_multiloop(child, stems, parent_mloop)
            # remove the last child index from the parent mloop
            parent_mloop.pop()
    else: # Not a mloop recursively search in the child nodes
        for child in node.children:
            find_stems_in_multiloop(child, stems, parent_mloop)

    return stems