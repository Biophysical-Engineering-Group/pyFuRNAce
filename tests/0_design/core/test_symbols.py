import random
import pytest

from pyfurnace.design.core.symbols import (
    # constants / sets / dicts
    nucleotides,
    iupac_code,
    base_pairing,
    db_pairs,
    all_pk_symbols,
    accept_symbol,
    bp_symbols,
    T7_PROMOTER,
    # translators
    pseudo_to_dot,
    pair_db_sym,
    nucl_to_none,
    symb_to_none,
    nucl_to_pair,
    only_nucl,
    horiz_flip,
    verti_flip,
    rotate_90,
    symb_to_road,
    # sequence funcs
    complement,
    reverse_complement,
    pair_nucleotide,
    mutate_nucleotide,
    gc_content,
    # structure funcs
    rotate_dot_bracket,
    dot_bracket_to_pair_map,
    pair_map_to_dot_bracket,
    dot_bracket_to_stacks,
    # trees
    Node,
    dot_bracket_to_tree,
    tree_to_dot_bracket,
    # distances
    hamming_distance,
    base_pair_difference,
    base_pair_distance,
)

# --------------------
# Constants / containers sanity
# --------------------


def test_constants_basic_sanity():
    # nucleotides include standard + some IUPAC + '&'
    for b in "AUCGWMRYKSVDBHNX&":
        assert b in nucleotides

    # iupac_code maps each symbol to a non-empty set of bases
    for k, v in iupac_code.items():
        assert isinstance(v, set) and v

    # base_pairing has entries for standard bases and some IUPACs
    for b in "AUCGN":
        assert b in base_pairing

    # db_pairs contains canonical bracket pairs
    assert db_pairs["("] == ")"
    assert db_pairs["A"] == "a" and db_pairs["Z"] == "z"

    # all_pk_symbols includes at least [] {} <> and letters
    for s in "[]{}<>AaZz":
        assert s in all_pk_symbols

    # accept_symbol contains '.', parentheses, ROAD glyphs, spaces, etc.
    for s in ".()&─│╭╮╰╯^*┼┊~◦↑↓⊗⊙•▂▄█-+|:=!\\/35":
        assert s in accept_symbol

    # bp_symbols are part of accept_symbol
    assert bp_symbols.issubset(accept_symbol)

    # T7 promoter literal
    assert T7_PROMOTER == "TAATACGACTCACTATA"


# --------------------
# Translators
# --------------------


def test_translators_examples():
    # pseudoknots -> dots
    s = "[{A}a](.)X"
    assert s.translate(pseudo_to_dot) == "......(.)."

    # pair_db_sym swaps open/close
    assert "(".translate(pair_db_sym) == ")"
    assert ")".translate(pair_db_sym) == "("
    assert "A".translate(pair_db_sym) == "a" and "a".translate(pair_db_sym) == "A"

    # nucl_to_none removes nucleotides; leaves other chars
    mixed = "AUCGRYNX&-()."
    stripped = mixed.translate(nucl_to_none)
    assert stripped == "-()."  # only non-nucleotide symbols remain

    # symb_to_none removes accepted symbols; leaves foreign chars
    weird = "AUC-G().  xyz"
    leftovers = weird.translate(symb_to_none)
    # accepted symbols are nukes + many glyphs; only foreign letters might remain
    assert set(leftovers) <= set("xyz")

    # nucl_to_pair complement via translate
    assert "AUGC".translate(nucl_to_pair) == "UACG"
    assert "RYSWKMBDHVN".translate(nucl_to_pair)  # doesn’t error

    # only_nucl removes non-nucleotide symbols
    assert "A-U|C".translate(only_nucl) == "AUC"

    # strand glyph flips/rotation
    assert "╭╮╰╯/\\".translate(horiz_flip) == "╮╭╯╰\\/"
    assert "╭╮╰╯/\\".translate(verti_flip) == "╰╯╭╮\\/"
    assert "╭╮╰╯│|─-".translate(rotate_90) == "╮╯╭╰──││"

    # ASCII to ROAD glyphs
    assert "-|+=:*!".translate(symb_to_road) == "─│┼┊┊┊┊"


# --------------------
# Sequence helpers
# --------------------


def test_complement_and_reverse_complement():
    assert complement("AUGC") == "UACG"
    assert reverse_complement("AUGC") == "GCAU"
    # Ambiguity: R <-> Y, S -> S, N passes through
    assert complement("RYSN") == "YRSN"


def test_pair_nucleotide_rules():
    # If constraint is a concrete base, return it
    assert pair_nucleotide("G", "A") == "A"
    # Special wobble for 'K' constraint
    assert pair_nucleotide("G", "K") == "U"
    assert pair_nucleotide("U", "K") == "G"
    # Otherwise return complement
    assert pair_nucleotide("A", "W") == "U"


def test_mutate_nucleotide_basic(monkeypatch):
    seq = "AUGC"
    cons = "NNNN"  # any base allowed everywhere
    # deterministic: pick smallest in set
    monkeypatch.setattr(random, "choice", lambda xs: sorted(xs)[0])

    # Simple hairpin: 0-3 paired, 1-2 paired
    pair_map = {0: 3, 1: 2, 2: 1, 3: 0}

    new0, paired0 = mutate_nucleotide(seq, cons, 0, pair_map)
    assert new0 in {"U", "C", "G"}  # anything but old 'A'
    assert paired0 in "AUCG"  # valid paired partner

    # If no alternatives in constraint, returns (None, None)
    seq2 = "AAAA"
    cons2 = "AAAA"  # only 'A' allowed
    n, p = mutate_nucleotide(seq2, cons2, 0, pair_map)
    assert (n, p) == (None, None)


def test_gc_content_modes():
    assert gc_content("") == 0
    # pure GC -> 1.0
    assert pytest.approx(gc_content("GGCC", extended_alphabet=False)) == 1.0
    # with ambiguity codes (S counts fully; K/N half; V/B at 2/3, D/H at 1/3)
    seq = "SSKKNNVVBBDDHH"  # mix
    val = gc_content(seq, extended_alphabet=True)
    assert 0.0 < val <= 1.0


# --------------------
# Structure helpers
# --------------------


def test_dot_bracket_pair_map_roundtrip():
    s = "((..)).(.)"
    pm = dot_bracket_to_pair_map(s)
    # basic properties
    pm_inds = set(pm.keys()).union(set(k for k in pm.values() if k is not None))
    assert len(pm_inds) == len(s)
    # round-trip
    s2 = pair_map_to_dot_bracket(pm, len(s))
    assert s2 == s


def test_rotate_dot_bracket_preserves_pairs():
    s = "((..))....(())"
    rot = rotate_dot_bracket(s, 3)  # rotate left by 3
    # Rotating back by len - 3 should return to original
    rot_back = rotate_dot_bracket(rot, len(s) - 3)
    assert rot_back == s


def test_dot_bracket_to_stacks_opening_and_all():
    s = "((..))..((..))"
    reduced_all, stacks_all = dot_bracket_to_stacks(s, only_opening=False)
    reduced_open, stacks_open = dot_bracket_to_stacks(s, only_opening=True)

    # Returns same number of stacks tuples as length of reduced string
    assert len(reduced_all) == len(stacks_all)
    assert len(reduced_open) == len(stacks_open)

    # Opening-only reduced string should contain just '(' and '.' segments
    assert set(reduced_open) <= {"(", "."}
    # Each (start,end) is within range and increasing
    for a, b in stacks_all:
        assert 0 <= a <= b < len(s)


# --------------------
# Tree conversion
# --------------------


def test_dot_bracket_tree_roundtrip_with_seq_constraints():
    s = "(.())."
    seq = "NNNNNN"
    root = dot_bracket_to_tree(s, sequence=seq)

    # Root is a Node; contains children
    assert isinstance(root, Node)
    assert root.children

    # Convert back to dot-bracket
    s_back = tree_to_dot_bracket(root)
    assert s_back == s

    # With sequence constraints
    s_back2, cons = tree_to_dot_bracket(root, seq_constraints=True)
    assert s_back2 == s
    # constraints are 'N' by default here (passed through)
    assert len(cons) == len(s)
    assert set(cons) <= set("NAUCGXRYKMSWBDHV&")


def test_node_search_behaviour():
    s = "(.())."
    root = dot_bracket_to_tree(s)

    # Find first child with index 0 '('
    node0 = root.search(0, "(")
    assert isinstance(node0, Node) and node0.index == 0 and node0.label == "("

    # Find an unpaired dot
    dot_idx = s.index(".")
    node_dot = root.search(dot_idx, ".")
    assert node_dot is not None and node_dot.label == "."


# --------------------
# Distances
# --------------------


def test_hamming_and_base_pair_distances():
    s1 = "((..)).."
    s2 = "(.().).."
    # Hamming across same length
    assert hamming_distance(s1, s2) >= 0

    # Base pair maps
    d = base_pair_difference(s1, s2)
    # symmetric-ish check via the wrapper count
    assert base_pair_distance(s1, s2) == len(d)

    # Ignoring an index can reduce differences
    if d:
        ign = {d[0]}
        d2 = base_pair_difference(s1, s2, ignore_ind=ign)
        assert len(d2) <= len(d)

    # Accepting an unpaired index in s2 should not increase the distance
    d3 = base_pair_difference(
        s1, s2, accept_unpaired_ind=[i for i, ch in enumerate(s2) if ch == "."]
    )
    assert len(d3) <= len(d)
