import pytest
from pyfurnace.design.core import (
    Origami,
    Motif,
    pseudo_to_dot,
    dot_bracket_to_tree,
    dot_bracket_to_pair_map,
)
from pyfurnace.design.motifs import Stem, TetraLoop
from pyfurnace.design.utils import simple_origami


@pytest.fixture
def sample_origami():
    """Fixture to create a sample Origami object with a basic structure."""
    motif1 = TetraLoop()
    motif2 = Stem(2)
    motif3 = Stem(3)
    motif4 = TetraLoop(open_left=True)
    return Origami(matrix=[[motif1, motif2], [motif3, motif4]])


def test_origami_add_hstack_and_flags_preserved():
    # Left operand with flags set
    o1 = Origami([Stem(3)])
    o1.align = "center"
    o1.ss_assembly = True

    # Right operand with multiple motifs
    o2 = Origami([TetraLoop(), Stem(2)])

    # Add (horizontal concatenation)
    o3 = o1 + o2

    # Result is a new Origami and originals are untouched
    assert isinstance(o3, Origami)
    assert o3 is not o1 and o3 is not o2
    assert len(o1[0]) == 1  # unchanged
    assert len(o2[0]) == 2  # unchanged

    # Number of rows is the max of the two; number of motifs is the sum
    assert o3.num_lines == max(o1.num_lines, o2.num_lines) == 1
    assert o3.num_motifs == o1.num_motifs + o2.num_motifs

    # __add__ keeps align and ss_assembly from the left operand
    assert o3.align == o1.align == "center"
    assert o3.ss_assembly is True

    # Type validation
    with pytest.raises(TypeError):
        _ = o1 + 123  # type: ignore


def test_positions_property_matches_assembled_and_sequence_length():
    o = Origami([Stem(4), TetraLoop(), Stem(3)])

    # positions proxies assembled.positions
    pos = o.positions
    assert pos == o.assembled.positions

    # Non-empty and well-formed (x, y) integer tuples
    assert len(pos) > 0
    for p in pos:
        assert isinstance(p, tuple) and len(p) == 2
        assert isinstance(p[0], int) and isinstance(p[1], int)

    # Count of positions equals sequence length (ignoring '&')
    seq_len = sum(len(strand) for strand in o.strands)
    assert len(pos) == seq_len


def make_simple_origami():
    # a tiny, always-available layout: one line with two stems and a loop
    return Origami([Stem(4), TetraLoop(), Stem(3)])


def test_construct_variants_and_properties():
    m = Stem(3)
    # single motif
    o1 = Origami(m)
    assert len(o1) == 1
    assert o1.num_lines == 1
    assert o1.num_motifs >= 1
    assert o1  # truthy

    # list of motifs (becomes one row)
    o2 = Origami([Stem(2), Stem(2)])
    assert len(o2) == 1
    assert o2.num_lines == 1
    assert o2.num_motifs >= 2

    # 2D matrix
    o3 = Origami([[Stem(2)], [Stem(2), TetraLoop()]])
    assert len(o3) == 2
    assert o3.num_lines == 2
    assert o3.num_motifs >= 3

    # align setter
    assert o3.align in {"left", "first", "center"}
    o3.align = "center"
    assert o3.align == "center"
    o3.align = "first"
    assert o3.align == "first"
    with pytest.raises(ValueError):
        o3.align = "weird"  # type: ignore


def test_matrix_getitem_and_setitem_shapes():
    o = Origami([[Stem(3), Stem(4)], [TetraLoop()]])
    # 2D index
    m = o[0, 1]
    assert isinstance(m, Motif)

    # row access
    row0 = o[0]
    assert isinstance(row0, list) and all(isinstance(x, Motif) for x in row0)

    # 2D slice -> submatrix (list[list[Motif]])
    sub = o[0:1, 0:2]
    assert isinstance(sub, list)
    assert isinstance(sub[0], list)
    assert all(isinstance(x, Motif) for x in sub[0])

    # function filter (keeps only stems)
    only_stems = o[lambda mot: isinstance(mot, Stem)]
    assert isinstance(only_stems, list)
    assert any(only_stems)

    # __setitem__: replace a single motif
    o[0, 0] = Stem(5)
    assert isinstance(o[0, 0], Stem)
    assert o[0, 0].length == 5

    # __setitem__: replace a 1D row region with a list of motifs
    o[0, 0:2] = [Stem(3), TetraLoop()]
    assert len(o[0]) >= 2
    assert isinstance(o[0, 0], Stem)
    assert o[0, 0].length == 3
    assert isinstance(o[0, 1], TetraLoop)

    # __setitem__: replace a 2D region with a 2D list
    o[0:2, 0:2] = [[Stem(5), TetraLoop()], [Stem(4), Stem(7)]]
    assert isinstance(o[0, 0], Stem)
    assert o[0, 0].length == 5

    # __setitem__: replace a 2D vertical region
    o[0:2, 0] = [[Stem(2)], [Stem(3)]]

    assert isinstance(o[0, 0], Stem)
    assert o[0, 0].length == 2
    assert isinstance(o[1, 0], Stem)
    assert o[1, 0].length == 3

    with pytest.raises(TypeError, match="Origami indexes can be: \n"):
        o["a key"] = Stem(2)  # type: ignore


def test_append_insert_pop_remove_index_copy_duplicate():
    o = Origami([Stem(3)])
    n_rows_start = len(o)

    # append single motif (adds to last row)
    o.append(TetraLoop())
    assert len(o[0]) >= 2

    # append list -> new row
    o.append([Stem(2), Stem(2)])
    assert len(o) == n_rows_start + 1

    # insert single at position
    o.insert((0, 1), Stem(4))
    assert "Stem" in type(o[0, 1]).__name__

    # insert list at position (extends row)
    o.insert((0, 2), [TetraLoop(), Stem(2)])
    assert len(o[0]) >= 4

    # index() by function
    indices = o.index(lambda m: "TetraLoop" in type(m).__name__)
    assert all(isinstance(t, tuple) and len(t) == 2 for t in indices)
    # index() by motif instance
    target = o[indices[0]]
    indices2 = o.index(target)
    assert indices2 and indices2[0] == indices[0]

    # pop single motif
    popped = o.pop((0, 0))
    assert isinstance(popped, Motif)

    # pop line
    popped_line = o.pop(0)
    assert isinstance(popped_line, list)
    assert all(isinstance(m, Motif) for m in popped_line)

    # remove (by identity)
    o2 = Origami([Stem(3), TetraLoop()])
    victim = o2[0, 1]
    o2.remove(victim)
    assert all(m is not victim for m in o2[0])

    # duplicate_line
    o3 = Origami([Stem(3), TetraLoop()], [Stem(2)])
    n_before = len(o3)
    o3.duplicate_line(0)
    assert len(o3) == n_before + 1

    # copy() should deep-copy motifs/matrix
    oc = o3.copy()
    # mutate original
    o3[0, 0] = Stem(5)
    assert "Stem" in type(oc[0, 0]).__name__  # still a stem


def test_getitem_single_motif(sample_origami):
    """Test __getitem__ for retrieving a single motif."""
    motif = sample_origami[0, 1]
    assert isinstance(motif, Motif), "Expected a Motif object."
    assert isinstance(motif, Stem), "Expected a Stem object."
    assert (
        sample_origami[0, 1] == sample_origami[0][1]
    ), "Expected equivalent results for matrix and row/column indexing."
    assert Stem(2) == motif, "Expected a Stem(2) motif."


def test_getitem_row(sample_origami):
    """Test __getitem__ for retrieving a row."""
    row = sample_origami[0]
    assert isinstance(row, list), "Expected a list of motifs."
    assert len(row) == 2, "Row length mismatch."
    assert all(isinstance(m, Motif) for m in row), "Row contains non-Motif elements."
    assert [TetraLoop(), Stem(2)] == row, "Row content mismatch."


def test_getitem_submatrix(sample_origami):
    """Test __getitem__ for retrieving a submatrix."""
    submatrix = sample_origami[0:2, 0:1]
    assert isinstance(submatrix, list), "Expected a list of lists."
    assert len(submatrix) == 2, "Submatrix row count mismatch."
    assert all(
        isinstance(row, list) for row in submatrix
    ), "Submatrix rows are not lists."
    assert all(
        isinstance(m, Motif) for row in submatrix for m in row
    ), "Submatrix contains non-Motif elements."
    assert [[TetraLoop()], [Stem(3)]] == submatrix, "Submatrix content mismatch."


def test_getitem_function_mask(sample_origami):
    """Test __getitem__ with a function to screen motifs."""
    submatrix = sample_origami[lambda m: isinstance(m, Stem)]
    assert isinstance(submatrix, list), "Expected a list of lists."
    assert len(submatrix) == 2, "Submatrix row count mismatch."
    assert all(
        isinstance(row, list) for row in submatrix
    ), "Submatrix rows are not lists."
    assert all(
        isinstance(m, Motif) for row in submatrix for m in row
    ), "Submatrix contains non-Motif elements."
    assert [[Stem(2)], [Stem(3)]] == submatrix, "Submatrix content mismatch."


def test_setitem_single_motif(sample_origami):
    """Test __setitem__ for setting a single motif."""
    new_motif = Stem(1)
    sample_origami[1, 0] = new_motif
    assert sample_origami[1, 0] == Stem(1), "Motif was not correctly set."


def test_setitem_single_motif_vaule_list(sample_origami):
    """Test __setitem__ for setting a row."""
    new_row = [Stem(1), Stem(2), Stem(3)]
    sample_origami[0, 1] = new_row
    assert sample_origami[0] == [
        TetraLoop(),
        Stem(1),
        Stem(2),
        Stem(3),
    ], "Row was not correctly set."
    assert len(sample_origami[0]) == 4, "Row length mismatch after setting."


def test_setitem_row(sample_origami):
    """Test __setitem__ for setting a row."""
    new_row = [Stem(1), Stem(3)]
    sample_origami[0] = new_row
    assert sample_origami[0] == [Stem(1), Stem(3)], "Row was not correctly set."
    assert len(sample_origami[0]) == 2, "Row length mismatch after setting."


def test_setitem_row_one_value(sample_origami):
    """Test __setitem__ for setting a row."""
    new_row = Stem(1)
    sample_origami[0] = new_row
    assert sample_origami[0] == [Stem(1)], "Row was not correctly set."
    assert len(sample_origami[0]) == 1, "Row length mismatch after setting."


def test_setitem_submatrix(sample_origami):
    """Test __setitem__ for setting a submatrix."""
    new_submatrix = [[Stem(1)], [Stem(2)]]
    sample_origami[0:2, 0:1] = new_submatrix
    assert len(sample_origami[0]) == 2, "Submatrix row count mismatch after setting."
    assert sample_origami[0][0] == Stem(1), "Submatrix element mismatch."
    assert sample_origami[1][0] == Stem(2), "Submatrix element mismatch."


def test_setitem_submatrix_function(sample_origami):
    """Test __setitem__ for setting a submatrix."""
    new_submatrix = Stem(1)
    sample_origami[lambda m: isinstance(m, TetraLoop)] = new_submatrix
    assert len(sample_origami[0]) == 2, "Submatrix row count mismatch after setting."
    assert sample_origami[0][0] == Stem(1), "Submatrix element mismatch."
    assert sample_origami[1][1] == Stem(1), "Submatrix element mismatch."

    # Put back the original tetraloop
    new_submatrix = [[TetraLoop()], [TetraLoop(open_left=True)]]
    sample_origami[lambda m: isinstance(m, Stem) and m.length == 1] = new_submatrix
    assert len(sample_origami[0]) == 2, "Submatrix row count mismatch after setting."
    assert sample_origami[0][0] == TetraLoop(), "Submatrix element mismatch."
    assert sample_origami[1][1] == TetraLoop(
        open_left=True
    ), "Submatrix element mismatch."


def test_invalid_getitem(sample_origami):
    """Test __getitem__ with invalid input."""
    with pytest.raises(TypeError):
        _ = sample_origami["invalid_key"]


def test_invalid_setitem(sample_origami):
    """Test __setitem__ with invalid input."""
    with pytest.raises(ValueError):
        sample_origami[0, 0] = "not_a_motif"


def test_delitem_single_row(sample_origami):
    origami = sample_origami
    # Deleting the first row
    del origami[0]
    assert len(origami) == 1
    assert origami[0] == [Stem(3), TetraLoop(open_left=True)]


def test_delitem_multiple_rows(sample_origami):
    origami = sample_origami
    # Deleting a slice of rows (first row)
    del origami[0:1]
    assert len(origami) == 1
    assert origami[0] == [Stem(3), TetraLoop(open_left=True)]


def test_delitem_single_motif(sample_origami):
    origami = sample_origami
    # Deleting a single motif at (0, 1) (first row, second motif)
    del origami[0, 1]
    assert origami[0] == [TetraLoop()]
    assert len(origami) == 2  # Should still have two rows


def test_delitem_row_slice(sample_origami):
    origami = sample_origami
    # Deleting a slice of the first row (motif2)
    del origami[0, 1:2]
    assert origami[0] == [TetraLoop()]
    assert len(origami) == 2  # Still 2 rows


def test_delitem_2d_submatrix(sample_origami):
    origami = sample_origami
    # Deleting a submatrix (first two rows, columns 1 to 2)
    del origami[0:2, 1:2]
    assert origami[0] == [TetraLoop()]
    assert origami[1] == [Stem(3)]
    assert len(origami) == 2


def test_delitem_column(sample_origami):
    origami = sample_origami
    # Deleting the second column (1st motif of each row)
    del origami[:, 1]
    assert origami[0] == [TetraLoop()]
    assert origami[1] == [Stem(3)]


def test_delitem_callable(sample_origami):
    origami = sample_origami
    # Deleting motifs where the motif is a Stem (2 bases)
    del origami[lambda m: isinstance(m, Stem) and len(m) == 2]
    assert origami[0] == [TetraLoop()]
    assert len(origami) == 2


def test_delitem_invalid_key(sample_origami):
    origami = sample_origami
    with pytest.raises(TypeError):
        # Invalid key type should raise an error
        del origami["invalid"]


def test_index_by_motif(sample_origami):
    """Test the index function with an exact match using a Motif."""
    # Find Stem(2) motif
    assert sample_origami.index(Stem(2)) == [(0, 1)]

    # Find TetraLoop motif (with open_left=True)
    assert sample_origami.index(TetraLoop(open_left=True)) == [(1, 1)]


def test_index_by_type(sample_origami):
    """
    Test the index function when given a Motif subclass (type)
    to match by isinstance.
    """
    # All TetraLoop instances: (0, 0) and (1, 1)
    assert sample_origami.index(TetraLoop) == [(0, 0), (1, 1)]

    # All Stem instances: (0, 1) and (1, 0)
    assert sample_origami.index(Stem) == [(0, 1), (1, 0)]


def test_index_by_type_matrix_format(sample_origami):
    """Test type-based indexing with return_matrix_format=True."""
    # TetraLoop positions in matrix format
    result = sample_origami.index(TetraLoop, return_matrix_format=True)
    # row 0: TetraLoop() at x=0
    # row 1: TetraLoop(open_left=True) at x=1
    assert result == [[0], [1]]

    # Stem positions in matrix format
    result = sample_origami.index(Stem, return_matrix_format=True)
    # row 0: Stem(2) at x=1
    # row 1: Stem(3) at x=0
    assert result == [[1], [0]]


def test_index_type_equivalent_to_lambda(sample_origami):
    """Passing a type should behave like using isinstance in a lambda."""
    # TetraLoop
    by_type = sample_origami.index(TetraLoop)
    by_lambda = sample_origami.index(lambda m: isinstance(m, TetraLoop))
    assert by_type == by_lambda

    # Stem
    by_type = sample_origami.index(Stem)
    by_lambda = sample_origami.index(lambda m: isinstance(m, Stem))
    assert by_type == by_lambda


def test_index_by_base_class(sample_origami):
    """If you pass the base Motif class, you should get all motifs."""
    indices = sample_origami.index(Motif)
    assert sorted(indices) == [(0, 0), (0, 1), (1, 0), (1, 1)]


def test_index_by_condition(sample_origami):
    """Test the index function with a condition (lambda function)."""
    # Find all TetraLoop motifs (without considering open_left)
    assert sample_origami.index(lambda m: isinstance(m, TetraLoop)) == [
        (0, 0),
        (1, 1),
    ]

    # Find all Stem motifs with length 2
    assert sample_origami.index(lambda m: isinstance(m, Stem) and m.length == 2) == [
        (0, 1)
    ]


def test_index_matrix_format(sample_origami):
    """Test the index function with return_matrix_format=True."""
    # Find Stem(2) motif and return in matrix format
    result = sample_origami.index(Stem(2), return_matrix_format=True)
    assert result == [[1], []], "Matrix format result is incorrect."

    # Find TetraLoop motifs and return in matrix format
    result = sample_origami.index(
        lambda m: isinstance(m, TetraLoop), return_matrix_format=True
    )
    assert result == [[0], [1]], "Matrix format result is incorrect."


def test_index_no_matches(sample_origami):
    """Test the index function when no motifs match."""
    # Searching for a motif that does not exist (no Stem of length 4)
    assert (
        sample_origami.index(lambda m: isinstance(m, Stem) and len(m) == 4) == []
    ), "No matches expected."

    # Searching for a non-existent TetraLoop
    fake_loop = TetraLoop()
    fake_loop.sequence = "AAAA"  # Different sequence
    assert sample_origami.index(fake_loop) == [], "No matches expected."


def test_index_empty_matrix():
    """Test the index function with an empty matrix."""
    origami = Origami(matrix=[[], []])  # Empty 2x2 matrix
    assert (
        origami.index(lambda m: isinstance(m, TetraLoop)) == []
    ), "Expected no matches in empty matrix."
    assert origami.index(Stem(2)) == [], "Expected no matches in empty matrix."


def test_index_invalid_condition(sample_origami):
    """Test the index function with invalid condition input."""
    with pytest.raises(ValueError, match="The condition must be a function"):
        sample_origami.index(123)  # Invalid condition type

    with pytest.raises(ValueError, match="The condition must be a function"):
        sample_origami.index(None)  # Invalid condition type


def test_from_structure():
    """Test the from_structure class method."""
    origami = Origami.from_structure(
        sequence=(
            "GCACAGUGCUAUGAGUGUGCACGGGAUCCCGACUGGCCGCAUCGCGAAAGUGGCCAGGUAAC"
            "GAAUGGAUCCUGUGCUGCACAUUAGAGUCGCUGUAUGACCCAUCGCGAAAGGGUCGUACAGCGGCUCUAGUG"
            "UGCUCGCGUGCCUCAGAGGACCUGUCACCAUCGCGAAAGGUGAUAGGUCCUUUGAGGUACGCGUCACUCGUA"
            "GCAUUGUGCCUGUCUCCAUCGCGAAAGGAGAUAG"
        )
    )
    assert origami.sequence == (
        "GCACAGUGCUAUGAGUGUGCACGGGAUCCCGACUGGCCGCAUCGCGAAAGUGGCCAGGUAACGAAUGG"
        "AUCCUGUGCUGCACAUUAGAGUCGCUGUAUGACCCAUCGCGAAAGGGUCGUACAGCGGCUCUAGUGUGCUCG"
        "CGUGCCUCAGAGGACCUGUCACCAUCGCGAAAGGUGAUAGGUCCUUUGAGGUACGCGUCACUCGUAGCAUUG"
        "UGCCUGUCUCCAUCGCGAAAGGAGAUAG"
    )
    assert origami.structure == (
        "(((((((((((((((((.(((((((((((((.((((((((.........))))))))....))...))"
        "))))))))).(((((((((((((((((((((((((.........))))))))))))))))))))))))).(("
        "(((((((((((((((((((((((.........))))))))))))))))))))))))).))))))))))))))"
        ")))((((((((.........))))))))"
    )
    origami2 = Origami.from_structure(
        sequence=(
            "GGGAGAUAUGGGGCUGGCCACGAACCCGAUACGUGGCUAGCGGGCUUUCGAGUCCGAUGCUGACGAACCCG"
            "AUACGUCAGUAUCUCCUGCCAACUUGCCAGGCGGGACAAGAGUAACCGUUCAACUUUUGCCCGUAUCUCCCU"
            "AAUGUGGCUAGGGGUCAAGAACGGAGACUCCUGACUCCAAAGGCAAGAUGGGGUCCACUGGUACGAACCCGA"
            "UACGUACCGGUGCAGCGUUCGCGUUGGCCUUGAACGAACCCGAUACGUUCAGGGCAGCCAUAUUACUGCAAG"
            "AGGAUCCCGACUGGCGAGAGCCAGGUAACGAAUGGAUCCUCUG"
        )
    )
    assert origami2.sequence == (
        "GGGAGAUAUGGGGCUGGCCACGAACCCGAUACGUGGCUAGCGGGCUUUCGAGUCCGAUGCUGACGAAC"
        "CCGAUACGUCAGUAUCUCCUGCCAACUUGCCAGGCGGGACAAGAGUAACCGUUCAACUUUUGCCCGUAUCUC"
        "CCUAAUGUGGCUAGGGGUCAAGAACGGAGACUCCUGACUCCAAAGGCAAGAUGGGGUCCACUGGUACGAACC"
        "CGAUACGUACCGGUGCAGCGUUCGCGUUGGCCUUGAACGAACCCGAUACGUUCAGGGCAGCCAUAUUACUGCA"
        "AGAGGAUCCCGACUGGCGAGAGCCAGGUAACGAAUGGAUCCUCUG"
    )
    assert origami2.structure == (
        "((((((((((((((((((((((.........))))))))))(((((....)))))((((((((((..."
        "......))))))))))(((((((.........)))))))(((((((.........)))))))))))))))))"
        "))(((((((((((((((((.........)))))))(((((((.........)))))))((((((((((...."
        ".....))))))))))(((((....)))))((((((((((.........))))))))))))))))))))...."
        ".(((((((((((.(((((....)))))....))...)))))))))."
    )

    origami = simple_origami([-2], align="center")
    assert len(origami.sequence) > 0
    assert len(origami.structure) > 0
    origami_from_str = Origami.from_structure(
        structure=origami.structure, sequence=origami.sequence
    )
    assert origami_from_str.sequence == origami.sequence
    assert origami_from_str.structure == origami.structure.translate(pseudo_to_dot)
    pk_from_str = origami_from_str.pseudoknots
    pk_values = list(origami.pseudoknots.values())
    for pk_dict1 in pk_from_str.values():
        hit = False
        for pk_dict2 in pk_values:
            if pk_dict1["ind_fwd"] == pk_dict2["ind_fwd"]:
                assert pk_dict1["ind_rev"] == pk_dict2["ind_rev"]
                hit = True
                break
            elif pk_dict1["ind_rev"] == pk_dict2["ind_fwd"]:
                assert pk_dict1["ind_fwd"] == pk_dict2["ind_rev"]
                hit = True
                break
        assert hit, "Pseudoknot not found"

    db = "...(...)....(((...)))"
    ori_from_db = Origami.from_structure(db)
    ori_from_pm = Origami.from_structure(dot_bracket_to_pair_map(db))
    ori_from_tree = Origami.from_structure(dot_bracket_to_tree(db))
    assert ori_from_db.structure == db
    assert ori_from_pm.structure == db
    assert ori_from_tree.structure == db

    # from structure with a split sequence
    seq = "GGGGGG&CCCCCCA"
    ori_from_seq = Origami.from_structure(sequence=seq)
    assert ori_from_seq.sequence == seq
    assert ori_from_seq.structure == "((((((&))))))."


def test_assembled_string_sequence_structure_are_accessible():
    o = make_simple_origami()
    # __str__ should render something
    s = str(o)
    assert isinstance(s, str) and s.strip() != ""
    assert repr(TetraLoop()) in repr(o)

    # sequence / structure accessible and lengths match (ignoring '&')
    seq = str(o.sequence).replace("&", "")
    struct = str(o.structure).replace("&", "")
    assert isinstance(seq, str) and isinstance(struct, str)
    assert len(seq) == len(struct) and len(seq) > 0


def test_sequence_setter_roundtrip():
    o = make_simple_origami()
    # new sequence must match length (ignoring '&')
    L = len(str(o.sequence).replace("&", ""))
    new_seq = "A" * L
    # setter accepts str
    o.sequence = new_seq
    assert str(o.sequence).replace("&", "") == new_seq


def test_save_text_and_fasta_return_text_only(tmp_path):
    o = make_simple_origami()
    # return_text=True -> just return strings, no file I/O required
    fasta_text = o.save_fasta(tmp_path / "foo.fasta", return_text=True)
    assert fasta_text.startswith(">foo")
    assert "\n" in fasta_text and len(fasta_text.splitlines()) >= 3

    txt = o.save_text(tmp_path / "bar.txt", to_road=False, return_text=True)
    # contains sections and blueprint
    assert ">bar" in txt
    assert "Sequence:" in txt
    assert "Structure:" in txt
    assert "Pseudoknots info:" in txt
    assert "Blueprint:" in txt
    assert "Folding Barriers:" in txt

    # When writing to disk (no return_text), files exist
    out_fa = tmp_path / "x.fasta"
    o.save_fasta(out_fa)
    assert out_fa.exists()
    out_txt = tmp_path / "y.txt"
    o.save_text(out_txt)
    assert out_txt.exists()


def test_to_road_replacements(monkeypatch):
    o = make_simple_origami()

    # Monkeypatch __str__ to include all special symbols the method replaces
    mock_body = "│ ┊┊┊┊┊┊ │  ┊┊  ┊  ↑  ↓"
    monkeypatch.setattr(Origami, "__str__", lambda self: mock_body, raising=True)

    out = o.to_road()
    # caret replaces arrows; spacer replacements
    assert "^" in out
    assert "******" in out  # for the │ ┊┊┊┊┊┊ │ replacement
    assert " !! " in out  # double spacer
    assert " ! " in out  # single spacer


def test_reload_and_ss_assembly_toggle():
    o = make_simple_origami()
    o.ss_assembly = not o.ss_assembly
    o.reload()  # should re-assemble without error
    now = str(o)
    assert isinstance(now, str) and now


def test_barrier_repr_shapes():
    o = make_simple_origami()
    lines = str(o).splitlines()
    out_list = o.barrier_repr(return_list=True)
    assert isinstance(out_list, list)
    assert len(out_list) == len(lines)
    # As a single string too
    out_str = o.barrier_repr(return_list=False)
    assert isinstance(out_str, str)
    assert out_str.count("\n") == len(lines) - 1


def test_lookup_helpers():
    o = make_simple_origami()
    # pick a valid sequence index somewhere in the middle
    seq_len = len(str(o.sequence).replace("&", ""))
    mid = seq_len // 2
    yx = o.get_slice_at_seq_index(mid)
    assert isinstance(yx, tuple) and len(yx) == 2
    m = o.get_motif_at_seq_index(mid)
    assert isinstance(m, Motif)

    # get_motif_type
    stems = o.get_motif_type(type(Stem(2)))
    assert isinstance(stems, list)
    assert all("Stem" in type(x).__name__ for x in stems)

    # get_motif_at_position: pick the first seq position
    pos0 = o.seq_positions[0]
    m2 = o.get_motif_at_position(pos0)
    assert isinstance(m2, Motif)

    # invalid inputs
    with pytest.raises(ValueError):
        o.get_motif_at_position(("x", "y"))  # type: ignore
    with pytest.raises(ValueError):
        o.get_slice_at_seq_index(10**9)  # out of range


def test_insert_variants():
    # Start with a 2x? matrix
    o = Origami([[Stem(3), TetraLoop()], [Stem(2)]])

    # 1) Insert a single motif at specific (row, col)
    before_row_len = len(o[0])
    o.insert((0, 1), Stem(4))  # insert between Stem(3) and TetraLoop()
    assert len(o[0]) == before_row_len + 1
    assert isinstance(o[0, 1], Motif)

    # 2) Insert a list of motifs at specific (row, col) (extends the row)
    prev_len = len(o[0])
    o.insert((0, 2), [TetraLoop(), Stem(2)])
    assert len(o[0]) == prev_len + 2
    # inserted order preserved
    assert isinstance(o[0, 2], Motif)
    assert isinstance(o[0, 3], Motif)

    # 3) Insert a single motif at row index (creates a NEW line at index 1)
    start_rows = len(o)
    o.insert(1, Stem(5))
    assert len(o) == start_rows + 1
    # the new row contains exactly one motif we inserted
    assert len(o[1]) == 1
    assert isinstance(o[1][0], Motif)

    # 4) Insert a whole row (list of motifs) at row index (creates a NEW line)
    start_rows = len(o)
    o.insert(2, [Stem(2), TetraLoop(), Stem(2)])
    assert len(o) == start_rows + 1
    assert len(o[2]) == 3
    assert all(isinstance(m, Motif) for m in o[2])

    # 5) Type validation
    with pytest.raises(ValueError):
        o.insert((0, 0), "not-a-motif")  # type: ignore
    with pytest.raises(ValueError):
        o.insert("bad-index", Stem(2))  # type: ignore


def test_pop_variants_and_effects():
    # Build an origami with two rows
    o = Origami([[Stem(3), TetraLoop(), Stem(2)], [Stem(4)]])

    # 1) Pop a single motif at (row, col)
    row0_len = len(o[0])
    popped_motif = o.pop((0, 1))
    assert isinstance(popped_motif, Motif)
    assert len(o[0]) == row0_len - 1

    # 2) Pop a whole line by row index
    rows_before = len(o)
    popped_line = o.pop(0)
    assert isinstance(popped_line, list)
    assert all(isinstance(m, Motif) for m in popped_line)
    assert len(o) == rows_before - 1

    # 3) Validation
    with pytest.raises(ValueError):
        o.pop("bad-index")  # type: ignore
    with pytest.raises(ValueError):
        o.pop((0,))  # not a 2-tuple


def test_improve_folding_pathway():
    origami = simple_origami([11, 11, 11], kl_columns=2)
    bar, penalty = origami.folding_barriers()
    assert penalty == 80
    origami = origami.improve_folding_pathway()
    new_bar, new_penalty = origami.folding_barriers()
    assert new_penalty < penalty
    assert new_penalty == 24
