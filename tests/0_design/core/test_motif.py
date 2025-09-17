import pytest
from pyfurnace.design.motifs import Motif
from pyfurnace.design.core.strand import Strand
from pyfurnace.design.core.basepair import BasePair
from pyfurnace.design.core import (
    RIGHT,
    LEFT,
    MotifStructureError,
    pseudo_to_dot,
    dot_bracket_to_tree,
    dot_bracket_to_pair_map,
)
from pyfurnace.design.utils import simple_origami


# Test Fixtures
@pytest.fixture
def m_two_strands():
    a = Strand("AGGGAUAAUAAUAUAUAU/U/-\\A\\G", start=(0, 2))
    b = Strand("UGCGUAUUAUUAUAUAUA\\AA\\U/GG/U", start=(0, 4), directionality="35")
    return Motif([a, b])


@pytest.fixture
def simple_motif():
    return Motif(Strand("AUCG-/||/AAAAAAAAAAAAAAA\\|/CUC/|\\A-", start=(0, 5)))


# -------- helpers --------


def h_strand(chars="──", start=(0, 0)):
    """Horizontal strand helper (rightward)."""
    return Strand(chars, start=start, direction=(1, 0), directionality="53")


def v_strand(chars="││", start=(0, 0)):
    """Vertical strand helper (downward)."""
    return Strand(chars, start=start, direction=(0, 1), directionality="53")


# =========================
# __contains__
# =========================


def test_motif_contains_strand_and_string_and_basepair():
    s = h_strand("──")
    m = Motif(s)

    # Strand membership
    assert s in m

    # String membership (substring within a Strand)
    assert "──" in m
    assert "xx" not in m

    # BasePair membership (set explicit basepair -> autopairing off)
    bp = BasePair({(0, 0): (2, 0)})
    m.basepair = bp
    assert (0, 0) in m.basepair
    assert bp in m


# =========================
# __eq__
# =========================


def test_motif_eq_same_strands_and_basepair():
    s1 = h_strand("──")
    s2 = h_strand("──")
    m1 = Motif(s1)
    m2 = Motif(s2)
    assert m1 == m2  # same strand geometry/content

    # Change one motif's strands -> not equal
    m3 = Motif(h_strand("─"))
    assert m1 != m3

    # Different basepair settings -> not equal
    m4 = Motif(h_strand("──"))
    m4.basepair = BasePair({(0, 0): (2, 0)})
    assert m1 != m4


# =========================
# autopairing setter
# =========================


def test_autopairing_setter_recomputes_basepair_no_error():
    # Two horizontal strands on y=0 and y=2 so vertical midpoint is free.
    top = h_strand("──", start=(0, 0))
    bot = h_strand("──", start=(0, 2))
    m = Motif(top, bot)

    # Give sequences (A/U pairs possible), but autopairing drives computation.
    top.sequence = "AU"
    bot.sequence = "UA"

    # Flip autopairing on — should (re)set basepair (may be empty but valid)
    m.autopairing = True
    assert m.autopairing is True
    assert isinstance(m.basepair, BasePair)  # existence/type check
    # Structure computed without raising
    _ = m.structure


# =========================
# sequence setter
# =========================


def test_sequence_setter_with_string_and_list_and_validation():
    s1 = h_strand("──")  # length 2
    s2 = h_strand("──")  # length 2
    m = Motif(s1, s2)

    # Set via single string split by '&'
    m.sequence = "AU&GC"
    assert s1.sequence == "AU"
    assert s2.sequence == "GC"
    assert str(m.sequence) == "AU&GC"

    # Set via list
    m.sequence = ["GG", "CC"]
    assert str(m.sequence) == "GG&CC"

    # Wrong number of segments -> ValueError
    with pytest.raises(ValueError):
        m.sequence = "AAA"  # only one segment


# =========================
# from_file
# =========================


def test_from_file_parses_simple_sketch(tmp_path):
    # Minimal one-strand horizontal motif
    content = "5──3\n"
    file_path = tmp_path / "sketch.txt"
    file_path.write_text(content, encoding="utf-8")

    m = Motif.from_file(str(file_path))
    assert isinstance(m, Motif)
    assert len(m) == 1
    # Drawable / non-empty blueprint
    assert str(m).strip() != ""


# =========================
# from_structure (non dot-bracket inputs)
# =========================


def test_from_structure_with_basepair_map():
    # Pair 0-1 (length 2)
    bp = BasePair({0: 1, 1: 0})
    m = Motif.from_structure(structure=bp, sequence="NN")
    assert isinstance(m, Motif)
    assert len(str(m.sequence).replace("&", "")) == 2
    # Should reflect a single paired unit in structure
    assert m.structure.replace("&", "").count("(") == 1
    assert m.structure.replace("&", "").count(")") == 1


# =========================
# _check_addition
# =========================


def test_check_addition_passes_when_junctions_exist():
    # Use strands without 5/3 terminals so junctions exist at both ends
    m1 = Motif(h_strand("──"))
    m2 = Motif(h_strand("──"))
    # Should NOT raise
    m1._check_addition(m2)  # right direction default
    with pytest.raises(MotifStructureError, match=r"in the direction Position\(0, 1\)"):
        m1._check_addition(m2, direction=(0, 1))  # can also ask for a direction


def test_check_addition_raises_when_junctions_missing():
    # '5' and '3' at ends -> no junctions
    m1 = Motif(Strand("5──3", start=(0, 0), direction=(1, 0)))
    m2 = Motif(Strand("5──3", start=(0, 0), direction=(1, 0)))

    with pytest.raises(MotifStructureError):
        m1._check_addition(m2)


# =========================
# extend_junctions (skip_axis)
# =========================


def test_extend_junctions_skip_axis_horizontal_only():
    # One short horizontal strand — has left/right junctions
    m = Motif(h_strand("──", start=(1, 1)))
    original_max = m.max_pos

    # Ask to extend only horizontally (skip vertical axis 0)
    m.extend_junctions(
        skip_axis=0, until=m.max_pos.replace(x=m.max_pos[0] + 4, y=m.max_pos[1])
    )
    # Horizontal extent should grow (x increases). y stays the same.
    assert m.max_pos[0] >= original_max[0] + 4
    assert m.max_pos[1] == original_max[1]


# =========================
# insert (index + join/copy combos + bad type)
# =========================


@pytest.mark.parametrize(
    "idx", [0, 1, 999]
)  # start, end, overshoot tolerated by list.insert
@pytest.mark.parametrize("join", [True, False])
@pytest.mark.parametrize("copy", [True, False])
def test_insert_variants(idx, join, copy):
    m = Motif(h_strand("──"))
    s_new = h_strand("──", start=(4, 0))

    m.insert(idx, s_new, join=join, copy=copy)

    # Motif grew
    assert len(m) >= 1

    # If copy=True, the inserted object identity should differ
    if copy:
        assert all(id(s) != id(s_new) for s in m)
    else:
        assert any(id(s) == id(s_new) for s in m)

    # Join may reduce strand count; we only check that motif remains valid
    assert isinstance(str(m), str) and str(m) != ""


def test_insert_rejects_non_strand():
    m = Motif(h_strand("──"))
    with pytest.raises(ValueError):
        m.insert(0, "not_a_strand")  # type: ignore


def test_motif_creation():
    """Test basic creation of Motif."""
    motif = Motif()
    assert motif is not None
    assert len(motif) == 0


def test_motif_creation_from_text(m_two_strands):

    # Using this weird string format because since the variable
    # 'text1' is indented, the string has a bad indentation
    # in the top and bottom lines. There are no problem with
    # triple quotes, if the variable is not indented.
    text1 = (
        "                   ╭─╮   \n"
        "                   U=A   \n"
        "5AGGGAUAAUAAUAUAUAU╯ ╰G  \n"
        " ┊ ┊ ┊┊┊┊┊┊┊┊┊┊┊┊┊┊   ┊  \n"
        " UGCGUAUUAUUAUAUAUA╮ ╭U5 \n"
        "                   A G   \n"
        "                   A G   \n"
        "                   ╰U╯   \n"
        "                         \n"
    )

    text2 = (
        "                   ╭─╮   \n"
        "                   U=A   \n"
        ">AGGGAUAAUAAUAUAUAU╯ ╰G  \n"
        " ┊ ┊ ┊┊┊┊┊┊┊┊┊┊┊┊┊┊   ┊  \n"
        " UGCGUAUUAUUAUAUAUA╮ ╭U< \n"
        "                   A G   \n"
        "                   A G   \n"
        "                   ╰U╯   \n"
        "                         \n"
    )

    motif1 = Motif.from_text(text1)
    motif2 = Motif.from_text(text2)
    assert motif1 == motif2
    motif1.strip()
    for s1, s2 in zip(motif1, m_two_strands):
        if s1.directionality == s2.directionality:
            assert s1 == s2
        else:
            assert s1 == s2.invert()


def test_m_two_strands(m_two_strands):
    """Test creating a motif with strands."""
    assert len(m_two_strands) == 2
    assert m_two_strands.sequence
    assert m_two_strands.structure is not None


def test_motif_addition(m_two_strands):
    """Test adding motifs."""
    motif1 = m_two_strands.copy()
    motif2 = m_two_strands.copy()
    combined = motif1 + motif2
    assert len(combined) == 2
    assert combined.structure == (
        "(.(.((((((((((((((()((.(.((((((((((((((()(&).."
        "...)))))))))))))).).)).....)))))))))))))).).)"
    )
    expect_seq = (
        str(m_two_strands[0].sequence) * 2
        + "&"
        + str(m_two_strands[1].sequence[::-1]) * 2
    )
    assert combined.sequence == expect_seq


def test_motif_in_place_addition(m_two_strands):
    """Test in-place addition of motifs."""
    motif1 = m_two_strands.copy()
    motif1 += m_two_strands.copy()
    assert len(motif1) == len(m_two_strands)
    expect_seq = (
        str(m_two_strands[0].sequence) * 2
        + "&"
        + str(m_two_strands[1].sequence[::-1]) * 2
    )
    assert motif1.sequence == expect_seq


def test_motif_rotation(m_two_strands):
    """Test rotating motifs."""
    motif = m_two_strands.copy()
    motif.rotate(times=1)
    assert motif.max_pos[0] == m_two_strands.max_pos[1]
    assert motif.max_pos[1] == m_two_strands.max_pos[0]


def test_motif_flip(m_two_strands):
    """Test flipping motifs horizontally and vertically."""
    motif = m_two_strands.copy()
    motif.flip(horizontally=True, vertically=False)
    assert motif[0].start == (21, 2)
    motif.flip(vertically=True, horizontally=False)
    assert motif[0].start == (21, 5)


def test_motif_shift(m_two_strands):
    """Test shifting motifs."""
    motif = m_two_strands.copy()
    motif.shift((3, 4))
    assert motif.min_pos == (3, 4)
    motif = m_two_strands.copy()
    motif.shift((3, 4), extend=True)
    assert motif.min_pos == (0, 4)


def test_motif_extend_junctions(m_two_strands):
    """Test extending junctions."""
    start_initial = m_two_strands[0].start
    motif = m_two_strands.copy()
    motif.shift((3, 4))
    motif.extend_junctions()
    assert motif[0].start == (0, start_initial[1] + 4)


def test_motif_basepair(m_two_strands):
    """Test basepair dictionary creation."""
    motif = m_two_strands.copy()
    basepair = motif.basepair
    assert isinstance(basepair, BasePair)
    assert len(basepair) >= 0


def test_motif_dot_bracket(m_two_strands):
    """Test dot-bracket notation."""
    motif = m_two_strands.copy()
    dot_bracket = motif.structure
    assert isinstance(dot_bracket, str)
    assert len(dot_bracket) > 0
    text = (
        "      ╭────A───╮                                       \n"
        " 5SSSA╯╭─NNNNNN╯╭─ASSS──────GGG──────CCC──────UUU───╮  \n"
        "  ┊┊┊  │ ┊┊┊┊┊┊ │  ┊┊┊      ┊┊┊      ┊┊┊      ┊┊┊   │  \n"
        "  SSSA─╯╭NNNNNN─╯╭ASSS─5   ╭CCC╮    ╭GGG╮  3──AAA╮  │  \n"
        "        ╰───A────╯         │   │    │   │        │  │  \n"
        "                           ╰───┼────┼───╯        │  │  \n"
        "                               │    ╰────────────╯  │  \n"
        "                               ╰────────────────────╯  \n"
    )
    motif = Motif.from_text(text)
    assert motif.structure == "(((..[[[[[[.)))&(((..]]]]]].)))((([[[{{{)))]]]}}}"
    assert motif.sequence == "SSSAANNNNNNASSS&SSSAANNNNNNASSSGGGCCCUUUCCCGGGAAA"


def test_motif_append_and_pop(m_two_strands):
    """Test appending and popping strands in a motif."""
    strand = Strand("A\\U|G/C", directionality="53", start=(0, 0), direction=RIGHT)
    motif = m_two_strands.copy()
    initial_len = len(motif)
    initial_structure = motif.structure
    motif.append(strand, join=False)
    assert len(motif) == initial_len + 1
    popped = motif.pop()
    assert popped == strand
    assert len(motif) == initial_len
    assert motif.structure == initial_structure


def test_motif_sequence_assignment(m_two_strands):
    """Test assigning sequences to a motif."""
    motif = m_two_strands.copy()
    sequences = ["AUCG", "GCGA"]
    motif.sequence = sequences
    assert motif.sequence == "&".join(sequences)


def test_motif_concat(m_two_strands):
    """Test concatenating motifs."""
    motif1 = m_two_strands.copy()
    motif2 = m_two_strands.copy()
    concatenated = Motif.concat(motifs=[motif1, motif2], axis=1, copy=True)
    combined = motif1 + motif2
    assert combined == concatenated


def test_motif_save_oxdna(m_two_strands, tmp_path):
    """Test saving motif in oxDNA format."""
    motif = m_two_strands.copy()
    filepath = tmp_path / "motif"
    motif.save_3d_model(
        str(filepath), forces="True", pk_forces=True, sequence=motif.sequence, pdb=True
    )
    assert filepath.with_suffix(".dat").exists()
    assert filepath.with_suffix(".top").exists()


def test_motif_save_text(m_two_strands, tmp_path):
    """Test saving motif in text format."""
    motif = m_two_strands.copy()
    filepath = tmp_path / "motif"
    motif.save_text(str(filepath))
    motif.save_fasta(str(filepath))
    assert filepath.with_suffix(".txt").exists()
    assert filepath.with_suffix(".fasta").exists()


def test_motif_copy(m_two_strands):
    """Test copying a motif."""
    motif = m_two_strands.copy()
    copy = motif.copy()
    assert motif == copy
    assert id(motif) != id(copy)


def test_motif_structure_manipulation(m_two_strands):
    """Test manually setting and getting structure."""
    motif = m_two_strands.copy()
    # set all bases to unpaired
    structure = "".join(["." if sym != "&" else sym for sym in m_two_strands.sequence])
    motif.structure = structure  # Reassigning should not cause an error
    assert motif.structure == structure  # Structure should be set correctly
    assert motif.basepair == BasePair()  # No basepairs should be present


def test_motif_strip(simple_motif):
    """Test stripping motifs."""
    motif = simple_motif.copy()
    stripped = motif.strip()
    assert stripped.min_pos == (0, 0)


def test_motif_bad_shift(simple_motif):
    """Test shifting motifs to negative positions."""
    motif = simple_motif.copy()
    with pytest.raises(MotifStructureError, match="The motif cannot be shifted."):
        motif.shift((-10, -10))  # Invalid shift tuple length


def test_strand_sorting():
    """Test sorting strands within a motif."""
    strand1 = Strand("AUGC", directionality="53", start=(0, 0), direction=RIGHT)
    strand2 = Strand("GAUG", directionality="35", start=(6, 2), direction=LEFT)
    strand3 = Strand("C\\|AU", directionality="53", start=(8, 3), direction=RIGHT)
    motif = Motif([strand1, strand2, strand3])
    motif.sort()
    assert motif[0] == strand1
    assert motif[1] == strand2
    assert motif[2] == strand3
    motif.sort(key=lambda s: 2 * (s[0] == "C") + 1 * (s[0] == "A"), reverse=True)
    assert motif[0] == strand3
    assert motif[1] == strand1
    assert motif[2] == strand2


def test_setitem_getitem(simple_motif):
    """Test __getitem__ and __setitem__ for Motif."""
    motif = simple_motif.copy()
    assert motif[0] == simple_motif[0], "__getitem__ did not return the correct strand."
    new_strand = Strand("AUGC", directionality="53", start=(0, 0), direction=RIGHT)
    motif[0] = new_strand
    assert motif[0] == new_strand, "__setitem__ did not set the strand correctly."
    with pytest.raises(ValueError, match="is not a Strand object."):
        motif[0] = "not_a_strand"  # Invalid type


def test_motif_from_structure():
    """Test creating a motif from structure."""
    structure = "....(((...)))..((..)).&..(((...)))"
    sequence = "AUAUGCCUUUGGCAACCUUGGA&AAGCGUAUCGC"
    motif = Motif.from_structure(sequence=sequence, structure=structure)
    assert len(motif) == 2
    assert motif.structure == structure
    assert motif.sequence == sequence
    motif = Motif.from_structure(sequence="GGGGUUCGCCCC")
    assert len(motif) == 1
    assert motif.structure == "((((....))))"
    assert motif.sequence == "GGGGUUCGCCCC"

    # convert an origami to motif
    origami = simple_origami([-2], align="center")
    assert len(origami.sequence) > 0
    assert len(origami.structure) > 0
    mot_from_str = Motif.from_structure(
        structure=origami.structure, sequence=origami.sequence
    )
    assert mot_from_str.sequence == origami.sequence
    assert mot_from_str.structure == origami.structure.translate(pseudo_to_dot)
    pk_from_str = mot_from_str[0].pk_info
    pk_values = list(origami.pseudoknots.values())
    for i, ind_fwd in enumerate(pk_from_str["ind_fwd"]):
        pk_id = pk_from_str["id"][i]
        if pk_id.endswith("'"):
            compl_id = pk_id[:-1]
        else:
            compl_id = pk_id + "'"
        ind_rev = pk_from_str["ind_fwd"][pk_from_str["id"].index(compl_id)]
        hit = False
        for pk_dict2 in pk_values:
            if ind_fwd == pk_dict2["ind_fwd"][0]:
                assert ind_rev == pk_dict2["ind_rev"][0]
                hit = True
                break
            elif ind_rev == pk_dict2["ind_fwd"][0]:
                assert ind_fwd == pk_dict2["ind_rev"][0]
                hit = True
                break
        assert hit, "Pseudoknot not found"

    # Build a simple "()" as a Node tree and pass it in.
    pm = dot_bracket_to_pair_map("...(...)....(((...)))")
    m = Motif.from_structure(structure=pm, sequence="NNNNNNNNNNNNNNUUUNNNN")
    assert isinstance(m, Motif)
    assert m.sequence == "NNNNNNNNNNNNNNUUUNNNN"
    assert m.structure == "...(...)....(((...)))"
    m = Motif.from_structure(structure=pm)
    assert isinstance(m, Motif)
    assert m.sequence == "N" * 21
    assert m.structure == "...(...)....(((...)))"

    # Build a simple "()" as a Node tree and pass it in.
    node = dot_bracket_to_tree(
        "...(...)....(((...)))", sequence="NNNNNNNNNNNNNNUUUNNNN"
    )
    m = Motif.from_structure(structure=node)
    assert isinstance(m, Motif)
    assert m.sequence == "NNNNNNNNNNNNNNUUUNNNN"
    assert m.structure == "...(...)....(((...)))"
    # Build a simple "()" as a Node tree and pass it in.
    node = dot_bracket_to_tree("...(...)....(((...)))")
    m = Motif.from_structure(structure=node)
    assert isinstance(m, Motif)
    assert m.sequence == "N" * 21
    assert m.structure == "...(...)....(((...)))"
