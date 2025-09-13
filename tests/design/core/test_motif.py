import pytest
from pyfurnace.design.motifs import Motif
from pyfurnace.design.core.strand import Strand
from pyfurnace.design.core.basepair import BasePair
from pyfurnace.design.core import RIGHT, LEFT, MotifStructureError


# Test Fixtures
@pytest.fixture
def m_two_strands():
    a = Strand("AGGGAUAAUAAUAUAUAU/U/-\\A\\G", start=(0, 2))
    b = Strand("UGCGUAUUAUUAUAUAUA\\AA\\U/GG/U", start=(0, 4), directionality="35")
    return Motif([a, b])


@pytest.fixture
def simple_motif():
    return Motif(Strand("AUCG-/||/AAAAAAAAAAAAAAA\\|/CUC/|\\A-", start=(0, 5)))


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
