import pytest
import os
import numpy as np
from pyfurnace.design import (
    AmbiguosStructure,
    MotifStructureError,
    Position,
    Sequence,
    StrandsBlock,
    Coords,
    Strand,
)
from pyfurnace.design import RIGHT, LEFT, UP, DOWN


# Test Fixture
@pytest.fixture
def strand():
    strand = Strand("A\\U|G/C", directionality="53", start=(0, 0), direction=RIGHT)
    return strand.copy()


# Test Fixture 2
@pytest.fixture
def strand_all_dir():
    strand = Strand("--A/U|", directionality="53", start=(5, 5), direction=RIGHT)
    return strand.copy()


def test_line_input():
    with pytest.raises(
        TypeError, match="The strand must be a string or a sequence object. "
    ):
        Strand(1)
    with pytest.raises(
        MotifStructureError,
        match="The strand can have only one start and one end symbol",
    ):
        Strand("A\\U|G/33")
    with pytest.raises(
        MotifStructureError, match="The '5' and '3' symbols are terminal symbols."
    ):
        Strand("5A\\U|3G/")
    with pytest.raises(
        ValueError, match="The string 'AUC__' contains symbols not allowed in ROAD"
    ):
        Strand("auC__")


# Test Cases
def test_strand_length(strand):
    assert len(strand) == 7


def test_strand_getitem(strand):
    assert strand[0] == "A"
    assert strand[1] == "\\"
    assert strand[2] == "U"
    assert strand[3] == "|"
    assert strand[4] == "G"
    assert strand[5] == "/"
    assert strand[6] == "C"


def test_strand_setitem(strand):
    strand[0] = "C"
    assert strand[0] == "C"
    assert str(strand) == "C\\U|G/C"


def test_strand_sequence(strand):
    assert strand.sequence == "AUGC"


def test_equivalence(strand):
    assert strand == "A\\U|G/C"
    assert strand == Sequence("AUGC")
    assert strand != 1
    assert strand.__repr__() == "A\\U|G/C"


def test_strand_addition(strand):
    strand = 0 + strand
    assert strand.start == Position.zero()
    end = strand.end
    strand = "-" + strand
    assert strand.end == end + RIGHT
    other_strand = Strand("CGUA", directionality="53", start=(0, 0), direction=RIGHT)
    result = strand + other_strand
    assert str(result) == "─A╮U│G╯CCGUA"
    other_strand.reverse()
    with pytest.raises(
        MotifStructureError,
        match="Cannot add two strands with different directionality",
    ):
        strand + other_strand
    with pytest.raises(
        MotifStructureError,
        match="Cannot add a strand with a Sequence with different directionality",
    ):
        strand + Sequence("UUG", directionality="35")
    strand2 = strand + Sequence("UUG")
    assert str(strand2) == "─A╮U│G╯CUUG"
    with pytest.raises(TypeError, match="1 is not a valid type for addition"):
        strand + 1


def test_strand_in_place_addition(strand):
    other_strand = Strand("CGUA", directionality="53", start=(0, 0), direction=RIGHT)
    strand += other_strand
    assert str(strand) == "A\\U|G/CCGUA"


def test_strand_contains(strand):
    assert "AU" in strand.sequence
    assert Strand("A\\U", directionality="53", start=(0, 0), direction=RIGHT) in strand
    assert "GC" in strand.sequence


def test_strand_equality():
    strand1 = Strand("A\\U|G/C", directionality="53", start=(0, 0), direction=RIGHT)
    strand2 = Strand("A\\U|G/C", directionality="53", start=(0, 0), direction=RIGHT)
    assert strand1 == strand2


def test_strand_inequality():
    strand1 = Strand("A\\U|G/C", directionality="53", start=(0, 0), direction=RIGHT)
    strand2 = Strand("CGUA", directionality="53", start=(0, 0), direction=RIGHT)
    assert strand1 != strand2


def test_strand_invert(strand):
    strand.invert()
    assert str(strand) == "C╯G│U╮A"


def test_strand_flip(strand):
    strand.flip(horizontally=True, flip_start=True)
    assert str(strand) == "A╭U│G╰C"
    assert strand.start == (1, 0)

    strand.flip(horizontally=False, vertically=True, flip_start=True)
    assert str(strand) == "A╰U│G╭C"
    assert strand.start == (1, 4)


def test_strand_reverse(strand):
    strand.reverse()
    assert strand.directionality == "35"


def test_strand_reverse_invert(strand):
    strand.reverse().invert()
    assert strand.directionality == "53"
    assert strand.sequence == "CGUA"


def test_strand_shift(strand):
    strand.shift((2, 3))
    assert strand.start == (2, 3)


def test_strand_insert(strand):
    strand.insert(1, "T")
    assert str(strand) == "AT\\U|G/C"


def test_strand_pop(strand):
    popped_value = strand.pop(1)
    assert popped_value == "\\"
    assert str(strand) == "AU|G/C"


def test_strand_joining():
    strand1 = Strand("A\\U|G/", directionality="53", start=(0, 0), direction=RIGHT)
    strand2 = Strand("/C\\U", directionality="53", start=(0, 4), direction=LEFT)
    strand1.join(strand2)
    assert str(strand1) == "A╮U│G╯╭C╰U"
    with pytest.raises(ValueError, match="The object to join is not a Strand object"):
        strand1.join("invalid")
    with pytest.raises(
        TypeError, match="The objects to join are not a Strand object, "
    ):
        Strand.join_strands(strand1, "invalid")
    # test the four join method, for inplace and not
    strand1 = Strand("-\\|||/", directionality="53", start=(0, 0), direction=RIGHT)
    strand2 = Strand("/C\\U", directionality="53", start=(0, 4), direction=LEFT)
    for s1 in [strand1.copy(), strand2.copy()]:
        for s2 in [strand2.copy(), strand1.copy()]:
            for invert in [True, False]:
                if invert:
                    s2.invert()
                joined = Strand.join_strands(s1, s2)
                if joined:
                    assert str(joined) == "─╮│││╯╭C╰U"
                    s1.join(s2)
                    assert str(s1) == "─╮│││╯╭C╰U"


def test_strand_get_base(strand):
    base = strand.get_base_at_pos(Position((1, 1)))
    assert base == "U"
    char = strand.get_char_at_pos(Position((1, 0)))
    assert char == "╮"
    pos_map = strand.get_position_map()
    assert pos_map[Position((1, 0))] == "╮"


def test_strand_directionality():
    strand = Strand(
        "A\\U|G/C",
        directionality="35",
        start=(0, 0),
        direction=RIGHT,
        strands_block=StrandsBlock(),
    )
    assert strand.directionality == "35"

    strand.directionality = "53"
    assert strand.directionality == "53"


def test_strand_properties(strand):
    orig_directions = strand.directions
    strand.direction = DOWN
    strand.direction = RIGHT
    assert strand.directions == orig_directions
    assert strand.max_pos == RIGHT + DOWN * 4
    strand.start = RIGHT
    assert strand.max_pos == RIGHT * 2 + DOWN * 4
    assert strand.min_pos == RIGHT
    strand.start = Position.zero()
    assert strand.min_pos == Position.zero()
    strand.start = Position.zero()
    assert strand.next_pos == LEFT + DOWN * 4
    assert strand.sequence_list == [
        Sequence("A"),
        Sequence("U"),
        Sequence("G"),
        Sequence("C"),
    ]
    assert strand.sequence_slice == [slice(0, 1), slice(2, 3), slice(4, 5), slice(6, 7)]
    assert strand.seq_positions == (
        Position((0, 0)),
        Position((1, 1)),
        Position((1, 3)),
        Position((0, 4)),
    )
    strand.start = Position.zero()
    assert strand.seq_positions == (
        Position((0, 0)),
        Position((1, 1)),
        Position((1, 3)),
        Position((0, 4)),
    )
    with pytest.raises(
        ValueError, match="The strands block must be a StrandsBlock object. "
    ):
        strand.strands_block = "invalid"
    with pytest.raises(
        ValueError, match=r"The 2D coordinates must be a tuple/list of \(x,y\) "
    ):
        strand.start = 1
    with pytest.raises(
        ValueError, match=r"2D direction not allowed. The allowed values are:"
    ):
        strand.direction = Position((1, 1))


def test_bad_strands(strand):
    strand.start = LEFT
    with pytest.raises(
        MotifStructureError, match=r"The strand reaches negative coordinates: "
    ):
        strand._calculate_positions()
    strand.start = Position.zero()
    strand.direction = DOWN
    with pytest.raises(
        MotifStructureError,
        match=r"The x direction of the strand \(dx: 1\) "
        r"is not compatible with the next symbol \(│\)",
    ):
        strand._calculate_positions()
    strand.strand = "-" + strand.strand[1:]
    with pytest.raises(
        MotifStructureError,
        match=r"The y direction of the strand \(dy: 1\) "
        r"is not compatible with the next symbol \(─\)",
    ):
        strand._calculate_positions()
    strand = Strand("-\\/\\|")
    with pytest.warns(
        AmbiguosStructure, match=r"The strand is doing a crossing not allowed, "
    ):
        strand._calculate_positions()


def test_strand_copy(strand):
    copied_strand = strand.copy()
    assert copied_strand == strand
    assert id(copied_strand) != id(strand)


def test_strand_extend(strand_all_dir):
    strand_all_dir.extend(LEFT)
    assert strand_all_dir.start == (0, 5)
    assert len(strand_all_dir.positions) == 11
    strand_all_dir.extend(UP, until=(1, 1))
    assert strand_all_dir.end == (8, 1)
    assert len(strand_all_dir.positions) == 13


def test_strand_draw(strand):
    canvas = strand.draw()
    assert "A" in canvas
    assert "C" in canvas


def test_save_3d_model(tmp_path, strand):
    conf_text, top_text = strand.save_3d_model(return_text=True)
    assert conf_text.startswith("t = 0\nb = 100 100 100\nE = 0 0 0\n")
    seq = strand.sequence
    def_text = f"{len(seq)} 1 5->3\n" + str(seq) + " type=RNA circular=false \n"
    assert top_text == def_text

    strand.reverse()
    strand.coords = strand.coords.reverse()
    tmp_path = str(tmp_path / "strand")
    strand.save_3d_model(str(tmp_path), pdb=True)
    seq = strand.sequence[::-1]
    def_text = f"{len(seq)} 1 5->3\n" + str(seq) + " type=RNA circular=false \n"

    pdb_path = tmp_path + ".pdb"
    dat_path = tmp_path + ".dat"
    top_path = tmp_path + ".top"
    assert os.path.exists(dat_path)
    assert os.path.exists(top_path)
    assert os.path.exists(pdb_path)
    with open(top_path, "r") as f:
        top_content = f.read()
    with open(dat_path, "r") as f:
        dat_content = f.read()
    assert def_text == top_content
    assert conf_text == dat_content


def test_minimal_dimension():
    strand = Strand(
        "A\\U|G/C", directionality="53", start="minimal", strands_block=None
    )
    assert strand.start == Position.zero()


def test_contains(strand):
    assert Position.zero() in strand
    assert Sequence("AUG") in strand
    assert Strand("N") not in strand
    assert 2 not in strand


def test_coords(strand):
    with pytest.raises(MotifStructureError, match=r"The number of oxDNA coordinates "):
        strand.coords = strand.coords[:-1]
    strand.coords = strand.coords._array
    assert isinstance(strand._coords, Coords)
    strand.coords = np.array(())
    # this will empty the coords
    assert strand._coords.is_empty()
    # if you need the coords, it recreates them
    assert not strand.coords.is_empty()
