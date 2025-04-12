import pytest
from pyfurnace.design import Strand
from pyfurnace.design import RIGHT, LEFT

# Test Fixture
@pytest.fixture
def strand():
    return Strand('A\\U|G/C', directionality='53', start=(0, 0), direction=RIGHT)

# Test Cases
def test_strand_length(strand):
    assert len(strand) == 7

def test_strand_getitem(strand):
    assert strand[0] == 'A'
    assert strand[1] == '\\'
    assert strand[2] == 'U'
    assert strand[3] == '|'
    assert strand[4] == 'G'
    assert strand[5] == '/'
    assert strand[6] == 'C'

def test_strand_setitem(strand):
    strand[0] = 'C'
    assert strand[0] == 'C'
    assert str(strand) == 'C\\U|G/C'

def test_strand_sequence(strand):
    assert strand.sequence == 'AUGC'

def test_strand_addition(strand):
    other_strand = Strand('CGUA', directionality='53', start=(0, 0), direction=RIGHT)
    result = strand + other_strand
    assert str(result) == 'A\\U|G/CCGUA'

def test_strand_in_place_addition(strand):
    other_strand = Strand('CGUA', directionality='53', start=(0, 0), direction=RIGHT)
    strand += other_strand
    assert str(strand) == 'A\\U|G/CCGUA'

def test_strand_contains(strand):
    assert 'AU' in strand.sequence
    assert Strand('A\\U', directionality='53', start=(0, 0), direction=RIGHT) in strand
    assert 'GC' in strand.sequence

def test_strand_equality():
    strand1 = Strand('A\\U|G/C', directionality='53', start=(0, 0), direction=RIGHT)
    strand2 = Strand('A\\U|G/C', directionality='53', start=(0, 0), direction=RIGHT)
    assert strand1 == strand2

def test_strand_inequality():
    strand1 = Strand('A\\U|G/C', directionality='53', start=(0, 0), direction=RIGHT)
    strand2 = Strand('CGUA', directionality='53', start=(0, 0), direction=RIGHT)
    assert strand1 != strand2

def test_strand_invert(strand):
    strand.invert()
    assert str(strand) == 'C╯G│U╮A'

def test_strand_flip(strand):
    strand.flip(horizontally=True, flip_start=True)
    assert str(strand) == 'A╭U│G╰C'
    assert strand.start == (1, 0)

    strand.flip(horizontally=False, vertically=True, flip_start=True)
    assert str(strand) =='A╰U│G╭C'
    assert strand.start == (1, 4)

def test_strand_reverse(strand):
    strand.reverse()
    assert strand.directionality == '35'

def test_strand_reverse_invert(strand):
    strand.reverse().invert()
    assert strand.directionality == '53'
    assert strand.sequence == 'CGUA'

def test_strand_shift(strand):
    strand.shift((2, 3))
    assert strand.start == (2, 3)

def test_strand_insert(strand):
    strand.insert(1, 'T')
    assert str(strand) == 'AT\\U|G/C'

def test_strand_pop(strand):
    popped_value = strand.pop(1)
    assert popped_value == '\\'
    assert str(strand) == 'AU|G/C'

def test_strand_joining():
    strand1 = Strand('A\\U|G/', directionality='53', start=(0, 0), direction=RIGHT)
    strand2 = Strand('/C\\U', directionality='53', start=(0, 4), direction=LEFT)
    strand1.join(strand2)
    assert str(strand1) == 'A╮U│G╯╭C╰U'

def test_strand_directionality():
    strand = Strand('A\\U|G/C', directionality='35', start=(0, 0), direction=RIGHT)
    assert strand.directionality == '35'

    strand.directionality = '53'
    assert strand.directionality == '53'

def test_strand_copy(strand):
    copied_strand = strand.copy()
    assert copied_strand == strand
    assert id(copied_strand) != id(strand)

def test_strand_draw(strand):
    canvas = strand.draw()
    assert 'A' in canvas
    assert 'C' in canvas
