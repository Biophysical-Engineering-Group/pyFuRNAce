import pytest
from pyfurnace.design.core import Origami
from pyfurnace.design.core import Motif
from pyfurnace.design.motifs import Stem, TetraLoop

@pytest.fixture
def sample_origami():
    """Fixture to create a sample Origami object with a basic structure."""
    motif1 = TetraLoop()
    motif2 = Stem(2)
    motif3 = Stem(3)
    motif4 = TetraLoop(open_left=True)
    return Origami(matrix=[[motif1, motif2], [motif3, motif4]])

def test_getitem_single_motif(sample_origami):
    """Test __getitem__ for retrieving a single motif."""
    motif = sample_origami[0, 1]
    assert isinstance(motif, Motif), "Expected a Motif object."
    assert type(motif) == Stem, "Expected a Stem object."
    assert sample_origami[0, 1] == sample_origami[0][1], "Expected equivalent results for matrix and row/column indexing."
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
    assert all(isinstance(row, list) for row in submatrix), "Submatrix rows are not lists."
    assert all(isinstance(m, Motif) for row in submatrix for m in row), "Submatrix contains non-Motif elements."
    assert [[TetraLoop()], [Stem(3)]] == submatrix, "Submatrix content mismatch."

def test_getitem_function_mask(sample_origami):
    """Test __getitem__ with a function to screen motifs."""
    submatrix = sample_origami[lambda m: type(m) == Stem]
    assert isinstance(submatrix, list), "Expected a list of lists."
    assert len(submatrix) == 2, "Submatrix row count mismatch."
    assert all(isinstance(row, list) for row in submatrix), "Submatrix rows are not lists."
    assert all(isinstance(m, Motif) for row in submatrix for m in row), "Submatrix contains non-Motif elements."
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
    assert sample_origami[0] == [TetraLoop(), Stem(1), Stem(2), Stem(3)], "Row was not correctly set."
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
    sample_origami[lambda m: type(m) == TetraLoop] = new_submatrix
    assert len(sample_origami[0]) == 2, "Submatrix row count mismatch after setting."
    assert sample_origami[0][0] == Stem(1), "Submatrix element mismatch."
    assert sample_origami[1][1] == Stem(1), "Submatrix element mismatch."

    # Put back the original tetraloop
    new_submatrix = [[TetraLoop()], [TetraLoop(open_left=True)]]
    sample_origami[lambda m: type(m) == Stem and m.length == 1] = new_submatrix
    assert len(sample_origami[0]) == 2, "Submatrix row count mismatch after setting."
    assert sample_origami[0][0] == TetraLoop(), "Submatrix element mismatch."
    assert sample_origami[1][1] == TetraLoop(open_left=True), "Submatrix element mismatch."

def test_invalid_getitem(sample_origami):
    """Test __getitem__ with invalid input."""
    with pytest.raises(TypeError):
        _ = sample_origami["invalid_key"]

def test_invalid_setitem(sample_origami):
    """Test __setitem__ with invalid input."""
    with pytest.raises(ValueError):
        sample_origami[0, 0] = "not_a_motif"

def test_index(sample_origami):
    """Test the index method."""
    assert sample_origami.index(Stem(2)) == [(0, 1)], "Index not found."
    assert sample_origami.index(lambda m: type(m) == TetraLoop) == [(0, 0), (1, 1)], "Index not found."
