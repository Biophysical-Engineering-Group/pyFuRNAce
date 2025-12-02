import pytest
import numpy as np

from pyfurnace import ProteinCoords, Coords, KissingLoop180, L7Ae


###
### TEST CACHING MECHANISMS IN Coords CLASS
###


def test_helix_caching():
    pos = (0, 0, 0)
    b = (0, 0, 1)
    n = (1, 0, 0)
    # the first cached helix is the kissing loop helix
    stem1 = Coords.compute_helix_from_nucl(pos, b, n, 8, double=False)

    assert Coords._CACHED_AE_T is None
    assert Coords._CACHED_AE_T_INV is None
    assert len(Coords._CACHED_HELICES) == 1
    cached_helices = iter(Coords._CACHED_HELICES.values())
    assert len(next(cached_helices)) == len(
        stem1
    )  # stem1 cached helix, single stranded

    KissingLoop180()  # create the kissing loop to cache its helix
    assert len(Coords._CACHED_HELICES) == 3
    cached_helices = iter(Coords._CACHED_HELICES.values())
    next(cached_helices)  # skip stem1
    stem2 = next(cached_helices)
    stem3 = next(cached_helices)
    assert len(stem2) == 9 * 2  # kissing loop helix, double stranded
    # side coordinates created for the kissing loop
    assert len(stem3) == 6  # 6 bases kl, single stranded

    # change helix parameters resets the coords
    Coords.set_helix_params(inclination=-18.0)
    assert Coords.inclination == -18.0
    assert len(Coords._CACHED_HELICES) == 0


def test_compute_AE_crossover_caching():
    assert Coords._CACHED_AE_T is None
    assert Coords._CACHED_AE_T_INV is None
    AE_T, AE_T_INV = Coords.compute_AE_crossover()

    assert Coords._CACHED_AE_T is not None
    assert Coords._CACHED_AE_T_INV is not None
    assert AE_T is Coords._CACHED_AE_T
    assert AE_T_INV is Coords._CACHED_AE_T_INV

    # calling again uses the cached versions
    AE_T_2, AE_T_INV_2 = Coords.compute_AE_crossover()
    assert AE_T_2 is Coords._CACHED_AE_T
    assert AE_T_INV_2 is Coords._CACHED_AE_T_INV

    # changing AE parameters resets the cached transformations
    assert Coords.ae_v_shift == 1.006
    Coords.set_AE_crossover_params(ae_v_shift=-1.2)
    # create a stem to be cached
    Coords.compute_helix_from_nucl((0, 0, 0), (0, 0, 1), (1, 0, 0), 5, double=False)
    assert Coords.ae_v_shift == -1.2
    assert Coords._CACHED_AE_T is None
    assert Coords._CACHED_AE_T_INV is None
    # this doesn't affect cached helices
    assert len(Coords._CACHED_HELICES) > 0


###
### Test cases for ProteinCoords class
###


def test_protein_coords_initialization():
    # Test with no parameters
    protein = ProteinCoords()
    assert protein.sequence == ""
    assert protein.coords.shape == (0,)

    # Test with sequence and coordinates
    sequence = "ACDE"
    coords = np.array(
        [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0], [10.0, 11.0, 12.0]]
    )
    protein = ProteinCoords(sequence=sequence, coords=coords)
    assert protein.sequence == sequence
    assert np.array_equal(protein.coords, coords)


def test_protein_coords_invalid_coords_type():
    # Test invalid coords type (should be numpy array)
    sequence = "ACDE"
    coords = "invalid_coords_type"
    with pytest.raises(ValueError, match="The coords must be a numpy array"):
        ProteinCoords(sequence=sequence, coords=coords)


def test_protein_coords_invalid_sequence_type():
    # Test invalid sequence type (should be a string)
    sequence = 1234  # Invalid type
    coords = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    with pytest.raises(ValueError, match="The sequence must be a string"):
        ProteinCoords(sequence=sequence, coords=coords)


def test_protein_coords_copy():
    # Test the copy method
    sequence = "ACDE"
    coords = np.array(
        [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0], [10.0, 11.0, 12.0]]
    )
    protein = ProteinCoords(sequence=sequence, coords=coords)
    protein_copy = protein.copy()
    assert protein_copy.sequence == protein.sequence
    assert np.array_equal(protein_copy.coords, protein.coords)


def test_protein_coords_setter_getter():
    # Test getter and setter of sequence and coords
    protein = ProteinCoords()
    protein.sequence = "ACDE"
    protein.coords = np.array(
        [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0], [10.0, 11.0, 12.0]]
    )

    assert protein.sequence == "ACDE"
    assert np.array_equal(
        protein.coords,
        np.array(
            [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0], [10.0, 11.0, 12.0]]
        ),
    )


def test_protein_coords_getitem_setitem():
    # Test __getitem__ and __setitem__
    sequence = "ACDE"
    coords = np.array(
        [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0], [10.0, 11.0, 12.0]]
    )
    protein = ProteinCoords(sequence=sequence, coords=coords)

    assert np.array_equal(
        protein[0], np.array([1.0, 2.0, 3.0])
    )  # Get the first coordinate
    protein[0] = np.array([13.0, 14.0, 15.0])  # Set the first coordinate
    assert np.array_equal(protein[0], np.array([13.0, 14.0, 15.0]))  # Verify the update


def test_protein_coords_len():
    # Test __len__ method
    sequence = "ACDE"
    coords = np.array(
        [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0], [10.0, 11.0, 12.0]]
    )
    protein = ProteinCoords(sequence=sequence, coords=coords)
    assert len(protein) == 4  # The length of the sequence


def test_transform_method():
    # Test transform method (simple identity matrix example)
    protein = L7Ae()[1].coords.proteins[0]

    T_matrix = np.eye(4)  # Identity matrix
    protein.transform(T_matrix)

    # We expect the coordinates to remain unchanged with identity transformation
    assert np.array_equal(protein.coords, L7Ae()[1].coords.proteins[0].coords)


###
### Test cases for Coords class
###


def test_coords_initialization():
    # Test with no parameters
    coords = Coords()
    assert coords.array.shape == (0,)

    # Test with some data
    coords_data = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    coords = Coords(input_array=coords_data)
    assert np.array_equal(coords.array, coords_data)


def test_coords_invalid_array_type():
    # Test invalid array type for coords
    with pytest.raises(
        ValueError,
        match="The Coordinates array must be a numpy array or a list of coordinates",
    ):
        Coords(input_array="invalid_array")


def test_coords_getitem_setitem():
    # Test __getitem__ and __setitem__ methods for Coords
    coords_data = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    coords = Coords(input_array=coords_data)

    assert np.array_equal(coords[0], np.array([1.0, 2.0, 3.0]))
    coords[0] = np.array([7.0, 8.0, 9.0])
    assert np.array_equal(coords[0], np.array([7.0, 8.0, 9.0]))


def test_coords_len():
    # Test __len__ method for Coords
    coords_data = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    coords = Coords(input_array=coords_data)
    assert len(coords) == 2  # Two coordinates


def test_coords_dummy_ends():
    # Test dummy ends setter and getter
    dummy_ends = (np.array([1.0, 1.0, 1.0]), np.array([2.0, 2.0, 2.0]))
    coords = Coords(dummy_ends=dummy_ends)
    assert np.array_equal(coords.dummy_ends[0], np.array([1.0, 1.0, 1.0]))
    assert np.array_equal(coords.dummy_ends[1], np.array([2.0, 2.0, 2.0]))


def test_coords_proteins():
    # Test proteins attribute
    protein1 = ProteinCoords(
        "ACDE", np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]])
    )
    coords = Coords(proteins=[protein1])
    assert len(coords.proteins) == 1
    assert coords.proteins[0] == protein1
