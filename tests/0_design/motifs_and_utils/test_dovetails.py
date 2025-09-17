import pytest
from pyfurnace.design.motifs import Dovetail


def test_dovetail_init_and_sign():
    with pytest.raises(ValueError, match="The length parameter must be an integer."):
        _ = Dovetail("test")

    positive = Dovetail(sequence="AUCG", sign=1)
    assert positive.length == 4
    assert "AUCG" in str(positive)
    assert "AUCG" in repr(positive)

    negative = Dovetail(sequence="AUCG", sign=-1)
    assert negative.length == -4
    assert "AUCG" in str(negative)
    assert "AUCG" in repr(negative)


def test_cross():
    dt = Dovetail(1, up_cross=True, down_cross=True)
    assert len(dt) == 4
    dt.up_cross = not dt.up_cross
    assert len(dt) == 3
    dt.down_cross = not dt.down_cross
    assert len(dt) == 2


def test_setup_sequence():
    dt = Dovetail(sequence="AUCG", sign=-1)
    with pytest.raises(TypeError, match="The sequence of a stem must be a"):
        dt.set_up_sequence(0)
    with pytest.raises(ValueError, match="The sign of the dovetail must be -1,"):
        dt.set_up_sequence("AC", sign="test")
    dt.set_up_sequence("AAUUCCGG")
    assert dt[1].sequence == "AAUUCCGG"
    assert dt.length < 0
    dt.length = 10
    assert dt.length == 10
    dt.set_up_sequence("AAUUCCGG")
    assert dt.length > 0

    dt.set_down_sequence("CCGGAAUU", sign=1)
    assert dt[0].sequence == "AAUUCCGG"
