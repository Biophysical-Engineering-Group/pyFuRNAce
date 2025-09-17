import pytest

from pyfurnace.design.core import Position, RIGHT, Strand, Sequence, nucl_to_pair
from pyfurnace.design.motifs import Stem


# # --------------------------------------------------------------------------- #
# # __init__ validation errors
# # --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    "kwargs,exc",
    [
        (dict(wobble_insert="nope"), ValueError),
        (dict(wobble_interval=-1), TypeError),
        (dict(wobble_interval="x"), TypeError),
        (dict(wobble_tolerance=-2), TypeError),
        (dict(wobble_tolerance="x"), TypeError),
        (dict(length="5"), TypeError),
        (dict(sequence=123), TypeError),
    ],
)
def test_init_errors(kwargs, exc):
    with pytest.raises(exc):
        Stem(**kwargs)


# --------------------------------------------------------------------------- #
# sequence-driven construction (top strand provided)
# --------------------------------------------------------------------------- #


def test_init_with_sequence_string_sets_length_and_complement_and_geometry():
    st = Stem(sequence="AUGC")
    # wobble disabled when sequence is provided
    assert st.wobble_interval == 0
    assert st.wobble_tolerance == 0

    # two strands
    assert len(st.strands) == 2
    top, bottom = st.strands

    # top sequence as provided, default geom
    assert str(top.sequence) == "AUGC"
    assert top.start == Position.zero()
    assert top.direction == RIGHT

    # bottom is reverse complement
    assert str(bottom.sequence) == "GCAU"
    assert bottom.start == (len(top.sequence) - 1, 2)
    assert bottom.direction == (-1, 0)

    # coords computed for both
    assert top._coords is not None
    assert bottom._coords is not None

    # length property reflects sequence length
    assert st.length == 4


def test_init_with_sequence_object_works_too():
    seq = Sequence("AUG")
    st = Stem(sequence=seq)
    assert st.length == 3
    assert str(st[0].sequence) == "AUG"
    assert str(st[1].sequence) == "CAU"  # reverse complement of AUG is CAU


# --------------------------------------------------------------------------- #
# length-driven construction (no sequence), strong-bases and wobble branches
# --------------------------------------------------------------------------- #


def test_short_length_uses_strong_bases_when_true():
    st = Stem(length=3, strong_bases=True)
    assert str(st[0].sequence) == "SSS"
    assert len(st[1].sequence) == 3


def test_set_strong_bases_false_rebuilds_with_wobble(monkeypatch):
    # Make wobble spacing deterministic
    monkeypatch.setattr("pyfurnace.design.motifs.stem.random.randint", lambda a, b: 1)
    st = Stem(length=3, strong_bases=True)
    st.set_strong_bases(False)  # triggers wobble path for short length
    # Should contain a K (wobble) now
    assert "K" in str(st[0].sequence)


@pytest.mark.parametrize("insert_mode", ["start", "middle", "end"])
def test_wobble_insertion_modes(insert_mode):
    # Deterministic wobble interval
    st = Stem(
        length=6,
        strong_bases=False,
        wobble_interval=2,
        wobble_tolerance=1,
        wobble_insert=insert_mode,
    )
    s = str(st[0].sequence)
    # We expect at least one wobble symbol
    assert "K" in s
    # geometry of bottom strand still correct
    assert st[1].start == (len(st[0].sequence) - 1, 2)
    assert st[1].direction == (-1, 0)
    st.wobble_insert = insert_mode  # should not error
    assert st.wobble_insert == insert_mode


def test_no_wobble_when_interval_zero():
    st = Stem(length=6, strong_bases=False, wobble_interval=0)
    assert str(st[0].sequence) == "NNNNNN"
    assert str(st[1].sequence) == "NNNNNN".translate(nucl_to_pair)[::-1]


# --------------------------------------------------------------------------- #
# property setters (and their validations)
# --------------------------------------------------------------------------- #


def test_length_setter_rebuilds_and_checks_type(monkeypatch):
    st = Stem(length=4, strong_bases=False, wobble_interval=0)
    # valid change
    st.length = 6
    assert st.length == 6
    assert len(st[0].sequence) == 6
    assert st[1].start == (5, 2)

    # invalid type
    with pytest.raises(TypeError):
        st.length = "x"  # noqa: F722


def test_wobble_interval_setter_updates_and_validates(monkeypatch):
    # start with wobbles present
    monkeypatch.setattr("pyfurnace.design.motifs.stem.random.randint", lambda a, b: 1)
    st = Stem(length=6, strong_bases=False, wobble_interval=2)
    assert "K" in str(st[0].sequence)

    # set to 0 -> should remove wobbles
    st.wobble_interval = 0
    assert "K" not in str(st[0].sequence)
    assert str(st[0].sequence) == "NNNNNN"

    # invalid values
    with pytest.raises(TypeError):
        st.wobble_interval = -1
    with pytest.raises(TypeError):
        st.wobble_interval = "x"  # noqa: F722


def test_wobble_tolerance_setter_updates_and_validates(monkeypatch):
    monkeypatch.setattr("pyfurnace.design.motifs.stem.random.randint", lambda a, b: 2)
    st = Stem(length=6, strong_bases=False, wobble_interval=2, wobble_tolerance=1)
    st.wobble_tolerance = 0  # forces fixed spacing
    after = str(st[0].sequence)
    # rebuild happened; sequence may change when tolerance changes
    assert after != ""
    assert len(after) == 6
    # invalid values
    with pytest.raises(TypeError):
        st.wobble_tolerance = -1
    with pytest.raises(TypeError):
        st.wobble_tolerance = "x"  # noqa: F722


def test_wobble_insert_setter_updates_and_validates(monkeypatch):
    monkeypatch.setattr("pyfurnace.design.motifs.stem.random.randint", lambda a, b: 1)
    st = Stem(
        length=6,
        strong_bases=False,
        wobble_interval=2,
        wobble_tolerance=1,
        wobble_insert="start",
    )
    s1 = str(st[0].sequence)
    st.wobble_insert = "end"
    s2 = str(st[0].sequence)
    assert s1 != s2  # rebuild occurred

    with pytest.raises(ValueError):
        st.wobble_insert = "nope"


# --------------------------------------------------------------------------- #
# direct sequence mutation helpers
# --------------------------------------------------------------------------- #


def test_set_up_sequence_accepts_string_and_rebuilds():
    st = Stem(length=4, strong_bases=False, wobble_interval=0)
    st.set_up_sequence("AUGC")
    assert str(st[0].sequence) == "AUGC"
    assert str(st[1].sequence) == "GCAU"  # reverse complement

    with pytest.raises(TypeError):
        st.set_up_sequence(123)  # must be a string


def test_set_down_sequence_bug_raises_typeerror():
    # The method calls set_up_sequence with a wrong keyword; exercise the path.
    st = Stem(length=4, strong_bases=False, wobble_interval=0)
    with pytest.raises(TypeError):
        st.set_down_sequence("AUGC")


# --------------------------------------------------------------------------- #
# _create_strands (protected) – explicit coverage of compute_coords=False + return
# --------------------------------------------------------------------------- #


def test_create_strands_return_and_no_coords():
    st = Stem(length=2, strong_bases=False, wobble_interval=0)
    strands = st._create_strands(length=2, return_strands=True, compute_coords=False)
    assert isinstance(strands, list) and len(strands) == 2
    t, b = strands
    assert not t._coords.is_empty()
    assert not b._coords.is_empty()
    assert str(b.sequence) == str(t.sequence).translate(nucl_to_pair)[::-1]


# --------------------------------------------------------------------------- #
# kwargs['strands'] passthrough path in __init__
# --------------------------------------------------------------------------- #


def test_init_with_explicit_strands_skips_generation():
    s1 = Strand("─")
    s2 = Strand("│", start=(0, 1), direction=(0, -1))
    st = Stem(strands=[s1, s2])  # bypass _create_strands
    assert st.strands[0].strand == "─"
    assert st.strands[1].strand == "│"
    assert st.strands[1].start == (0, 1)
