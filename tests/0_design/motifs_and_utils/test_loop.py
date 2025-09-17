import pathlib
import pytest

from pyfurnace.design.motifs import Loop, TetraLoop, CONFS_PATH, Coords, Strand


# -------- Loop --------


def test_loop_no_sequence_creates_no_strands_and_sets_up(cleanup_tmpdir=None):
    """Loop without sequence should not create any strands."""
    m = Loop()
    # The public contract (per docstring) is that motifs expose `.strands`
    assert hasattr(m, "strands")
    assert isinstance(m.strands, list)
    assert len(m.strands) == 0


@pytest.mark.parametrize("sequence", ["A", "AC", "ACG", "GGAU"])
def test_loop_with_sequence_builds_single_strand(sequence):
    """
    With a sequence, Loop should create exactly one Strand whose raw 'strand'
    text matches the constructor logic:
        '─' * len(sequence) + '╰│╭' + sequence
    and with start=(len(sequence), 2), direction=(-1, 0).
    """
    m = Loop(sequence=sequence)
    assert len(m.strands) == 1
    s = m.strands[0]
    # Check the literal structure string passed into Strand
    expected = "─" * len(sequence) + "╰│╭" + sequence
    assert isinstance(s, Strand)
    assert getattr(s, "strand") == expected
    assert getattr(s, "start") == (len(sequence), 2)
    assert getattr(s, "direction") == (-1, 0)


def test_loop_open_left_calls_flip_without_error():
    """
    When open_left=True, Loop flips horizontally and vertically.
    We don't rely on Motif internals — just ensure construction succeeds.
    """
    m = Loop(open_left=True, sequence="ACG")
    # Must still have the single created strand
    assert len(m.strands) == 1
    assert isinstance(m.strands[0], Strand)


# -------- TetraLoop --------


def test_tetraloop_rejects_wrong_length():
    with pytest.raises(ValueError):
        TetraLoop(sequence="AAA")  # not 4 nt
    with pytest.raises(ValueError):
        TetraLoop(sequence="AAAAA")  # not 4 nt


def test_tetraloop_default_constructs_strand_and_loads_coords():
    """
    Default TetraLoop(): sequence 'UUCG' => strand text 'UU╰│╭CG',
    start=(2,2), direction=(-1,0), and coords loaded from
    CONFS_PATH / 'TetraLoop.dat'.
    """
    assert isinstance(CONFS_PATH, pathlib.Path)
    t = TetraLoop()  # default "UUCG"

    assert len(t.strands) == 1
    s = t.strands[0]
    assert isinstance(s, Strand)

    # Strand layout constructed in __init__
    assert getattr(s, "strand") == "UU╰│╭CG"
    assert getattr(s, "start") == (2, 2)
    assert getattr(s, "direction") == (-1, 0)

    # The constructor sets s._coords = Coords.load_from_file(...)
    coords = getattr(s, "_coords", None)
    assert coords is not None, "Expected coordinates to be loaded on the strand"
    assert isinstance(
        coords, Coords
    ), "Expected a Coords instance attached to the strand"


def test_tetraloop_open_left_builds_and_flips_without_error():
    t = TetraLoop(open_left=True, sequence="UUCG")
    assert len(t.strands) == 1
    assert isinstance(t.strands[0], Strand)


def test_tetraloop_respects_custom_strands_and_does_not_override_them():
    """
    If 'strands' is supplied via kwargs, TetraLoop must use those and skip
    creating its own strand (and thus not touch _coords on the user-provided one).
    """
    strand = TetraLoop(sequence="GNRA").strands[0]

    t = TetraLoop(sequence="AACG", strands=[strand])
    assert len(t.strands) == 1
    assert t.strands[0] is strand


def test_tetraloop_set_sequence_happy_and_unhappy_paths():
    t = TetraLoop(sequence="UUCG")

    # valid
    t.set_sequence("AAAA")
    # The implementation sets t[0].sequence = new_sequence; with real classes,
    # this should be accessible as 'sequence' on the first strand.
    first = t.strands[0]
    assert getattr(first, "sequence") == "AAAA"

    # invalid
    with pytest.raises(ValueError):
        t.set_sequence("AAA")
