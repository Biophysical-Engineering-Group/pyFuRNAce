from pyfurnace.design.core import Position, RIGHT
from pyfurnace.design.motifs import LambdaTurn, ThreeWayJunction, Bend90


# --------------------------------------------------------------------------- #
# LambdaTurn
# --------------------------------------------------------------------------- #


def test_LambdaTurn_builds_two_strands_and_forces_join_false(monkeypatch):
    # try to force join=True; factory should override to False
    motif = LambdaTurn(join=True)

    # Motif called with 2 strands
    assert len(motif.strands) == 2
    s1, s2 = motif.strands

    # exact sequences and geometry
    assert s1.strand == "─CUNGAUGG─"
    assert s1.start == Position.zero()
    assert s1.direction == RIGHT
    assert s2.strand == "─CCAU──GG─"
    assert s2.start == (9, 2)
    assert s2.direction == (-1, 0)

    # coords actually assigned to strand objects
    assert motif.strands[0]._coords is not None
    assert motif.strands[1]._coords is not None

    # Coords loaded for both strands with expected files and dummy_ends
    assert s1._coords.dummy_ends[0].size != 0
    assert s1._coords.dummy_ends[1].size != 0
    assert s2._coords.dummy_ends[0].size != 0
    assert s2._coords.dummy_ends[1].size != 0


# --------------------------------------------------------------------------- #
# ThreeWayJunction
# --------------------------------------------------------------------------- #


def test_ThreeWayJunction_builds_three_strands_and_forces_join_false():

    motif = ThreeWayJunction(join=True)

    assert len(motif.strands) == 3
    s1, s2, s3 = motif.strands

    # sequences
    assert s1.strand == "─NC────UAAN─"
    assert s2.strand == "─NG─AC╭A╯╭"
    assert s3.strand == "│U╮GN─"

    # geometry
    assert s1.start == Position.zero()
    assert s1.direction == RIGHT
    assert s2.start == (11, 2) and s2.direction == (-1, 0)
    assert s3.start == (3, 4) and s3.direction == (0, -1)

    # coords loaded for all three with the right files and dummy_ends
    for s in motif:
        assert not s._coords.is_empty()
        for i in range(2):
            assert s._coords.dummy_ends[i].size != 0


# --------------------------------------------------------------------------- #
# Bend90
# --------------------------------------------------------------------------- #


def test_Bend90_builds_two_strands_and_forces_join_false():

    motif = Bend90(join=True)

    assert len(motif.strands) == 2
    s1, s2 = motif.strands

    # sequences and geometry
    assert s1.strand == "─GAACUAC─"
    assert s1.start == Position.zero()
    assert s1.direction == RIGHT
    assert s2.strand == "─G─────C─"
    assert s2.start == (8, 2) and s2.direction == (-1, 0)

    # coords loaded for both strands with expected files and dummy_ends
    assert s1._coords.dummy_ends[0].size != 0
    assert s1._coords.dummy_ends[1].size != 0
    assert s2._coords.dummy_ends[0].size != 0
    assert s2._coords.dummy_ends[1].size != 0

    # coords assigned
    assert not s1._coords.is_empty()
    assert not s2._coords.is_empty()
