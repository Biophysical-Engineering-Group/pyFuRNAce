import pytest

from pyfurnace.design.core import Position, RIGHT, LEFT, Strand
from pyfurnace.design.utils.motif_lib import (
    Utils,
    start_end_stem,
    vertical_link,
    vertical_double_link,
    stem_cap_link,
)

# --------------------------------------------------------------------------- #
# Utils
# --------------------------------------------------------------------------- #


def test_utils_init_without_transforms_creates_motif():
    # minimal Utils with one strand, no flips/rotation
    s = Strand("─")
    m = Utils(strands=[s])
    assert len(m.strands) == 1
    s0 = m.strands[0]
    assert s0.strand == "─"
    assert s0.start == Position.zero()
    assert s0.direction == RIGHT


def test_utils_init_calls_flip_and_rotate(monkeypatch):
    called = {}

    def fake_flip(self, *, horizontally=False, vertically=False):
        called["flip"] = (horizontally, vertically)

    def fake_rotate(self, angle):
        called["rotate"] = angle

    monkeypatch.setattr(Utils, "flip", fake_flip)
    monkeypatch.setattr(Utils, "rotate", fake_rotate)

    s = Strand("─")
    m = Utils(strands=[s], hflip=True, vflip=True, rotate=180)

    # still a motif with the same strand object
    assert len(m.strands) == 1
    assert called["flip"] == (True, True)
    assert called["rotate"] == 180


# --------------------------------------------------------------------------- #
# start_end_stem
# --------------------------------------------------------------------------- #


def test_start_end_stem_default_labels_and_geometry():
    # defaults: up_left='3', up_right='5', down_left='-', down_right='-'
    m = start_end_stem()
    assert len(m.strands) == 3
    uR, uL, down = m.strands

    # top-left: "-3" at origin, RIGHT
    assert "3" in uL
    assert uL.start == Position.zero()
    assert uL.direction == RIGHT

    # top-right: "5-" at (3, 0), RIGHT
    assert "5" in uR
    assert uR.start == (3, 0)
    assert uR.direction == RIGHT

    # down strand
    assert len(set(down.strand)) == 1
    assert down.start == Position((4, 2))
    assert down.direction == LEFT


def test_start_end_stem_top_rule_appends_bar_when_both_top_are_dashes():
    m = start_end_stem(up_left="─", up_right="-")
    # both top entries are in "─-", so code appends "─" to up_left only
    # up_left strand is "-" + "──" = "-──"
    assert "─────" in m[0]
    assert "5" not in m


def test_start_end_stem_none_and_empty_suppress_strands():
    # None/"" are normalized to "", then skipped
    m = start_end_stem(up_left=None, up_right="", down_left="─", down_right=None)
    # we expect only one line strands at the top
    assert any(s.strand == "─-" for s in m.strands)  # down-left: "─-"
    # ensure we didn't create top strands
    assert all(s.strand not in ("-3", "5-") for s in m.strands)


def test_start_end_stem_rejects_invalid_values():
    with pytest.raises(ValueError):
        start_end_stem(up_left="X")  # invalid token


def test_start_end_stem_respects_explicit_strands_kwarg():
    custom = [Strand("X"), Strand("Y", start=(1, 1))]
    m = start_end_stem(strands=custom)
    # should use the provided list without adding more
    assert (
        m.strands is not custom
    )  # may copy internally, but content must match order/values
    assert [s.strand for s in m.strands] == ["X", "Y"]
    assert m.strands[1].start == (1, 1)


# --------------------------------------------------------------------------- #
# vertical links
# --------------------------------------------------------------------------- #


def test_vertical_link_single_strand_geometry_and_transforms(monkeypatch):
    called = {}
    monkeypatch.setattr(
        Utils, "flip", lambda self, **k: called.setdefault("flip", True)
    )
    monkeypatch.setattr(
        Utils, "rotate", lambda self, angle: called.setdefault("rot", angle)
    )

    m = vertical_link(hflip=True, rotate=90)
    assert len(m.strands) == 1
    s = m.strands[0]
    assert s.strand == "│"
    assert s.start == Position.zero()
    assert s.direction == (0, -1)
    # transform hooks triggered
    assert called["flip"] is True
    assert called["rot"] == 90


def test_vertical_double_link_two_strands_geometry():
    m = vertical_double_link()
    assert len(m.strands) == 2
    s1, s2 = m.strands

    assert s1.strand == "│"
    assert s1.start == Position.zero()
    assert s1.direction == (0, -1)

    assert s2.strand == "│"
    assert s2.start == (1, 0)
    assert s2.direction == (0, 1)


# --------------------------------------------------------------------------- #
# stem_cap_link
# --------------------------------------------------------------------------- #


def test_stem_cap_link_two_strands_with_curves_and_verticals():
    m = stem_cap_link()
    assert len(m.strands) == 2
    s1, s2 = m.strands

    assert s1.strand == "││╭─"
    assert s1.start == (0, 2)
    assert s1.direction == (0, -1)

    assert s2.strand == "╭"
    assert s2.start == (1, 2)
    assert s2.direction == (0, -1)
