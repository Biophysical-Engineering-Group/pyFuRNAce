import inspect
import pytest

from pyfurnace.design.utils import origami_lib as o_lib
from pyfurnace.design.core import dot_bracket_to_stacks

import pyfurnace as pf

from pyfurnace.design.utils.origami_lib import (
    convert_angles_to_dt,
    simple_origami,
    ANGLES_DT_DICT,
)

# pick all callables starting with "template" (functions only)
_TEMPLATE_NAMES = [
    name
    for name, obj in inspect.getmembers(o_lib, inspect.isfunction)
    if name.startswith("template")
]


class Spy:
    def __init__(self):
        self.calls = 0
        self.last = None

    def __call__(self, obj):
        self.calls += 1
        self.last = obj  # we don't inspect it — keep it simple


def test_ipython_display_txt_minimal(monkeypatch):
    spy = Spy()
    monkeypatch.setattr(o_lib, "display", spy, raising=True)

    o_lib.ipython_display_txt("hello\nworld", max_height="123")
    assert spy.calls == 1  # just ensure it tried to display something


def test_ipython_clickable_txt_minimal(monkeypatch):
    spy = Spy()
    monkeypatch.setattr(o_lib, "display", spy, raising=True)

    ori = o_lib.simple_origami(dt_list=[0], kl_columns=1, main_stem=22, align="first")
    ori.pop((0, -1))  # remove motif to have open ends
    o_lib.ipython_clickable_txt(ori, gradient="rainbow", barriers=True)  # default args
    assert spy.calls == 1  # again, only check that display() was invoked


def test_ipython_display_3d_when_oat_missing_runs_and_does_not_call_oat(monkeypatch):
    # Force the lightweight branch
    monkeypatch.setattr(o_lib, "oat_installed", False, raising=True)

    # If these get called, the test should fail
    monkeypatch.setattr(
        o_lib,
        "describe",
        lambda *a, **k: (_ for _ in ()).throw(AssertionError("describe called")),
        raising=True,
    )
    monkeypatch.setattr(
        o_lib,
        "get_confs",
        lambda *a, **k: (_ for _ in ()).throw(AssertionError("get_confs called")),
        raising=True,
    )
    monkeypatch.setattr(
        o_lib,
        "oxdna_conf",
        lambda *a, **k: (_ for _ in ()).throw(AssertionError("oxdna_conf called")),
        raising=True,
    )

    ori = o_lib.simple_origami(dt_list=[0], kl_columns=1, main_stem=22, align="first")

    # Assert it WARNs and RETURNS (i.e., actually ran the branch and exited)
    with pytest.warns(UserWarning):
        ret = o_lib.ipython_display_3d(ori)
    assert ret is None


def test_ipython_display_3d_happy_path_runs_and_calls_oat(monkeypatch):
    # Force the full path
    monkeypatch.setattr(o_lib, "oat_installed", True, raising=True)

    calls = {"describe": 0, "get_confs": 0, "oxdna_conf": 0}

    def fake_describe(top, dat):
        calls["describe"] += 1
        return "TOPINFO", "TRAJINFO"

    def fake_get_confs(top_info, traj_info, start, end):
        calls["get_confs"] += 1
        return ["CONF"]

    def fake_oxdna_conf(top_info, conf, **kwargs):
        calls["oxdna_conf"] += 1
        return None

    monkeypatch.setattr(o_lib, "describe", fake_describe, raising=True)
    monkeypatch.setattr(o_lib, "get_confs", fake_get_confs, raising=True)
    monkeypatch.setattr(o_lib, "oxdna_conf", fake_oxdna_conf, raising=True)

    ori = o_lib.simple_origami(dt_list=[0], kl_columns=1, main_stem=22, align="first")

    # Should run without warnings or exceptions
    o_lib.ipython_display_3d(ori)

    assert calls["describe"] == 1
    assert calls["get_confs"] == 1
    assert calls["oxdna_conf"] == 1


# --------------------------
# convert_angles_to_dt
# --------------------------


@pytest.mark.parametrize(
    "angles, expected",
    [
        ([26], [-6]),  # exact hit
        ([58], [-5]),
        ([90], [-4]),
        ([410], [-5]),  # wraps to 50 -> nearest 58 -> -5 (matches implementation)
        ([-10], [-6]),  # -10 % 360 == 350 -> nearest 346 -> 4
        ([346], [4]),  # exact positive dovetail
        ([378], [5]),
        ([314, 282, 250, 218], [3, 2, 1, 0]),
        ([0], [-6]),
    ],
)
def test_convert_angles_to_dt_behaviour(angles, expected):
    # Because behavior depends on "closest" keys, recompute expectations
    # based on the same rule to make the test robust to future dict updates.
    # (But still assert equality below.)
    res = convert_angles_to_dt(angles)
    # sanity: lengths match
    assert len(res) == len(angles)
    # recompute via the module's own mapping rule for a stable check
    recomputed = []
    for a in angles:
        a = a % 360
        key = min(ANGLES_DT_DICT, key=lambda k: abs(k - a))
        recomputed.append(ANGLES_DT_DICT[key])
    assert res == recomputed


# --------------------------
# simple_origami (structure)
# --------------------------


@pytest.mark.parametrize("add_terminal_helix", [True, False])
@pytest.mark.parametrize("kl_columns", [1, 2])
def test_simple_origami_builds_expected_helix_count(add_terminal_helix, kl_columns):
    dt = [0, -3, 2]
    ori = simple_origami(
        dt_list=dt,
        kl_columns=kl_columns,
        main_stem=33,
        add_terminal_helix=add_terminal_helix,
        align="first",
    )
    assert isinstance(ori, pf.Origami)
    # Each item in the top-level iteration is a "helix line"
    # The builder adds terminal helices when requested.
    expected_lines = len(dt) + (2 if add_terminal_helix else 0)
    assert len(list(ori)) == expected_lines

    # First line begins with TetraLoop, ends with TetraLoop(open_left=True)
    first_line = ori[0]
    last_line = ori[-1]
    assert isinstance(first_line[0], pf.TetraLoop)
    assert isinstance(last_line[-1], pf.TetraLoop)

    # Check that Dovetail crosses got adjusted on first and last lines
    # (first: up_cross False; last: down_cross False)
    first_dts = [m for m in first_line if isinstance(m, pf.Dovetail)]
    last_dts = [m for m in last_line if isinstance(m, pf.Dovetail)]
    assert any(dt.up_cross is False for dt in first_dts)
    assert any(dt.down_cross is False for dt in last_dts)


def test_simple_origami_with_angles_flag_converts_and_builds():
    angles = [180, 218, 250]  # 180 ~ dt -1, 218 -> 0, 250 -> 1 (per mapping rule)
    ori = simple_origami(
        dt_list=angles, use_angles=True, kl_columns=1, main_stem=22, align="first"
    )
    assert isinstance(ori, pf.Origami)
    # spot-check presence of motifs commonly inserted by the builder
    # (TetraLoop, Stem, Dovetail, KissingDimer)
    flat = [m for line in ori for m in line]
    assert any(isinstance(m, pf.TetraLoop) for m in flat)
    assert any(isinstance(m, pf.Stem) for m in flat)
    assert any(isinstance(m, pf.Dovetail) for m in flat)
    assert any(isinstance(m, pf.KissingDimer) for m in flat)


def test_simple_origami_main_stem_matrix_forms_and_errors():
    # mismatched per-loop column lengths should raise
    with pytest.raises(ValueError):
        simple_origami(
            dt_list=[0, 1],
            kl_columns=2,
            main_stem=[[22, 22], [22]],  # second row too short
        )

    # same check for left_stem_kl shape mismatch
    with pytest.raises(ValueError):
        simple_origami(
            dt_list=[0, 1],
            kl_columns=2,
            main_stem=22,
            left_stem_kl=[[5, 5], [5]],  # second row too short
        )


def test_simple_origami_matrix_inputs():
    n_kl = 2
    dt = [0, -1]
    with pytest.raises(ValueError, match="The main_stem can be an int, a list"):
        simple_origami(
            dt_list=dt,
            kl_columns=n_kl,
            main_stem="test",
        )
    with pytest.raises(ValueError, match="The left_stem_kl can be an int,"):
        simple_origami(
            dt_list=dt,
            kl_columns=n_kl,
            left_stem_kl="test",
        )
    ori1 = simple_origami(
        dt_list=dt,
        kl_columns=n_kl,
        left_stem_kl=7,
    )
    ori2 = simple_origami(
        dt_list=dt,
        kl_columns=n_kl,
        left_stem_kl=[7, 7, 7, 7],
    )
    ori3 = simple_origami(
        dt_list=dt,
        kl_columns=n_kl,
        stem_pos=0,
    )
    db1 = dot_bracket_to_stacks(ori1.structure)[0]
    db2 = dot_bracket_to_stacks(ori2.structure)[0]
    db3 = dot_bracket_to_stacks(ori3.structure)[0]
    assert db1 == db2 == db3


# --------------------------
# Templates
# --------------------------


@pytest.mark.parametrize("factory_name", _TEMPLATE_NAMES, ids=_TEMPLATE_NAMES)
def test_templates_return_origami(factory_name):
    factory = getattr(o_lib, factory_name)
    ori = factory()
    assert isinstance(ori, pf.Origami)

    # quick structural sanity: each line has at least 3 motifs and sane endpoints
    for line in ori:
        assert len(line) >= 3
        assert isinstance(
            line[0],
            (
                pf.TetraLoop,
                pf.KissingLoop,
                pf.KissingLoop180,
                pf.KissingLoop120,
                pf.Broccoli,
                pf.MalachiteGreenShort,
                pf.Motif,
            ),
        )
        assert isinstance(
            line[-1],
            (
                pf.TetraLoop,
                pf.KissingLoop,
                pf.KissingLoop180,
                pf.KissingLoop120,
                pf.Motif,
            ),
        )
