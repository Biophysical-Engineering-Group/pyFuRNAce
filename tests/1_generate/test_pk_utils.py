# tests/test_pseudoknots.py
import pytest

# Adjust this import if your module path differs:
from pyfurnace.generate.pk_utils import (
    parse_pseudoknots,
    add_untracked_pseudoknots,
)


# ----------------------------- parse_pseudoknots ----------------------------- #


def test_parse_pseudoknots_happy_path_and_defaults():
    # Includes commas inside list values to exercise the "attr_fixed" join logic.
    s = (
        "id:pk1,ind_fwd:[(10,20),(21,31)],ind_rev:[(40,50),(51,61)],E:-7.5,dE:0.5 ; "
        " id:pk2, ind_fwd:[(0,5)], ind_rev:[(95,100)] "  # uses default E/dE
    )

    with pytest.warns(UserWarning, match="Skipping pseudoknot with id pk1"):
        parse_pseudoknots("id:pk1,ind_fwd:[(20,10)],ind_rev:[(10, 20)],E:-7.5,dE:0.5 ;")

    out = parse_pseudoknots(s)
    assert set(out.keys()) == {"pk1", "pk2"}

    pk1 = out["pk1"]
    assert pk1["ind_fwd"] == [(10, 20), (21, 31)]
    assert pk1["ind_rev"] == [(40, 50), (51, 61)]
    assert pk1["E"] == -7.5
    assert pk1["dE"] == 0.5

    pk2 = out["pk2"]
    # defaults
    assert pk2["E"] == -9.0
    assert pk2["dE"] == 1.0
    assert pk2["ind_fwd"] == [(0, 5)]
    assert pk2["ind_rev"] == [(95, 100)]


def test_parse_pseudoknots_missing_id_raises():
    s = "ind_fwd:[(0,5)],ind_rev:[(10,15)]"
    with pytest.raises(ValueError, match="Pseudoknot id is missing"):
        _ = parse_pseudoknots(s)


def test_parse_pseudoknots_missing_indices_warns_and_skips():
    s = "id:x,E:-8.0"  # no ind_fwd / ind_rev
    with pytest.warns(
        UserWarning, match="Skipping pseudoknot with id x due to missing indices"
    ):
        out = parse_pseudoknots(s)
    assert "x" not in out


def test_parse_pseudoknots_inconsistent_lengths_warns_and_skips():
    # (0,5) length 5 vs (10,21) length 11 -> inconsistent -> warn + skip
    s = "id:bad,ind_fwd:[(0,5)],ind_rev:[(10,21)]"
    with pytest.warns(
        UserWarning, match=r"Skipping pseudoknot with id bad due to invalid indices"
    ):
        out = parse_pseudoknots(s)
    assert "bad" not in out


def test_parse_pseudoknots_literal_eval_error_prints_and_skips(capsys):
    # Make one of the lists unparsable; should print an error,
    # then skip for missing indices.
    s = "id:oops,ind_fwd:NOT_A_LIST,ind_rev:[(1,2)]"
    with pytest.warns(
        UserWarning, match="Skipping pseudoknot with id oops due to missing indices"
    ):
        out = parse_pseudoknots(s)
    assert "oops" not in out


def test_parse_pseudoknots_ignores_extra_whitespace_and_empty_entries():
    s = " ;  ; id:A, ind_fwd:[(1,6)], ind_rev:[(15,20)] ;  "
    out = parse_pseudoknots(s)
    print(out)
    assert set(out.keys()) == {"A"}
    assert out["A"]["ind_fwd"] == [(1, 6)]
    assert out["A"]["ind_rev"] == [(15, 20)]


# -------------------------- add_untracked_pseudoknots ------------------------ #


def test_add_untracked_pseudoknots_adds_missing_and_skips_known(monkeypatch):
    """
    We stub out design helpers so the test is independent of the real implementation.
    - db_pairs: define symbols that we treat as pseudoknot openers.
    - dot_bracket_to_stacks: returns a 'reduced_db' and matching 'stacks'.
    - dot_bracket_to_pair_map: returns the base-pair map when pair_map is None.
    """
    # 1) Make sure the module under test sees our stubs
    import pyfurnace.generate.pk_utils as mod

    monkeypatch.setattr(mod, "db_pairs", {"A", "B"}, raising=True)

    # We'll present three "symbols": 'A', '.', 'B'. Only A and B should count.
    reduced_db = "A.B"
    stacks = [
        (5, 6),  # A -> to be ADDED unless already in pk_dict
        (42,),  # '.' -> ignored by the function
        (20, 21),  # B -> to be ADDED
    ]
    monkeypatch.setattr(
        mod,
        "dot_bracket_to_stacks",
        lambda s, only_opening=True: (reduced_db, stacks),
        raising=True,
    )

    # Pair map for indices in stacks
    pair_map = {
        5: 95,
        6: 94,
        20: 80,
        21: 79,
        42: 999,  # irrelevant; '.' ignored
    }
    # Ensure the builder is used if pair_map=None
    monkeypatch.setattr(
        mod, "dot_bracket_to_pair_map", lambda s: pair_map, raising=True
    )

    # Existing pk has the (5,6) indexes already (so that one is skipped)
    pk_dict = {
        "pre": {
            "id": "pre",
            "ind_fwd": [(5, 6)],  # use tuples so the set(...) in impl works
            "ind_rev": [(94, 95)],  # existing reverse tuple
            "E": -1.0,
            "dE": 0.1,
        }
    }

    out = add_untracked_pseudoknots(
        pk_dict=pk_dict,
        structure="does-not-matter-here",
        energy=-3.5,
        energy_tolerance=0.7,
        pair_map=None,  # force the stubbed dot_bracket_to_pair_map to be used
    )

    # Only one new PK should be added: for symbol 'B' / indexes (20,21)
    assert set(out.keys()) == {"pre", "extra_0"}

    new_pk = out["extra_0"]
    assert new_pk["id"] == "extra_0"
    # ind_fwd is exactly what comes from stacks
    assert new_pk["ind_fwd"] == [(20, 21)]
    # ind_rev is computed from pair_map using reversed forward order
    # indexes[::-1] = (21, 20) -> [pair_map[21], pair_map[20]] == [79, 80]
    assert new_pk["ind_rev"] == [[79, 80]] or new_pk["ind_rev"] == [
        (79, 80)
    ]  # tolerate list/tuple
    assert new_pk["E"] == -3.5
    assert new_pk["dE"] == 0.7


def test_add_untracked_pseudoknots_uses_provided_pair_map_and_never_calls_builder(
    monkeypatch,
):
    import pyfurnace.generate.pk_utils as mod

    monkeypatch.setattr(mod, "db_pairs", {"X"}, raising=True)
    monkeypatch.setattr(
        mod,
        "dot_bracket_to_stacks",
        lambda s, only_opening=True: ("X", [(100, 101)]),
        raising=True,
    )

    # If the function tries to call this, we want to know.
    def _boom(_):
        raise AssertionError(
            "dot_bracket_to_pair_map should NOT be called when pair_map is provided"
        )

    monkeypatch.setattr(mod, "dot_bracket_to_pair_map", _boom, raising=True)

    supplied_pair_map = {100: 200, 101: 199}

    out = add_untracked_pseudoknots(
        pk_dict={},
        structure="ignored",
        energy=-8.8,
        energy_tolerance=0.2,
        pair_map=supplied_pair_map,  # provided -> builder must not be used
    )

    assert "extra_0" in out
    pk = out["extra_0"]
    assert pk["ind_fwd"] == [(100, 101)]
    assert pk["ind_rev"] == [[199, 200]] or pk["ind_rev"] == [(199, 200)]
    assert pk["E"] == -8.8 and pk["dE"] == 0.2


def test_add_untracked_pseudoknots_ignores_symbols_not_in_db_pairs_or_in_open_paren(
    monkeypatch,
):
    import pyfurnace.generate.pk_utils as mod

    # db_pairs only recognizes 'Q'. We'll feed '.', '(' and 'Q'
    monkeypatch.setattr(mod, "db_pairs", {"Q"}, raising=True)
    reduced_db = ".(Q"
    stacks = [
        (7, 8),  # '.' -> ignored
        (30, 31),  # '(' -> ignored by explicit check
        (50, 51),  # 'Q' -> added
    ]
    monkeypatch.setattr(
        mod,
        "dot_bracket_to_stacks",
        lambda s, only_opening=True: (reduced_db, stacks),
        raising=True,
    )
    monkeypatch.setattr(
        mod,
        "dot_bracket_to_pair_map",
        lambda s: {50: 10, 51: 9, 7: 1000, 8: 999, 30: 900, 31: 899},
        raising=True,
    )

    out = add_untracked_pseudoknots(pk_dict={}, structure="x")
    assert set(out.keys()) == {"extra_0"}
    assert out["extra_0"]["ind_fwd"] == [(50, 51)]


def test_add_untracked_pseudoknots_multiple_new_ids_increment(monkeypatch):
    import pyfurnace.generate.pk_utils as mod

    monkeypatch.setattr(mod, "db_pairs", {"A", "B"}, raising=True)
    monkeypatch.setattr(
        mod,
        "dot_bracket_to_stacks",
        lambda s, only_opening=True: ("AB", [(1, 2), (10, 11)]),
        raising=True,
    )
    monkeypatch.setattr(
        mod,
        "dot_bracket_to_pair_map",
        lambda s: {1: 99, 2: 98, 10: 90, 11: 89},
        raising=True,
    )

    out = add_untracked_pseudoknots(pk_dict={}, structure="y")
    assert set(out.keys()) == {"extra_0", "extra_1"}
    assert out["extra_0"]["ind_fwd"] == [(1, 2)]
    assert out["extra_1"]["ind_fwd"] == [(10, 11)]
