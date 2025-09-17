import random
import pytest
from pyfurnace.design import Sequence, AmbiguosStructure

# --------------------
# Core construction & dunders
# --------------------


def test_init_and_str_repr_len_bool_iter():
    s = Sequence("ACGU")
    assert str(s) == "ACGU"
    assert repr(s) in {"5 ACGU 3", "3 ACGU 5"}  # default '53' -> "5 ACGU 3"
    assert len(s) == 4
    assert bool(s) is True
    assert list(iter(s)) == list("ACGU")


def test_invalid_directionality_raises_on_init_and_setter():
    with pytest.raises(ValueError):
        Sequence("ACGU", directionality="xx")  # type: ignore

    s = Sequence("ACGU")
    with pytest.raises(ValueError):
        s.directionality = "99"  # type: ignore


def test_getitem_and_directionality_with_negative_step():
    s = Sequence("AUGC", directionality="53")
    sub1 = s[1:3]  # "UG", dir '53'
    assert str(sub1) == "UG"
    assert sub1.directionality == "53"
    sub2 = s[3:0:-1]  # "CGU" (reverse slice), dir flips to '35'
    assert str(sub2) == "CGU"
    assert sub2.directionality == "35"


def test_setitem_index_and_slice():
    s = Sequence("AUGC")
    s[0] = "C"
    assert str(s) == "CUGC"
    s[1:3] = "AA"
    assert str(s) == "CAAC"
    s[-1] = "U"
    assert str(s) == "CAAU"


def test_add_iadd_radd_mul_and_mismatch_addition():
    s = Sequence("AU", "53")
    t = Sequence("GC", "53")
    assert str(s + t) == "AUGC"
    s2 = Sequence("AU", "53")
    s2 += "GC"
    assert str(s2) == "AUGC"

    assert str("GG" + Sequence("AA")) == "GGAA"
    assert str(3 * Sequence("AU")) == "AUAUAU"
    assert str(Sequence("AU") * 2) == "AUAU"
    with pytest.raises(ValueError):
        _ = Sequence("AU") * 2.5  # type: ignore

    # different directionality cannot add
    with pytest.raises(ValueError):
        _ = Sequence("AU", "53") + Sequence("GC", "35")


def test_contains_and_eq_respect_directionality():
    host = Sequence("AUGC", "53")
    # other with opposite directionality: reversed string is tested for containment
    other = Sequence("GU", "35")  # reversed -> "UG" which *is* in host
    assert other in host

    # equality: same content and directionality
    assert Sequence("ACGU", "53") == Sequence("ACGU", "53")
    # equality with opposite directionality compares reversed
    assert Sequence("AUGC", "53") == Sequence("CGUA", "35")
    # equality with plain str ignores directionality
    assert Sequence("ACGU", "35") == "ACGU"
    assert (Sequence("ACGU") == 123) is False
    assert "UGC" in host
    assert "CUG" not in host


def test_hash_is_stable_and_uses_repr():
    s1 = Sequence("ACGU", "53")
    s2 = Sequence("ACGU", "53")
    assert hash(s1) == hash(s2)


# --------------------
# Biochem utilities
# --------------------


def test_complement_and_reverse_complement():
    s = Sequence("AUGC")
    assert str(s.complement()) == "UACG"  # A-U, U-A, G-C, C-G
    assert (
        str(s.reverse_complement()) == str(s.complement())[::-1]
    )  # equals complement reversed
    # Explicit check
    assert str(s.reverse_complement()) == "GCAU"


def test_copy_is_independent():
    s = Sequence("ACGU")
    c = s.copy()
    assert str(c) == "ACGU"
    c[0] = "U"
    assert str(s) == "ACGU"
    assert str(c) == "UCGU"


def test_distance_basic_and_errors():
    s = Sequence("AUGC")
    # identical -> distance 0
    assert s.distance("AUGC") == 0
    # 1 mismatch
    assert s.distance("AUGU") == 1
    # length mismatch -> error
    with pytest.raises(ValueError):
        s.distance("AU")


def test_find_and_split_lower_upper_translate():
    with pytest.warns(
        AmbiguosStructure,
        match="Warning: The string 'AUGC GC' contains nucleotides not"
        " allowed in ROAD that will be removed.",
    ):
        s = Sequence("AUGC GC")
    assert str(s) == "AUGCGC"
    assert s.find("GC") == 2
    parts = s.split("U")
    assert [str(p) for p in parts] == ["A", "GCGC"]

    assert s.lower() == "augcgc"
    assert s.upper() == "AUGCGC"

    t = s.translate({"A": "G"})  # remove spaces, non-inplace
    assert str(t) == "GUGCGC"
    assert str(s) == "AUGCGC"  # original unchanged

    s.translate({"G": "A"}, inplace=True)
    assert str(s) == "AUACAC"
    # More robust direct check: every 'G' became 'A'
    assert "G" not in str(s)


def test_gc_content_basic_and_extended():
    s = Sequence("GGCCAAUU")  # 4/8 are GC -> 50%
    assert pytest.approx(s.gc_content(), rel=1e-6) == 50.0

    # Extended alphabet: 'S' counts fully as GC, 'N' contributes 1/4
    s2 = Sequence(
        "SSNN"
    )  # extended adds 1 (S) * 2 + 0.25 (N) * 2 = 2.5 GC out of 4 -> 62.5%
    assert pytest.approx(s2.gc_content(extended_alphabet=True), rel=1e-6) == 62.5


def test_molecular_weight_sum():
    # from table in code: A=347.2, G=363.2, C=323.2, U=324.2
    s = Sequence("AGCU")
    expected = 347.2 + 363.2 + 323.2 + 324.2
    assert pytest.approx(s.molecular_weight(), rel=1e-6) == expected


def test_pop_and_replace_and_reverse():
    s = Sequence("AUGC")
    popped = s.pop(1)
    assert popped == "U"
    assert str(s) == "AGC"

    s.replace("G", "U")
    assert str(s) == "AUC"

    # reverse (not inplace): directionality flips, sequence unchanged
    r = s.reverse(inplace=False)
    assert str(r) == str(s)
    assert r.directionality == "35"
    # reverse (inplace): flips directionality on self
    s.reverse(inplace=True)
    assert s.directionality == "35"


def test_directionality_property_and_concat_callbacks_noop():
    s = Sequence("ACGU", directionality="53")
    s = 0 + s  # radd with int no-op
    assert str(s) == "ACGU"
    s.directionality = "35"
    assert s.directionality == "35"
    # __iadd__ triggers callbacks (stubbed as noop); ensure sequence updates
    s += "AC"
    assert str(s) == "ACGUAC"


# --------------------
# Randomization & structure-aware pairing
# --------------------


def test_get_random_sequence_length_and_chars_no_structure():
    s = Sequence("NNNN")
    rnd = s.get_random_sequence()
    assert len(rnd) == len(s)
    assert set(str(rnd)) <= set("ACGU")  # N expands to A/C/G/U
    assert rnd.directionality == s.directionality


def test_get_random_sequence_raises_on_length_mismatch_structure():
    s = Sequence("NN")
    with pytest.raises(ValueError):
        s.get_random_sequence(structure="()()")  # length mismatch


def test_get_random_sequence_respects_pairing_basic(monkeypatch):
    # Set choice to deterministic but still allow pairing logic to run
    def _first(xs):
        return sorted(xs)[0]

    monkeypatch.setattr(random, "choice", _first)

    # Dot-bracket "()" pairs 0 with 1
    s = Sequence("NN")
    rnd = s.get_random_sequence(structure="()")
    seq = str(rnd)
    # Check complementary pairing using Sequence's own complement mapping
    comp0 = Sequence(seq[0]).complement()
    assert str(comp0) == seq[1]  # e.g., A <-> U or C <-> G
