import pytest
from RNA import fold

# Import the public surface (as you said)
from pyfurnace.design import (
    KissingLoop,
    KissingLoop120,
    KissingLoop180,
    BranchedKissingLoop,
    KissingDimer,
    KissingDimer120,
    BranchedDimer,
)

# --- base KissingLoop ---------------------------------------------------------


def test_default_init_and_properties():
    kl = KissingLoop(seq_len=4)  # empty sequence path -> 'N' * 4
    assert kl.pk_index == "0"
    assert kl.energy == -8.5
    assert kl.energy_tolerance == 1.0

    # structure created: one strand with pk_info
    assert len(kl) == 1
    s = kl[0]
    assert hasattr(s, "pk_info")
    assert s.pk_info["id"] == ["0"]
    assert s.pk_info["E"] == [-8.5]
    assert s.pk_info["dE"] == [1.0]


def test_pk_index_validation_and_normalization():
    kl = KissingLoop(seq_len=1)
    # int -> string
    assert kl._check_pk_index(3) == "3"
    # negative int -> abs + "'"
    assert kl._check_pk_index(-2) == "2'"
    # None -> "0"
    assert kl._check_pk_index(None) == "0"
    # bad type
    with pytest.raises(ValueError):
        kl._check_pk_index(3.14)  # neither int nor str


def test_complementary_pk_index():
    kl = KissingLoop(seq_len=1, pk_index="5")
    assert kl.complementary_pk_index() == "5'"
    assert kl.complementary_pk_index(-4) == "4"
    assert kl.complementary_pk_index("7") == "7'"
    assert kl.complementary_pk_index("9'") == "9"


def test_pk_index_setter_rebuilds_strands_and_id_updates():
    kl = KissingLoop(seq_len=2)
    old_id = kl[0].pk_info["id"][0]
    kl.pk_index = "X"
    assert kl.pk_index == "X"
    # strands recreated with new pk id
    assert kl[0].pk_info["id"] == ["X"]
    assert old_id != "X"


def test_energy_setter_updates_both_value_and_pkinfo():
    kl = KissingLoop(seq_len=3)
    kl.energy = -3.14159
    assert kl.energy == -3.14  # rounded in setter
    assert kl[0].pk_info["E"] == [-3.14]


def test_energy_tolerance_setter_validation_and_pkinfo_update():
    kl = KissingLoop(seq_len=2)
    with pytest.raises(ValueError):
        kl.energy_tolerance = -1  # must be positive number

    kl.energy_tolerance = 2.5
    assert kl.energy_tolerance == 2.5
    assert kl[0].pk_info["dE"] == [2.5]


def test_sequence_len_mismatch_raises():
    with pytest.raises(ValueError, match="sequence length doesn't match"):
        KissingLoop(seq_len=3, sequence="AAAA")  # 4 vs 3


def test_sequence_sets_energy_via_rna_fold_when_all_acgu():
    # triggers the `from RNA import fold` path and sets energy + zero tolerance
    kl = KissingLoop(seq_len=4, sequence="ACGU")
    kl = KissingLoop(
        seq_len=4,
        sequence="ACGU",
        strands=kl.strands,
        energy=kl.energy,
        energy_tolerance=kl.energy_tolerance,
    )
    assert kl.sequence == "ACGU"
    assert kl.energy == pytest.approx(fold("ACGU&ACGU")[1], 0.1)  # from fold
    assert kl.energy_tolerance == 0
    # pk_info reflects new energy/tolerance
    assert kl[0].pk_info["E"][0] == pytest.approx(fold("ACGU&ACGU")[1], 0.1)
    assert kl[0].pk_info["dE"] == [0]


def test_get_and_set_sequence_roundtrip():
    kl = KissingLoop(seq_len=3)
    assert str(kl.get_kissing_sequence()) == "NNN"
    kl.set_sequence("ACG")
    assert str(kl.get_kissing_sequence()) == "ACG"


# --- KissingLoop120 -----------------------------------------------------------


def test_kl120_basic_shape_and_len():
    kl = KissingLoop120()  # seq_len fixed to 7 by class
    assert len(kl) == 1
    assert len(kl[0].sequence) == 7


# --- KissingLoop180 -----------------------------------------------------------


def test_kl180_sequence_energy_and_trimmed_getter():
    seq = "ACGUAC"  # length 6 mandated by class
    kl = KissingLoop180(sequence=seq)  # uses fold_180kl=True
    # energy computed via RNA.fold (the 180 path still uses fold)
    assert kl.energy == pytest.approx(fold("AAACGUACA&AAGUACGUA")[1], 0.1)  # from fold
    assert kl.energy_tolerance == 0

    # get_kissing_sequence trims [2:-1]
    trimmed = str(kl.get_kissing_sequence())
    assert trimmed == seq


def test_kl180_overridden_create_strands_affects_pk_info_indices():
    kl = KissingLoop180(sequence="ACGUAC", pk_index=1)
    # overridden _create_strands sets pk_info["ind_fwd"] = [(2, 7)]
    assert kl[0].pk_info["ind_fwd"] == [(2, 7)]
    assert kl[0].pk_info["id"] == ["1"]


# --- BranchedKissingLoop ------------------------------------------------------


def test_branched_kissing_loop_two_strands_and_getter():
    bkl = BranchedKissingLoop()  # seq_len fixed to 6
    assert len(bkl) == 2  # main + connector
    # getter finds strand with len 7 and drops last base
    ks = str(bkl.get_kissing_sequence())
    assert len(ks) == 6


# --- KissingDimer (180 dimer) -------------------------------------------------


def test_kissing_dimer_builds_two_complementary_strands():
    kd = KissingDimer(sequence="ACGUAC", pk_index="Z")
    assert len(kd) == 2
    top, bottom = kd[0], kd[1]

    # top uses pk_index as given; bottom created with complementary pk index
    assert top.pk_info["id"] == ["Z"]
    assert bottom.pk_info["id"] == ["Z'"]


def test_kissing_dimer_set_up_and_down_sequence_updates():
    kd = KissingDimer(sequence="ACGUAC", pk_index="Z")
    new_seq = "UGCAUG"
    kd.set_up_sequence(new_seq)
    assert kd[0].sequence == "AA" + new_seq + "A"  # respects KissingLoop180 layout

    kd.set_down_sequence(new_seq)  # setter takes reverse complement internally
    # bottom is stored as AA + revcomp(top core) + A
    assert str(kd[1].sequence).startswith("AA") and str(kd[1].sequence).endswith("A")
    assert len(kd[1].sequence) == 9  # AA + 6 + A


# --- KissingDimer120 ----------------------------------------------------------


def test_kissing_dimer_120_builds_two_strands():
    kd120 = KissingDimer120(sequence="ACGUUCA", pk_index="Q")  # len 7
    assert len(kd120) == 2
    top, bottom = kd120[0], kd120[1]
    assert top.pk_info["id"] == ["Q"]
    assert bottom.pk_info["id"] == ["Q'"]


# --- BranchedDimer ------------------------------------------------------------


def test_branched_dimer_three_strands_and_offsets():
    bd = BranchedDimer(sequence="ACGUAC", pk_index="B")
    assert len(bd) == 3  # top + branched KL (2 strands)

    # sanity checks on pk ids
    assert bd[0].pk_info["id"] == ["B"]
    # one of the other strands (branched side) holds complementary id
    ids = [
        s.pk_info.get("id", [None])[0]
        for s in bd
        if hasattr(s, "pk_info") and s.pk_info
    ]
    assert "B'" in ids
