from pyfurnace import Coords


def test_helix_caching():
    pos = (0, 0, 0)
    b = (0, 0, 1)
    n = (1, 0, 0)
    # the first cached helix is the kissing loop helix
    stem1 = Coords.compute_helix_from_nucl(pos, b, n, 8, double=False)

    assert Coords._CACHED_AE_T is None
    assert Coords._CACHED_AE_T_INV is None
    assert len(Coords._CACHED_HELICES) == 2
    cached_helices = iter(Coords._CACHED_HELICES.values())
    assert len(next(cached_helices)) == 9 * 2  # KL cached helix, double stranded
    assert len(next(cached_helices)) == len(
        stem1
    )  # stem1 cached helix, single stranded

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
