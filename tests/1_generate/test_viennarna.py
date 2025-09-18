from pyfurnace.generate.viennarna import fold_p


def test_fold_p_long_sequence():
    seq = "G" * 300 + "A" * 298 + "GGUUCGCC" + "U" * 298 + "C" * 300
    # no warnings or errors
    data = fold_p(seq)
    (
        mfe_struct,
        mfe,
        frequency_mfe_ensemble,
        centroid_struct,
        centroid_en,
        frequency_centr_ensemble,
        ensemble_diversity,
    ) = data

    assert len(mfe_struct) == len(seq)
    assert len(centroid_struct) == len(seq)
    assert mfe_struct == "(" * 600 + "...." + ")" * 600
    assert centroid_struct == mfe_struct
    assert frequency_mfe_ensemble > 0.01
    assert frequency_centr_ensemble > 0.01
    assert mfe < 0
    assert ensemble_diversity > 1
