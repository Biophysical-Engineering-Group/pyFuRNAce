import sys
import RNA


def distance(seq1, seq2):
    print(RNA.bp_distance(seq1, seq2))


def duplex(seq1, seq2):
    dot_bracket, energy = RNA.fold(f"{seq1}&{seq2}")
    dot_bracket = dot_bracket[: len(seq1)] + "&" + dot_bracket[len(seq1) :]
    print(dot_bracket, round(energy, 2))


def ensemble_diversity(seq):
    # create model details
    md = RNA.md()
    # adjust the scaling factor for long sequences
    if len(seq) > 1000:
        md.sfact = 1.01
    # create fold_compound data structure
    # (required for all subsequently applied  algorithms)
    fc = RNA.fold_compound(seq, md)
    # compute MFE and MFE structure
    (mfe_struct, mfe) = fc.mfe()
    # rescale Boltzmann factors for partition function computation
    fc.exp_params_rescale(mfe)
    # compute partition function (NECESSARY STEP)
    (pp, pf) = fc.pf()
    # calculate Boltzmann factors
    ensemble_diversity = fc.mean_bp_distance()
    print(f"Ensemble diversity: {round(ensemble_diversity, 2)}")


def fold(seq):
    print(RNA.fold(seq)[0])


if __name__ == "__main__":
    run_function = sys.argv[1]
    match run_function:
        case "distance":
            distance(sys.argv[2], sys.argv[3])
        case "duplex":
            duplex(sys.argv[2], sys.argv[3])
        case "diversity":
            ensemble_diversity(sys.argv[2])
        case "fold":
            fold(sys.argv[2])
