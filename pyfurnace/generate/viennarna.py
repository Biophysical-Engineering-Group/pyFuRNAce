from typing import Tuple
import RNA


def fold(seq: str) -> str:
    """
    Compute the minimum free energy (MFE) structure for an RNA sequence
    using ViennaRNA.

    Parameters
    ----------
    seq : str
        RNA sequence composed of standard bases (A, U, C, G).

    Returns
    -------
    str
        Dot-bracket notation representing the MFE secondary structure.
    """
    # create fold_compound data structure
    fc = RNA.fold_compound(seq)
    # compute MFE and MFE structure
    (mfe_struct, mfe) = fc.mfe()
    return mfe_struct


def fold_p(seq: str) -> Tuple[str, float, float, str, float, float, float]:
    """
    Compute detailed thermodynamic properties of an RNA sequence using ViennaRNA.

    This includes MFE structure, partition function, centroid structure,
    and Boltzmann ensemble statistics.

    Parameters
    ----------
    seq : str
        RNA sequence composed of standard bases (A, U, C, G). For long sequences,
        a small scaling factor is applied for numerical stability.

    Returns
    -------
    tuple
        A 7-tuple containing:
        - str: Dot-bracket notation of MFE structure.
        - float: Free energy of the MFE structure.
        - float: Boltzmann frequency of the MFE structure.
        - str: Dot-bracket notation of the centroid structure.
        - float: Free energy of the centroid structure.
        - float: Boltzmann frequency of the centroid structure.
        - float: Ensemble diversity (mean base pair distance).

    Notes
    -----
    Requires ViennaRNA Python bindings (`import RNA`).
    """
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

    # compute centroid structure
    (centroid_struct, dist) = fc.centroid()

    # compute free energy of centroid structure
    centroid_en = fc.eval_structure(centroid_struct)

    # calculate Boltzmann factors
    frequency_mfe_ensemble = fc.pr_structure(mfe_struct)
    frequency_centr_ensemble = fc.pr_structure(centroid_struct)
    ensemble_diversity = fc.mean_bp_distance()

    return (
        mfe_struct,
        mfe,
        frequency_mfe_ensemble,
        centroid_struct,
        centroid_en,
        frequency_centr_ensemble,
        ensemble_diversity,
    )
