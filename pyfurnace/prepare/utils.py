from typing import Union, Tuple, List, Callable
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import gc_fraction
from functools import partial

# https://www.bioinformatics.org/sms/iupac.html
# dictionary of melting temperature methods:
# name as key and function as value
tm_methods = {
    "Nearest Neighbor": mt.Tm_NN,
    "Empirical GC content": mt.Tm_GC,
    "Wallace, rule of thumb": mt.Tm_Wallace,
}

# dictionary of melting temperature models:
#  name as key and function as value
tm_models = {
    "Nearest Neighbor": [
        "DNA_NN4",
        "DNA_NN1",
        "DNA_NN2",
        "DNA_NN3",
        "RNA_NN1",
        "RNA_NN2",
        "RNA_NN3",
        "R_DNA_NN1",
    ],
    "Empirical GC content": [1, 2, 3, 4, 5, 6, 7, 8],
    "Wallace, rule of thumb": [],
}

# dictionary of Nearest Neighbour melting temperature models:
#  name as key and function as value
NN_models = {
    "DNA_NN4": mt.DNA_NN4,
    "DNA_NN1": mt.DNA_NN1,
    "DNA_NN2": mt.DNA_NN2,
    "DNA_NN3": mt.DNA_NN3,
    "RNA_NN1": mt.RNA_NN1,
    "RNA_NN2": mt.RNA_NN2,
    "RNA_NN3": mt.RNA_NN3,
    "R_DNA_NN1": mt.R_DNA_NN1,
}

# dictionary of default values for the primer energy parameters

default_values = {
    "Na": 0,
    "K": 50,
    "Tris": 20,
    "Mg": 1.5,
    "dNTPs": 0.2,
    "Method": 7,
    "DMSO (%)": 0,
    "mt_method": list(tm_methods)[0],
    "mt_model": tm_models[list(tm_methods)[0]][0],
    "Primer": 500,
}


def make_tm_calculator(
    method_name: str = "Nearest Neighbor",
    model_name: str = "DNA_NN4",
    primer_conc: float = 500,
    tm_kwargs: dict = default_values.copy(),
) -> callable:
    """
    Creates a melting temperature (Tm) calculator function for DNA primers using
    specified calculation method and model.
    Parameters
    ----------
    method_name : str, optional
        The method used for Tm calculation.
        Options include "Nearest Neighbor",
        "Empirical GC content",
        and "Wallace, rule of thumb".
        Default is "Nearest Neighbor".
    model_name : str, optional
        The nearest neighbor model to use for Tm calculation.
        Only relevant if `method_name` is "Nearest Neighbor".
        Default is "DNA_NN4".
    primer_conc : float, optional
        The concentration of the primer (in nM) used in the calculation.
        Only relevant for "Nearest Neighbor" method.
        Default is 500.
    tm_kwargs : dict, optional
        Additional keyword arguments for the Tm calculation method.
        May include chemical corrections (e.g., DMSO, dNTPs).
        Default is a copy of `default_values`.
    Returns
    -------
    callable
        A function that takes a DNA sequence (str) and returns its melting
        temperature (float), possibly corrected for chemical additives.

    Notes
    -----
    - If DMSO is specified in `tm_kwargs`, the returned function
        applies a chemical correction to the calculated Tm.
    - For "Nearest Neighbor" method, dNTPs and DMSO are removed
        from the keyword arguments before calculation.
    - For "Empirical formulas based on GC content",
        DMSO is removed from the keyword arguments before calculation.
    - For "Wallace" method, no additional arguments are processed.
    Examples
    --------
    >>> tm_calc = make_tm_calculator(method_name="Nearest Neighbor", primer_conc=250)
    >>> tm_calc("ATCGATCG")
    62.5
    """
    method = tm_methods[method_name]

    first_word = method_name.split()[0].lower()

    if first_word == "nearest":
        model = NN_models[model_name]
        method_kwargs = tm_kwargs.copy()
        if "dNTPs" in method_kwargs:
            method_kwargs.pop("dNTPs")
        if "DMSO" in method_kwargs:
            method_kwargs.pop("DMSO")
        partial_func = partial(
            method, nn_table=model, dnac1=primer_conc, dnac2=0, **method_kwargs
        )

    elif first_word == "empirical":
        method_kwargs = tm_kwargs.copy()
        if "DMSO" in method_kwargs:
            method_kwargs.pop("DMSO")
        partial_func = partial(method, valueset=model_name, **method_kwargs)

    else:  # Wallace
        partial_func = method

    if tm_kwargs.get("DMSO", 0) > 0:
        return lambda seq: mt.chem_correction(partial_func(seq), DMSO=tm_kwargs["DMSO"])

    return partial_func


def calculate_gc(seq: str) -> float:
    """
    Calculate the GC content percentage of a nucleotide sequence.

    Parameters
    ----------
    seq : str
        The nucleotide sequence for which to calculate the GC content.

    Returns
    -------
    float
        The GC content of the sequence as a percentage, rounded to one decimal place.

    Notes
    -----
    Ambiguous nucleotides are ignored in the calculation.
    """
    return round(gc_fraction(seq, ambiguous="ignore") * 100, 1)


def annealing_temp(
    mts: List[float], seq: str, tm_kwargs: dict, method: str = "IDT"
) -> float:
    """
    Calculate the annealing temperature for PCR primers using either
    the IDT or Phusion method.

    Parameters
    ----------
    mts : List[float]
        List of melting temperatures (Tm) for the primers.
    seq : str
        DNA sequence for which the annealing temperature is to be calculated.
    tm_kwargs : dict
        Dictionary of keyword arguments for Tm calculation, such as
        salt concentration, DMSO, etc.
    method : str, optional
        Method to use for calculation. Either 'IDT' or 'Phusion'.
        Default is 'IDT'.

    Returns
    -------
    float
        Calculated annealing temperature, rounded to two decimal places.

    Notes
    -----
    - The 'IDT' method applies a chemical correction and a weighted average formula.
    - The 'Phusion' method uses the minimum Tm, with an adjustment for primer length.
    - Certain keys in `tm_kwargs` (e.g., 'dNTPs', 'DMSO') are handled specifically
        for the 'IDT' method.
    """
    if method == "IDT":
        method_kwargs = tm_kwargs.copy()
        if "dNTPs" in method_kwargs:
            method_kwargs.pop("dNTPs")
        if "DMSO" in method_kwargs:
            method_kwargs.pop("DMSO")
        t_anneal = (
            0.3 * min(mts)
            + 0.7
            * mt.chem_correction(
                mt.Tm_GC(seq, valueset=7, **method_kwargs), DMSO=tm_kwargs["DMSO"]
            )
            - 14.9
        )
    else:  # Phusion method
        t_anneal = min(mts)
        if all(len(p) > 20 for p in seq):
            t_anneal += 3
    return round(t_anneal, 2)


def check_dimer(
    seq1: str,
    seq2: str,
    dict_format: bool = False,
    basepair: dict = {"A": "T", "T": "A", "C": "G", "G": "C"},
) -> Union[str, dict]:
    """
    Check the dimerization of two sequences and return the best dimer found.
    The dimerization is checked in an extremely simple way, by aligning the
    two sequences and checking for WC basepairing at the same index.
    If dict_format is True, return a dictionary with all the dimers found.

    Parameters:
    -----------
    seq1 (str):
        The first sequence to check.
    seq2 (str):
        The second sequence to check.
    dict_format (bool):
        If True, return a dictionary with all the dimers found.
    basepair (dict):
        A dictionary with the basepairing rules.
        Default is {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'} for RNA.
        For DNA, use {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}.

    Returns:
    --------
    str or dict:
        If dict_format is True, return a dictionary with the dimers found.
        If dict_format is False, return the best dimer found as a string.
    """
    dimers_dict = {}
    basepair = str.maketrans(basepair)
    top_seq = " " * (len(seq2) - 1) + str(seq1) + " " * len(seq2)
    bottom_seq = str(seq2[::-1])
    max_match = [top_seq, " " * len(top_seq), bottom_seq]

    for start_ind in range(len(seq1) + len(seq2) - 2):
        bp = ""
        for ind, b in enumerate(bottom_seq):
            if (
                b.translate(basepair) == top_seq[ind]
                and b != " "
                and top_seq[ind] != " "
            ):
                bp += "┊"
            else:
                bp += " "
        num_bp = bp.count("┊")
        index_num = f"[{max(0, start_ind - len(seq2) + 1)}]"
        if num_bp > max_match[1].count("┊"):
            max_match = [
                index_num + top_seq,
                " " * len(index_num) + bp,
                " " * len(index_num) + bottom_seq,
            ]
        if dict_format:
            if num_bp not in dimers_dict:
                dimers_dict[num_bp] = []
            dimers_dict[num_bp].append("\n".join(max_match))
        top_seq = top_seq[1:]

    return dimers_dict if dict_format else "\n".join(max_match)


def auto_design_primers(
    seq: Seq, target_temp: float = 65.0, tm_func: Callable = make_tm_calculator()
) -> Tuple[List[str], List[float]]:
    """
    Automatically designs forward and reverse primers for a given DNA sequence,
      optimizing for melting temperature (Tm) and primer quality.

    Parameters
    ----------
    seq : Seq
        The DNA sequence for which primers are to be designed.
    target_temp : float, optional
        The target melting temperature (Tm) for the primers, by default 65.0°C.
    tm_func : Callable, optional
        A function to calculate the melting temperature of a primer sequence,
        by default uses `make_tm_calculator()`.

    Returns
    -------
    final_primers : list of str
        A list containing the best forward and reverse primer sequences.
    final_mts : list of float
        A list containing the melting temperatures (Tm) of the selected primers.

    Notes
    -----
    - The function evaluates primer candidates based on Tm, GC content, length,
        and dimerization potential.
    - Returns empty string and 0 Tm if no suitable primer is found for a direction.
    """
    final_primers = []
    final_mts = []

    for direction in (1, -1):
        prim_length = 10
        primers_info = []

        to_prime = seq if direction == 1 else str(Seq(seq).reverse_complement())

        while prim_length < len(seq):
            primer = to_prime[:prim_length]
            tm = round(tm_func(primer), 2)

            if tm > (target_temp + 2.5):
                break
            if tm < (target_temp - 2.5):
                prim_length += 1
                continue

            score = 0
            if primer[-1] in "GC":
                score += 1
            if primer[-2] in "GC":
                score += 1
            if 18 < prim_length < 30:
                score += 1
            elif prim_length < 18:
                score -= 17 - prim_length
            if 40 < gc_fraction(primer, ambiguous="ignore") * 100 < 60:
                score += 1
            score -= abs(tm - target_temp) / 2

            dimer_score = max(check_dimer(primer, primer, dict_format=True).keys())
            score -= dimer_score / (len(primer) // 2)

            primers_info.append((score, tm, str(primer)))
            prim_length += 1

        primers_info.sort(reverse=True)
        if primers_info:
            final_mts.append(primers_info[0][1])
            final_primers.append(primers_info[0][2])
        else:
            final_primers.append("")
            final_mts.append(0)

    return final_primers, final_mts
