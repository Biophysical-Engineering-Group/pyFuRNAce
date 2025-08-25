import streamlit as st

symbols = {
    "A",
    "T",
    "U",
    "C",
    "G",
    "R",
    "Y",
    "S",
    "W",
    "K",
    "M",
    "B",
    "D",
    "H",
    "V",
    "N",
    "5",
    "3",
}


def sanitize_input(text):
    """Remove spaces and convert to upper case the text input"""
    return text.upper().replace(" ", "")


def upload_setting_button():
    """Allow to upload setting"""
    st.session_state["upload_setting"] = True
    return


def check_dimer(
    seq1: str,
    seq2: str,
    dict_format: bool = False,
    basepair={"A": "U", "U": "A", "C": "G", "G": "C"},
):
    """
    Check the dimer between two sequences and return the best dimer found.
    If dict_format is True, return a dictionary with all the dimers found."""
    dimers_dict = {}
    basepair = str.maketrans(basepair)

    # add spaces to the top sequence so it can be aligned with the bottom sequence
    top_seq = " " * (len(seq2) - 1) + str(seq1) + " " * len(seq2)

    # take the second sequence and reverse it to align it with the first sequence
    bottom_seq = str(seq2[::-1])

    # initialize the best dimer found
    max_match = [top_seq, " " * len(top_seq), bottom_seq]

    # iterate over the sequences to find the best dimer
    for start_ind in range(len(seq1) + len(seq2) - 2):
        bp = ""  # initialize the basepair line

        # iterate over the bottom sequence to check simple basepairing
        #  with the top sequence at the same index
        for ind, b in enumerate(bottom_seq):
            if (
                b.translate(basepair) == top_seq[ind]
                and b != " "
                and top_seq[ind] != " "
            ):
                bp = bp + "┊"
            else:
                bp = bp + " "

        # count the basepairs found
        num_bp = bp.count("┊")

        # create a index_num variable, which is the minimum index
        # of the top sequence that is aligned with the bottom sequence
        if start_ind - len(seq2) + 1 <= 0:
            index_num = "[0]"
        else:
            index_num = f"[{start_ind - len(seq2) + 1}]"

        # check if the current dimer is better than the previous best dimer
        if num_bp > max_match[1].count("┊"):
            max_match = [
                index_num + top_seq,
                " " * len(index_num) + bp,
                " " * len(index_num) + bottom_seq,
            ]

        # add the dimer to the dictionary with the number of basepairs as key
        if dict_format:
            if num_bp not in dimers_dict:
                dimers_dict[num_bp] = list()
            dimers_dict[num_bp].append(
                "\n".join(
                    [
                        index_num + top_seq,
                        " " * len(index_num) + bp,
                        " " * len(index_num) + bottom_seq,
                    ]
                )
            )
        top_seq = top_seq[1:]

    # if the dictionary format is selected, return the dictionary
    if dict_format:
        return dimers_dict

    # if the dictionary format is not selected, return the best dimer found
    return "\n".join(max_match)


def reference(primers=False):
    """Add the references to the app"""
    for _ in range(4):
        st.write("\n")
    with st.popover("References"):
        # st.markdown("**References**")
        st.markdown(
            "[1] [Biopython documentation](https://biopython.org/docs"
            "/1.75/api/Bio.SeqUtils.MeltingTemp.html). Cock, P.J.A. et"
            " al. Biopython: freely available Python tools for computational"
            " molecular biology and bioinformatics. Bioinformatics 2009 "
            " Jun 1; 25(11) 1422-3 https://doi.org/10.1093/bioinformatics/"
            "btp163 pmid:19304878"
        )
        if primers:
            st.markdown(
                "[2] [IDT suggested annealing temperature]"
                "(https://eu.idtdna.com/pages/support/faqs/how-do-you"
                "-calculate-the-annealing-temperature-for-pcr-). "
                "Rychlik W, Spencer WJ, Rhoads RE. Optimization of the"
                " annealing temperature for DNA amplification in vitro."
                " Nucleic Acids Res. 1990;18(21):6409-6412."
            )
            st.markdown(
                "[3] [Phusion suggested annealing temperature](https://"
                "www.thermofisher.com/de/de/home/brands/thermo-scientific/"
                "molecular-biology/molecular-biology-learning-center/"
                "molecular-biology-resource-library/spotlight-articles/"
                "optimizing-tm-and-annealing.html)"
            )
