import streamlit as st
from streamlit_option_menu import option_menu
import re
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import gc_fraction, molecular_weight
import warnings
from utils import load_logo, main_menu_style, copy_to_clipboard
from utils.template_functions import symbols, write_format_text, check_dimer, sanitize_input, reference

def convert_tab(seq):
    """ Calculate the main property of the sequence (gc content, molecular weight, nucleotide composition, melting temperature) 
        and convert the sequence to the different formats: reverse, complement, reverse complement, RNA transcribed, DNA template.
    """
    # format the input sequence
    # st.write(f"Your sequence ({len(seq)} nt):")
    st.write(f'Sequence type: {seq_type}. Length: {len(seq)} bases')
    # write_format_text(seq)
        
    # calculate the main properties of the sequence
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.write("GC content (%)")
        write_format_text(round(gc_fraction(seq, ambiguous='ignore')*100, 2))
    with col2:
        st.write("Molecular weight (Dalton)")
        write_format_text(f'{molecular_weight(seq, seq_type):.3e}')
    with col3:
        st.write("Nucleotide composition")
        bases = set(str(seq))
        t = ''
        for b in bases:
            t = t + f"{b}: {str(seq).count(b)}; "
        write_format_text(t)
    with col4:
        st.write("Melting temperature (°C)")
        if seq_type == 'RNA':
            nn_table = mt.RNA_NN3
            complement = seq.complement_rna()
        else:
            nn_table = mt.DNA_NN4
            complement = seq.complement()
        write_format_text(round(mt.Tm_NN(seq, 
                                         c_seq=complement, 
                                         nn_table=nn_table),2))
    
    # convert the sequence to the different formats: reverse, complement, reverse complement, RNA transcribed (for DNA), DNA template (for RNA)
    col1, col2, col3 = st.columns(3)
    with col1:
        copy_to_clipboard(seq[::-1], "Reverse:")
        write_format_text(seq[::-1])
    with col2:
        copy_to_clipboard(complement, "Complement:")
        write_format_text(complement)
    with col3:
        if seq_type == 'DNA':
            rev_complement = seq.reverse_complement()
        else:
            rev_complement = seq.reverse_complement_rna()
        copy_to_clipboard(rev_complement, "Reverse complement:")
        write_format_text(rev_complement)

    if seq_type == 'DNA':
        st.write("RNA transcribed:")
        write_format_text(seq.transcribe())

    # if the sequence is RNA, create a DNA template 
    else:
        col1, col2 = st.columns([1, 5])
        # st.write("#### DNA template:")
        promoter = st.text_input("**Select a promoter** (default: T7 promoter)", value='TAATACGACTCACTATA')
        # take a specific promoter sequence for the DNA template and check that the RNA sequence starts with G
        if promoter == 'TAATACGACTCACTATA' and seq[0] != 'G':
            st.warning("The RNA sequence doesn't start with G", icon="⚠️")
        dna_template = Seq(promoter) + seq.back_transcribe()
        col1, col2 = st.columns([1, 5])    

        # coding strand
        col1, col2 = st.columns([1, 9], vertical_alignment='center')
        with col1:
            copy_to_clipboard(dna_template, "Coding strand:")
        with col2:
            write_format_text(dna_template)

        ### Save the DNA template in the session state and add a link to the primer page
        st.session_state["dna_template"] = str(dna_template)
        st.session_state.prepare_ind = 0
        st.page_link("pages/4_Prepare.py", label=":orange[Prepare the PCR primers]", icon=":material/sync_alt:")

def search_tab(seq):
    """ Search a subsequence and highlight the bases"""
    subseq_list = []
    i = 0
    # take a subsequence
    subseq = sanitize_input(st.text_input("Search for the subsequence:", key=f'seq{i}'))
    subseq_list.append(subseq)
    while subseq_list[-1]:
        # take the last subsequence
        sub_seq = subseq_list[-1]
        # calculate the indexes of the subsequence in the sequence
        indexes = [substr.start()+1 for substr in re.finditer(str(sub_seq), str(seq))]
        if not indexes:
            st.error('Subsequence not found', icon=":material/personal_injury:")
        else:
            # highlight in red the subsequence in the sequence
            highlighted = str(seq).replace(sub_seq, f':red[{sub_seq}]')
            st.write(f"Sequence found at index: {indexes}")
            st.markdown(highlighted)
        # add the space for another subsequence
        st.divider()
        i+=1 # add another subsequence with a different widget key
        subseq = sanitize_input(st.text_input("Search for the subsequence:", key=f'seq{i}'))
        # add the subsequence to the subsequence list
        subseq_list.append(subseq)
        # if the last subsequences is empty, stop the loop

def dimer_tab(seq):
    """ Check the dimer between the main sequence and a list of subsequences and highlight the aligned bases"""
    subseq_list = []
    i = 0
    # take a subsequence
    subseq = sanitize_input(st.text_input("Check dimer with the sequence:", key=f'dim{i}'))
    subseq_list.append(subseq)
    while subseq_list[-1]:
        # take the last subsequence
        sub_seq = subseq_list[-1]
        st.write('Best dimer found:')
        # calculate the best dimer between the sequence and the subsequence
        dimer_subseq = check_dimer(seq, sub_seq, basepair=basepair)
        # write the dimer that is formed
        write_format_text(dimer_subseq)
        # add expander with all dimers found for the subsequence
        with st.expander("All dimers"):
            # calculate all the dimers between the sequence and the subsequence
            dimer_dict = check_dimer(seq, sub_seq, dict_format=True, basepair=basepair)
            # add a slider to select the number of basepairs   
            n_pairs = st.slider('Dimer with n° of basepairs:', min_value=min(dimer_dict), max_value=max(dimer_dict), value=max(dimer_dict), key=f'dimers{i}')
            if n_pairs in dimer_dict:
                for other_dimer in dimer_dict[n_pairs]:
                    write_format_text(other_dimer)
            else: st.write(f"No Dimer found for {n_pairs} basepairs")
        # add the space for another subsequence
        st.divider()
        i+=1
        subseq = sanitize_input(st.text_input("Search for the subsequence:", key=f'dim{i}'))
        subseq_list.append(subseq)
        # if the last subsequences is empty, stop the loop

if __name__ == "__main__":
    ### set the logo of the app
    load_logo()
    warnings.filterwarnings("ignore") # ignore warnings
    
    if "rna_origami_seq" not in st.session_state:
        st.session_state["rna_origami_seq"] = ''
    st.header('Convert', help='Prepare the DNA template for you RNA Origami, align structures and search for dimers.')
    # take the input sequence and sanitize it
    seq = sanitize_input(st.text_input("Input sequence (DNA or RNA):", value = st.session_state["rna_origami_seq"]))
    # check the symbols in the sequence
    if set(seq) - symbols:
        st.warning('The sequence contains symbols not included in the [IUPAC alphabet](https://www.bioinformatics.org/sms/iupac.html).', icon="⚠️")
    
    if seq != st.session_state["rna_origami_seq"]:
        st.session_state["rna_origami_seq"] = seq
        st.rerun()
    
    seq = Seq(seq)
    seq_type = 'DNA'
    if 'T' in seq and 'U' in seq:
        st.error('Both T and U found in the sequence', icon=":material/personal_injury:")
        # create a porper biopython symbol
    elif 'T' in seq : 
        seq_type = 'DNA'
    elif "U" in seq:
        seq_type = 'RNA'
    elif seq:
        seq_type = st.radio(
            "Choose the sequence type: (usually auto-detected)",
            ["DNA", "RNA"],
            horizontal = True,
            )

    # according to the sequence type, create the basepair dictionary
    if seq_type == 'DNA':
        basepair = str.maketrans('TACG', 'ATGC')
    elif seq_type == 'RNA':
        basepair = str.maketrans('UACG', 'AUGC')
    
    # create the tabs with the functions
    st.write("\n") # add space between initial menu and motif menu
    option_data = {'Convert': "bi bi-arrow-repeat",
                   'Search': "bi bi-align-center",
                   'Dimer': "bi bi-bar-chart-steps"}

    selected_operation = option_menu(None, 
                                    list(option_data.keys()),
                                    icons=list(option_data.values()),
                                    menu_icon="cast", 
                                    orientation="horizontal",
                                    styles=main_menu_style)

    if not seq:
        st.stop()
    elif selected_operation == "Convert":
        convert_tab(seq)
    elif selected_operation == "Search":
        search_tab(seq)
    elif selected_operation == "Dimer":
        dimer_tab(seq)
    
    # add bibliography
    reference()