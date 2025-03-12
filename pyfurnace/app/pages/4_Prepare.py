import streamlit as st
from pathlib import Path
from functools import partial
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import gc_fraction
from streamlit_option_menu import option_menu
import json
import warnings
from pathlib import Path
### import the template functions
from utils import load_logo, main_menu_style, second_menu_style, copy_to_clipboard
from utils.template_functions import symbols, write_format_text, check_dimer, reference, sanitize_input

# https://www.bioinformatics.org/sms/iupac.html
# dictionary of melting temperature methods with name as key and function as value
tm_methods = {"Nearest Neighbor": mt.Tm_NN, 
             "Empirical formulas based on GC content": mt.Tm_GC, 
             "Wallace, 'Rule of thumb'": mt.Tm_Wallace}
# dictionary of melting temperature models with name as key and function as value
tm_models = {"Nearest Neighbor": ['DNA_NN4', 'DNA_NN1', 'DNA_NN2', 'DNA_NN3', 'RNA_NN1', 'RNA_NN2', 'RNA_NN3', 'R_DNA_NN1'],
             "Empirical formulas based on GC content": [1,2,3,4,5,6,7,8], 
             "Wallace, 'Rule of thumb'": []}
# dictionary of Nearest Neighbour melting temperature models with name as key and function as value
NN_models = { "DNA_NN4": mt.DNA_NN4, "DNA_NN1": mt.DNA_NN1, "DNA_NN2": mt.DNA_NN2, "DNA_NN3": mt.DNA_NN3, "RNA_NN1": mt.RNA_NN1, "RNA_NN2": mt.RNA_NN2, "RNA_NN3": mt.RNA_NN3, "R_DNA_NN1": mt.R_DNA_NN1}
# dictionary of default values for the primer energy parameters
default_values = {"Na": 0, "K": 50, "Tris": 20, "Mg": 1.5, "dNTPs": 0.2, "Method": 7, "DMSO (%)": 0, "mt_method": list(tm_methods)[0], "mt_model": list(tm_models[list(tm_methods)[0]])[0], "Primer": 500}

def upload_setting_button():
    """Allow to upload setting"""
    st.session_state['upload_setting'] = True
    return

def primers_tab(seq):
    """ Calculate the melting temperature of the primers and check the dimer between the primers and the sequence"""
    
    ###
    # Melting Temperature for PCR settings
    ###
    mcol1, mcol2 = st.columns(2, vertical_alignment='center') 
    with mcol1:
        with st.popover("Melting temperature parameters"):
            # load settings if you want to upload custom settings
            load_settings = st.checkbox("Upload previous energy parameter:",
                                        help='''If you have saved a json file with the energy parameters, you can upload it and use the same settings.''')
            if load_settings:
                saved_setting = st.file_uploader("Upload previous energy parameter (optional):", on_change=upload_setting_button,
                                                type='json',
                                                help='''If you have saved a json file with the energy parameters, you can upload it and use the same settings.''')
                # if load the settings, upload the valuse from the json file
                if saved_setting and st.session_state['upload_setting']:
                    # To read file as bytes:
                    session_old = json.load(saved_setting)
                    for key, value in session_old.items():
                        st.session_state[key] = value
                    st.session_state['upload_setting'] = False
                    st.rerun()


            # initialize session state values for the energy model
            for key, val in default_values.items():
                if key not in st.session_state:
                    st.session_state[key] = val
        
            # select the energy model
            st.write('Energy Model:')
            col1, col2, col3 = st.columns([2, 2, 1])
            with col1:
                # select the melting temperature method
                mt_type = st.selectbox("Melting Temperature Method", list(tm_methods), help=(mt.__doc__), key='mt_method', label_visibility='collapsed')
                method = tm_methods[mt_type]
            with col2:
                # select the melting temperature model
                mt_model = st.selectbox("Energy Correction", list(tm_models[mt_type]), key='mt_model', label_visibility='collapsed')
            with col3:
                # add a button to show the help of the selected method
                help_model = st.button("Show info", key='energy_model_info')
            if help_model:
                st.help(method)
                
            # buffer correction parameters
            st.divider()
            if method != list(tm_methods)[2]:
                st.write('Buffer corrections (mM):')
                cols = st.columns(8)
                with cols[0]: na = st.number_input('Na', key="Na")
                with cols[1]: k = st.number_input('K', key='K')
                with cols[2]: tris = st.number_input('Tris', key='Tris')
                with cols[3]: mg = st.number_input('Mg', value=st.session_state["Mg"], key='Mg')
                with cols[4]: dntps = st.number_input('dNTPs', value=st.session_state["dNTPs"], key='dNTPs')
                with cols[5]: dmso = st.number_input('DMSO (%)', value=st.session_state["DMSO (%)"], key='DMSO (%)')
                with cols[6]: correction_met = st.selectbox('Method', [0, 1, 2, 3, 4, 6, 7], key='Method')
                with cols[7]: help_correction = st.button("Show info", key='energy_correction_info')

            if help_correction:
                # add a button to show the help of the selected method
                st.help(mt.salt_correction)
                st.help(mt.chem_correction)

            # create the function to calculate the TM
            if mt_type == list(tm_models)[0]:
                # primer correction
                st.divider()
                # take into account the primer concentration
                primer_conc = st.number_input('Primer conc. (nM)', key="Primer")
                calculate_mt = partial(method, nn_table=NN_models[mt_model], Na=na, K=k, Tris=tris, Mg=mg, dNTPs=dntps, dnac1=primer_conc, dnac2=0, saltcorr=correction_met)
            elif mt_type == list(tm_models)[0]:
                calculate_mt = partial(method, valueset=mt_model, Na=na, K=k, Tris=tris, Mg=mg, dNTPs=dntps, saltcorr=correction_met)
            else:
                calculate_mt = method
            
            def mt_correct(seq):
                # calculate the melting temperature and correct it with the primer correction
                return mt.chem_correction(calculate_mt(seq), DMSO=dmso)
                
            # save settings and allow download
            with open('energy_parameters.json', 'w') as settings:
                # cannot save session state as it is, I have to convert it to a dictionary
                session_dict = {key: st.session_state[key] for key in default_values}
                json.dump(session_dict, settings)
            with open("energy_parameters.json", "rb") as file:
                btn = st.download_button(
                    label="Download energy parameters",
                    data=file,
                    file_name="energy_parameters.json",
                    mime="application / json",
                    help='''Save the current settings (e.g. ions concentratio, Tm model), 
                            so you can easily reload them in you refresh the page!'''
                    )
    
    with mcol2:
        with st.popover('Add restriction sites'):
            restric_site_1 = st.text_input("Restriction site sequence (5'->3') before the promoter:")
            restric_site_2 = st.text_input("Restriction site sequence (5'->3') after the fragment:")

    seq = Seq(restric_site_1 + str(seq) + Seq(restric_site_2).reverse_complement())

    option_data = {'Manual': "bi bi-vector-pen",
                   'Automatic': "bi bi-person-workspace"}

    selected_operation = option_menu(None, 
                                    list(option_data.keys()),
                                    icons=list(option_data.values()),
                                    menu_icon="cast", 
                                    orientation="horizontal",
                                    styles=second_menu_style)
    
    if selected_operation == 'Automatic':
        st.write('Coming soon')
        st.stop()

    # show the settings for the two primers: choose the number of bases and show the gc content and melting temperature
    col1, col2 = st.columns(2)
    mts = [0, 0]
    # settings for the coding primer
    with col1:
        c_bases = st.slider('Forward primer length:', min_value=1, max_value=50, value = 21, key="coding_primer")
        c_primer = seq[:c_bases]
        copy_to_clipboard(c_primer, "Coding primer:")
        write_format_text(c_primer)
        mts[0] = round(mt_correct(c_primer),2)
        st.write(f"GC content: {round(gc_fraction(c_primer, ambiguous='ignore') * 100, 1)}%; Tm: {round(mts[0], 1)}°C")
    # settings for the non-coding primer
    with col2:
        nc_bases = st.slider('Reverese primer length:', min_value=1, max_value=50, value = 21, key="non_coding_primer")
        nc_primer = seq[-nc_bases:].reverse_complement()
        copy_to_clipboard(nc_primer, "Reverse primer:")
        write_format_text(nc_primer)
        mts[1] = round(mt_correct(nc_primer),2)
        st.write(f"GC content: {round(gc_fraction(nc_primer, ambiguous='ignore') * 100, 1)}%; Tm: {round(mts[1], 1)}°C")

    # show the primers preview, check the self-dimerization of the primer and check the dimer between the primers and the sequence
    with st.expander("Primers Preview"):
        st.markdown(f"""<div style="text-align: center;"><span style="color: #FF5733">{c_primer}</span>{seq[c_bases]}[...]{seq[-nc_bases-1:]}</div>""", unsafe_allow_html=True)
        st.markdown(f"""<div style="text-align: center;">{seq[:c_bases+1].complement()}[...]{seq[-nc_bases]}<span style="color: #FF5733">{nc_primer[::-1]}</span></div>""", unsafe_allow_html=True)
        st.divider()
        col1, col2 = st.columns(2)
        with col1:
            self_dimer1 = check_dimer(c_primer, c_primer)
            write_format_text('Self-Dimer coding primer\n' + self_dimer1)
        with col2:
            self_dimer2 = check_dimer(nc_primer, nc_primer)
            write_format_text('Self-Dimer non-coding primer\n' + self_dimer2)
        self_dimer12 = check_dimer(c_primer, nc_primer)
        write_format_text('Dimer between primers\n' + self_dimer12)
        
            

    st.divider()
    # add a warning if the melting temperature of the primers is too different
    if abs(mts[0] - mts[1]) >= 5:
        st.warning('The difference of Tm should be below 5°C', icon="⚠️")
                    
    # take into account the two main method to calculate the PCR annealing temperature: IDT method and Phusion Buffer method
    col1, col2, col3 = st.columns([2, 1, 3])
    ann_methods = ['IDT method [2]', 'Phusion method [3]']
    with col1:
        annealing_model = st.selectbox("Calculate annealing:", ann_methods)
    # the IDT method takes into account the GC melting temperature of the whole DNA sequence 
    if annealing_model == ann_methods[0]:
        t_anneal = round(0.3 * min(mts) + 0.7 * mt.chem_correction(mt.Tm_GC(seq, valueset=7, Na=na, K=k, Tris=tris, Mg=mg, dNTPs=dntps, saltcorr=correction_met), DMSO=dmso) - 14.9, 2)
    # the Phusion method takes into account the Melting temperature of the primers
    else:
        if len(c_primer) <= 20 or len(nc_primer) <= 20:
            t_anneal = min(mts)
        else:
            t_anneal = min(mts) + 3
        if t_anneal < 50:
            st.warning('Is suggested a temperature of annealing higher than 50°C.', icon="⚠️")
    # write the annealing temperature in big
    with col3:
        st.write('\n')
        st.subheader(f"T annealing: {t_anneal}")

# def auto_primer(sequence, target_temp):


if __name__ == "__main__":
    ### set the logo of the app
    load_logo()
    warnings.filterwarnings("ignore") # ignore warnings

    # create the tabs with the functions
    st.header('Prepare', help='Design primers for your DNA template or prepare the Origami for OxDNA simulation.')
    option_data = {'Primers': "bi bi-arrow-left-right",
                   'MD simulations': "bi bi-cpu"}

    selected_operation = option_menu(None, 
                                    list(option_data.keys()),
                                    icons=list(option_data.values()),
                                    menu_icon="cast", 
                                    orientation="horizontal",
                                    styles=main_menu_style)
    
    if selected_operation == 'MD simulations':
        st.write('Coming soon')
        st.stop()

    if 'streamlit_cwd' not in st.session_state:
        st.session_state.streamlit_cwd = str(Path(__file__).resolve().parent.parent)
    
    if "dna_template" not in st.session_state:
        st.session_state["dna_template"] = ''
    # st.header('Primer design')

    # take the input sequence and sanitize it
    seq = sanitize_input(st.text_input("Input sequence:", value = st.session_state["dna_template"]))
    # check the symbols in the sequence
    if set(seq) - symbols:
        st.warning('The sequence contains symbols not included in the [IUPAC alphabet](https://www.bioinformatics.org/sms/iupac.html).', icon="⚠️")
    if 'U' in seq:
        st.error('The DNA template contains U', icon=":material/personal_injury:")
    # create a porper biopython symbol
    seq = Seq(seq)

    
    if not seq:
        st.stop()
    elif len(seq) < 40:
        st.error('The sequence is too short', icon=":material/personal_injury:")
        st.stop()
    primers_tab(seq)

    # add bibliography
    reference(True)