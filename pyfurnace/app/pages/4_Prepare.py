import streamlit as st
import os
from functools import partial
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import gc_fraction
from streamlit_option_menu import option_menu
import json
import warnings
### import the template functions
from utils import load_logo, main_menu_style, second_menu_style, copy_to_clipboard
from utils.template_functions import symbols, write_format_text, check_dimer, reference, sanitize_input
from pyfurnace.prepare import oxdna_simulations
from utils.design_functions import origami_build_view

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

def calculate_annealing(seq, mts, c_primer, nc_primer, tm_kwargs):

    col1, col2 = st.columns(2, gap='large')
    with col1:
        subcol1, subcol2 = st.columns([11,1], vertical_alignment='center')
        with subcol1:
            write_format_text(c_primer)
        with subcol2:
            copy_to_clipboard(c_primer, "")

        subcol1, subcol2, subcol3 = st.columns(3)
        with subcol1:
            st.markdown(f"GC content: {round(gc_fraction(c_primer, ambiguous='ignore') * 100, 1)}%")
        with subcol2:
            st.markdown(f"**:green[Tm: {round(mts[0], 1)}°C]**")
        with subcol3:
            st.markdown(f"Length: {len(c_primer)}")

    with col2:
        
        subcol1, subcol2 = st.columns([11,1], vertical_alignment='center')
        with subcol1:
            write_format_text(nc_primer)
        with subcol2:
            copy_to_clipboard(nc_primer, "")

        subcol1, subcol2, subcol3 = st.columns(3)
        with subcol1:
            st.markdown(f"GC content: {round(gc_fraction(nc_primer, ambiguous='ignore') * 100, 1)}%")
        with subcol2:
            st.markdown(f"**:green[Tm: {round(mts[1], 1)}°C]**")
        with subcol3:
            st.markdown(f"Length: {len(nc_primer)}")

    # show the primers preview, check the self-dimerization of the primer and check the dimer between the primers and the sequence
    c_bases = len(c_primer)
    nc_bases = len(nc_primer)
    with st.expander("Dimerization preview"):
        st.markdown(f"""<div style="text-align: center;"><span style="color: #FF5733">{c_primer}</span>{seq[c_bases]}[...]{seq[-nc_bases-1:]}</div>""", unsafe_allow_html=True)
        st.markdown(f"""<div style="text-align: center;">{seq[:c_bases+1].complement()}[...]{seq[-nc_bases]}<span style="color: #FF5733">{nc_primer[::-1]}</span></div>""", unsafe_allow_html=True)
        st.divider()
        col1, col2 = st.columns(2, gap='large')
        with col1:
            self_dimer1 = check_dimer(c_primer, c_primer, basepair={'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'})
            write_format_text('Self-Dimer forward primer\n' + self_dimer1)
        with col2:
            self_dimer2 = check_dimer(nc_primer, nc_primer, basepair={'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'})
            write_format_text('Self-Dimer reverse primer\n' + self_dimer2)
        self_dimer12 = check_dimer(c_primer, nc_primer, basepair={'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'})
        write_format_text('Dimer between primers\n' + self_dimer12)

    # add a warning if the melting temperature of the primers is too different
    anneal_color = 'green'
    if abs(mts[0] - mts[1]) >= 5:
        st.warning('The difference of Tm should be below 5°C', icon="⚠️")
        anneal_color = 'red'
                    
    # take into account the two main method to calculate the PCR annealing temperature: IDT method and Phusion Buffer method
    col1, col2 = st.columns(2, gap='large', vertical_alignment='bottom')
    ann_methods = ['IDT method [2]', 'Phusion method [3]']
    with col1:
        annealing_model = st.selectbox("Calculate annealing:", ann_methods)
        
    # the IDT method takes into account the GC melting temperature of the whole DNA sequence 
    if annealing_model == ann_methods[0]:
        method_kwargs = tm_kwargs.copy()
        method_kwargs.pop('dNTPs')
        method_kwargs.pop('DMSO')
        t_anneal = round(0.3 * min(mts) + 0.7 * mt.chem_correction(mt.Tm_GC(seq, valueset=7, **method_kwargs), DMSO=tm_kwargs['DMSO']) - 14.9, 2)
    # the Phusion method takes into account the Melting temperature of the primers
    else:
        if len(c_primer) <= 20 or len(nc_primer) <= 20:
            t_anneal = min(mts)
        else:
            t_anneal = min(mts) + 3
        if t_anneal < 50:
            st.warning('Is suggested a temperature of annealing higher than 50°C.', icon="⚠️")
            anneal_color = 'red'
    # write the annealing temperature in big
    with col2:
        st.markdown(f"### Anneal at: :{anneal_color}[{t_anneal}°C]")

def primers_tab(seq, mt_correct, tm_kwargs):
    """ Calculate the melting temperature of the primers and check the dimer between the primers and the sequence"""

    # show the settings for the two primers: choose the number of bases and show the gc content and melting temperature
    st.write('\n')
    col1, col2 = st.columns(2, gap='large')
    mts = [0, 0]
    # settings for the coding primer
    with col1:
        with st.columns(5)[2]:
            st.markdown('###### Forward')

        c_bases = st.slider('Primer length:', min_value=5, max_value=50, value = 21, 
                            key="coding_primer")
        c_primer = seq[:c_bases]

        mts[0] = round(mt_correct(c_primer),2)

    # settings for the non-coding primer
    with col2:
        with st.columns(5)[2]:
            st.markdown('###### Reverse')

        nc_bases = st.slider('Primer length:', min_value=5, max_value=50, value = 21, 
                             key="non_coding_primer")
        nc_primer = seq[-nc_bases:].reverse_complement()

        mts[1] = round(mt_correct(nc_primer),2)

    calculate_annealing(seq, mts, c_primer, nc_primer, tm_kwargs)

def auto_primer(seq, mt_correct, tm_kwargs):
    """ Automatically design the primers for the sequence"""
    # check if the sequence is a DNA sequence
    target_temp = st.number_input('Target melting temperature (°C)', value=65, min_value=0, max_value=100, step=1)

    if st.button("Design primers", key='auto_primer'):
        status = st.empty()
        status.info("Designing primers...", icon=":material/precision_manufacturing:")

        final_mts = []
        primers = []

        # check the melting temperature of the primers
        for direct in (1, -1):
            prim_length = 10
            primers_info = list()

            while True:

                to_prime = seq
                if direct == -1:
                    to_prime = seq.reverse_complement()
                primer = to_prime[:prim_length]
                tm = round(mt_correct(primer), 2)

                if tm <= (target_temp - 2.5):
                    prim_length += 1
                    continue
                elif tm >= (target_temp + 2.5):
                    break

                score = 0
                if primer[-1] in 'GC':
                    score += 1
                if primer[-2] in 'GC':
                    score += 1
                if 18 < prim_length < 30:
                    score += 1
                elif prim_length < 18:
                    score -= 17 - prim_length
                if 40 < round(gc_fraction(primer, ambiguous='ignore') * 100, 1) < 60:
                    score += 1
                score -= abs(tm - target_temp) / 2

                dimers = check_dimer(primer, 
                                     primer, 
                                     dict_format=True,
                                     basepair={'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'})
                max_dimer = max(dimers.keys())
                score -= max_dimer / (len(primer) // 2)


                primers_info.append((score, tm, primer))

                prim_length += 1
                if prim_length > len(seq) - 1 or not tm:
                    st.error('No primer found for the sequence', icon=":material/personal_injury:")
                    break

            # sort the primers by score
            primers_info.sort(reverse=True)
            # st.write(primers_info)
            final_mts.append(primers_info[0][1])
            primers.append(primers_info[0][2])

            status.success(f"Primer designed!", icon=":material/precision_manufacturing:")

        col1, col2 = st.columns(2, gap='large')
        # settings for the coding primer
        with col1:
            with st.columns(5)[2]:
                st.markdown('###### Forward')
        with col2:
            with st.columns(5)[2]:
                st.markdown('###### Reverse')

        # check the dimerization of the primers
        calculate_annealing(seq, final_mts, primers[0], primers[1], tm_kwargs)
        
def primers_setup():
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
    
    # create a proper biopython sequence
    seq = Seq(seq)

    
    if not seq:
        st.stop()
    elif len(seq) < 40:
        st.error('The sequence is too short', icon=":material/personal_injury:")
        st.stop()
    
    ###
    # Melting Temperature for PCR settings
    ###

    mcol1, mcol2 = st.columns(2, vertical_alignment='center') 
    with mcol1:
        with st.popover("**Melting temperature parameters**"):
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
            st.write('Energy model:')
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

                tm_kwargs = {'Na': na, 
                             'K': k, 
                             'Tris': tris, 
                             'Mg': mg, 
                             'saltcorr': correction_met,
                             'dNTPs': dntps, 
                             'DMSO': dmso}

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
                method_kwargs = tm_kwargs.copy()
                method_kwargs.pop('dNTPs')
                method_kwargs.pop('DMSO')
                calculate_mt = partial(method, nn_table=NN_models[mt_model], dnac1=primer_conc, dnac2=0, **method_kwargs)
            
            elif mt_type == list(tm_models)[1]:
                method_kwargs = tm_kwargs.copy()
                method_kwargs.pop('DMSO')
                calculate_mt = partial(method, valueset=mt_model, **method_kwargs)
            else:
                calculate_mt = method
            
            def mt_correct(seq):
                # calculate the melting temperature and correct it with the primer correction
                try:
                    return mt.chem_correction(calculate_mt(seq), DMSO=dmso)
                except ValueError:
                    st.error('Error in the melting temperature calculation. Try to switch page and come back.',
                             icon=":material/personal_injury:")
                return 0
                
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
    
    if selected_operation == 'Manual':
        primers_tab(seq, mt_correct, tm_kwargs)
    elif selected_operation == 'Automatic':
        auto_primer(seq, mt_correct, tm_kwargs)

    # add bibliography
    reference(True)

def md_setup():
    if (not st.session_state.get('origami') or
        not all(n in 'AUCG' for n in st.session_state.origami.sequence)):

        st.error("The origami blueprint is empty or it doesn't contain a sequence",
                 icon=":material/personal_injury:")
        
        col1, col2 = st.columns(2)
        with col1:
            st.page_link("pages/1_Design.py", 
                         label="**:orange[Have a look to the blueprint]**", 
                         icon=":material/draw:")
        
        with col2:
            st.page_link("pages/2_Generate.py", 
                         label="**:orange[Generate the sequence]**", 
                         icon=":material/network_node:")
        st.stop()
    
    with st.expander('**Origami preview:**', expanded=True):
        origami_build_view('Origami split view')
    md_input()

@st.fragment
def md_input():
    st.divider()
    col1, col2 = st.columns(2)
    with col1:
        oxdna_dir = st.text_input('OxDNA directory (e.g. "~/Documents/software/oxDNA")')
    with col2:
        temp = st.number_input('Temperature (°C)', value=37, min_value=0, max_value=100, step=1)

    help_mc_rel = ("Use Monte Carlo relaxation to relax the origami before the MD "
                  "simulation, using external forces to keep the basepairing.")
    help_md_rel = ("Use Molecular Dynamics to relax the structure, using forces to"
                   "maintain the basepairing")
    help_md_equ = ("An extra Molecular Dynamics simulation to make sure the Origami"
                   "is fully relaxed and not bias the starting configuration of the " 
                   "production. **External forces for pseudoknots only** are used "
                   "in this simulation, to make sure the pseudoknots are paired "
                   "at the start of the production")
    help_md_run = ("The final Molecular Dynamics simulation, to simulate the "
                   "structure without any external force.")
    
    num_input = st.toggle("Use numerical input", key='num_input')

    col1, col2, col3, col4 = st.columns(4)
    with col1:
        if num_input:
            step_mc_rel = st.number_input('MC relaxation steps', 
                                          min_value=0,
                                          value=int(5e3),
                                          help=help_mc_rel,
                                          )
        else:
            step_mc_rel = st.slider('MC relaxation steps', 
                                    min_value=0,
                                    value=int(5e3),
                                    max_value=int(10e3),
                                    step=int(1e3),
                                    help=help_mc_rel,
                                    format="%0.0e"
                                    )
    with col2:
        if num_input:
            step_md_rel = st.number_input('MD relaxation steps', 
                                          min_value=0,
                                          value=int(1e7),
                                          help=help_md_rel,
                                          )
        else:
            step_md_rel = st.slider('MD relaxation steps', 
                                    min_value=0,
                                    value=int(1e7),
                                    max_value=int(1e8),
                                    step=int(1e6),
                                    help=help_md_rel,
                                    format="%0.0e"
                                    )
    with col3:
        if num_input:
            step_md_equ = st.number_input('MD equilibration steps', 
                                          min_value=0,
                                          value=int(1e8),
                                          help=help_md_equ,
                                          )
        else:
            step_md_equ = st.slider('MD equilibration steps', 
                                    min_value=0,
                                    value=int(1e8),
                                    max_value=int(1e9),
                                    step=int(1e7),
                                    help=help_md_equ,
                                    format="%0.0e"
                                    )
    with col4:
        if num_input:
            step_md_run = st.number_input('MD production steps', 
                                          min_value=0,
                                          value=int(1e9),
                                          help=help_md_run,
                                          )
        else:
            step_md_run = st.slider('MD production steps', 
                                    min_value=0,
                                    value=int(1e9),
                                    max_value=int(1e10),
                                    step=int(1e8),
                                    help=help_md_run,
                                    format="%0.0e"
                                    )
            
    if not oxdna_dir or not temp:
        st.stop()

    zip_path = oxdna_simulations(origami=st.session_state.origami,
                                oxdna_directory=oxdna_dir,
                                temperature=temp,
                                mc_relax_steps=step_mc_rel,
                                md_relax_steps=step_md_rel,
                                md_equil_steps=step_md_equ,
                                md_prod_steps=step_md_run
                                )
    # add a button to download the zip file
    st.divider()
    col1, col2 = st.columns(2, vertical_alignment='bottom')
    with col1:
        name = st.text_input('Name of the simulation:', value='origami_simulation')
    with col2:
        with open(zip_path, "rb") as file:
            st.download_button(
                label="Download MD simulation files",
                data=file,
                file_name=f"{name}.zip",
                mime="application/zip",
                help='''Download the MD simulation files. 
                        The zip file contains the input files for the MD simulation.'''
                )
        
    # delete the zip file
    if os.path.exists(st.session_state.zip_path):
        os.remove(zip_path)
    return

if __name__ == "__main__":
    ### set the logo of the app
    load_logo()
    warnings.filterwarnings("ignore") # ignore warnings

    if "prepare_ind" not in st.session_state:
        st.session_state.prepare_ind = 0

    # create the tabs with the functions
    st.header('Prepare', help='Design primers for your DNA template or prepare the Origami for OxDNA simulation.')
    option_data = {'Primers': "bi bi-arrow-left-right",
                   'MD simulations': "bi bi-cpu"}

    selected_operation = option_menu(None, 
                                    list(option_data.keys()),
                                    icons=list(option_data.values()),
                                    menu_icon="cast", 
                                    orientation="horizontal",
                                    default_index=st.session_state.prepare_ind,
                                    styles=main_menu_style)
    
    if selected_operation == 'MD simulations':
        if st.session_state.prepare_ind == 0:
            st.session_state.prepare_ind = 1
            st.rerun()
        md_setup()
    elif selected_operation == 'Primers':
        if st.session_state.prepare_ind == 1:
            st.session_state.prepare_ind = 0
            st.rerun()
        primers_setup()
    