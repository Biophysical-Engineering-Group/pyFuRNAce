from functools import partial
import streamlit as st
from colour import Color
import os
import warnings

### My modules
from utils import (check_import_pyfurnace, 
                   load_logo, 
                   save_origami, 
                   copy_to_clipboard)
check_import_pyfurnace()
from pyfurnace.generate import generate_road, fold_p

def format_text(text):
    return "```\n" + text + "\n```"

def forna_options(lenght):
    ### checkboxes for the options
    col1, col2, col3 = st.columns(3, vertical_alignment='bottom')
    with col1:
        st.session_state.zoomable = st.checkbox("Zoomable", 
                                                value=True, 
                help = """Enable zooming in and out the structure""")
    with col2:
        st.session_state.animation = st.checkbox("Interact", 
                                                 value=False, 
                help = """Enable interaction with the structure""")
    with col3:
        st.session_state.node_label = st.checkbox("Label", 
                                                  value=True, 
                help = """Show the nucleotide label""")
    
    ### second row of options
    col1, col2, col3 = st.columns(3, vertical_alignment='center')
    with col1:
        ### color scheme options
        scheme_options = ["sequence", 
                          "structure", 
                          "positions", 
                          "color range", 
                          "custom"]
        st.session_state.color_scheme = st.selectbox("Color scheme", 
                                                     scheme_options, 
                                                     index=1)

    with col2:
        st.session_state.height = st.slider("Frame height", 
                                            min_value=10, 
                                            max_value=800, 
                                            value=300, 
                help = """Set the height of the frame""")

    with col3:
        st.session_state.label_interval = st.slider("Label every", 
                                                    min_value=0, 
                                                    max_value=max(lenght, 10), 
                                                    value=min(lenght, 10), 
                help = """Show the nucleotide number every n nucleotides""")

    st.session_state.color_text = ''
    st.session_state.colors = {}

    ### color scheme based on the sequence
    if st.session_state.color_scheme == "color range":
        st.session_state.color_scheme = "custom"

        with col2: # start color
            first = st.color_picker("Start color", "#ff0000")

        with col3: # end color
            last = st.color_picker("End color", "#00ff00")

        # create the colors range
        first = Color(first)
        last = Color(last)

        color_range = [first] + list(first.range_to(last, lenght))

        # save the colors in the session state and create the colors string
        for i, c in enumerate(color_range):
            st.session_state.colors[i] = c.hex
            st.session_state.color_text += str(i) + ":" + c.hex + " "

    ### custom colors for each nucleotide
    elif st.session_state.color_scheme == "custom":

        with col2: # select nucleotide index
            index = st.number_input("Select nucleotide index", 1, lenght, 1)

        with col3: # select color
            color = st.color_picker("Select a color", "#ffffff")

        # save the color in the session state
        st.session_state.colors[index] = color

        # create the colors string
        for i, c in st.session_state.colors.items():
            st.session_state.color_text += str(i) + ":" + c + " "

    ### display the custom colors
    else:
        st.session_state.color_text = None

def generate_sequence():
    if not structure:
        st.error('No structure input given')
        return
    
    if not sequence_constraint:
        st.error('No sequence constraint input given')
        return
        
    # Initialize the UI elements
    progress_bar = st.empty()
    output_status = st.empty()
    
    # Initialize the progress bar
    progress_bar.progress(0, "Initializing...")

    def callback_forna(structure, sequence, step, stage_name, n_stage):
        # Update the progress bar
        progress_bar.progress(n_stage, text = stage_name)
        output_status.markdown(format_text(f"Structure:\n{structure}"
                                            f"\nSequence:\n{sequence}"
                                            f"\nSpool:\n{step}"))

    ori_txt = '\n\n'.join(st.session_state.code)

    if st.button("Stop optimization"):
        return
    
    opti_sequence, zip_out = generate_road(structure, 
                                           sequence_constraint, 
                                           pseudoknot_info, 
                                           name=filename,
                                           callback=callback_forna,
                                           zip_directory=True,
                                           origami_code=ori_txt,
                                           )


    # Clear the progress bar
    progress_bar.progress(1.0, text = f"Generation finished")
    output_status.empty()


    if not opti_sequence:
        st.warning("Optimization failed." 
                   "Please check the input parameters for Revoler.")
        with open(zip_out, "rb") as fp:
            st.download_button("Download failed optimization files",
                                data=fp,
                                file_name = f"{filename}.zip",
                                mime = "application/zip",
                                on_click='ignore',
                                )

        return
    
    # if there was already a zip file, remove it
    if 'zip_path' in st.session_state:
        if os.path.exists(st.session_state.zip_path):
            os.remove(st.session_state.zip_path)
    st.session_state.rna_origami_seq = opti_sequence
    st.session_state.rna_origami_folds = fold_p(opti_sequence)
    # save the new zip file
    st.session_state.zip_path = zip_out

    if ('origami' in st.session_state 
            and st.session_state.origami
            and len(st.session_state.origami.sequence) == len(st.session_state.rna_origami_seq)):
        # try:
        st.session_state.origami.sequence = st.session_state.rna_origami_seq
        if 'code' in st.session_state:
            code_text = 'origami.sequence = "' + st.session_state.rna_origami_seq + '"'
            st.session_state.code.append(code_text)
        st.success("Sequence loaded to the origami design!")
        # except Exception as e:
        #     st.error(f"Error while loading the sequence into the origami design: {e}")


if __name__ == "__main__":
    # somehow st components cause `Thread 'MainThread': missing ScriptRunContext!` 
    # warning when using multiprocessing
    from st_forna_component import forna_component

    load_logo()

    ### ignore warnings
    warnings.filterwarnings("ignore")

    ### initialize the session state

    ### Forna options
    if 'zoomable' not in st.session_state:
        st.session_state.zoomable = True
    if 'animation' not in st.session_state:
        st.session_state.animation = False
    if 'editable' not in st.session_state:
        st.session_state.editable = False
    if 'labels' not in st.session_state:
        st.session_state.node_label = True
    if 'height' not in st.session_state:
        st.session_state.height = 300
    if 'label_interval' not in st.session_state:
        st.session_state.label_interval = 10
    if 'color_scheme' not in st.session_state:
        st.session_state.color_scheme = "structure"
    if 'color_text' not in st.session_state:
        st.session_state.color_text = ''

    ### RNA origami options
    if 'generate_structure' not in st.session_state:
        st.session_state.generate_structure = ""
    if 'generate_sequence' not in st.session_state:
        st.session_state.generate_sequence = ""
    if 'generate_pseudoknots' not in st.session_state:
        st.session_state.generate_pseudoknots = ""
    if 'rna_origami_seq' not in st.session_state:
        st.session_state.rna_origami_seq = ""
    if 'rna_origami_folds' not in st.session_state:
        st.session_state.rna_origami_folds = ()
    
    st.header('Generate', help='Generate the RNA sequence that matches the desired dot-bracket notation (and sequence constraints) for the nanostructure.')
    structure = st.text_input("RNA strucutre (dot-bracket notation)", value=st.session_state.generate_structure)
    sequence_constraint = st.text_input("Sequence constraints", value=st.session_state.generate_sequence)
    pseudoknot_info = st.text_input("Pseudoknot constraints (semicolon-separated)", value=st.session_state.generate_pseudoknots)
    
    if not sequence_constraint:
        sequence_constraint = "N" * len(structure.replace("&", ""))

    # check if the dot-bracket notation contains multiple strands   
    if "&" in structure:
        st.warning("Experimental: the dot-bracket notation contains multiple strands")

    with st.columns([1.3, 3])[0]:
        with st.popover("RNA display options", use_container_width=True):
            forna_options(len(structure))

    partial_forna = partial(forna_component,
                           height = st.session_state.height,
                           animation = st.session_state.animation,
                           zoomable = st.session_state.zoomable,
                           label_interval = st.session_state.label_interval,
                           node_label = st.session_state.node_label,
                           editable = st.session_state.editable,
                           color_scheme = st.session_state.color_scheme,
                           colors = st.session_state.color_text)
    

    # # Initialize the FORNA link for the target structure
    if structure:
    #     st.markdown("#### Target Structure")
    #     st.markdown(format_text(structure))
        edited = partial_forna(structure = structure, 
                                sequence = sequence_constraint,
                                key='target_forna')

    if sequence_constraint and sequence_constraint[0] != 'G':
        col1, col2, col3 = st.columns(3, vertical_alignment='bottom')
        with col1:
            st.warning("The RNA sequence doens't start with G.")
        with col2:
            new_start = st.text_input("Start the sequence with", value='GGGA',
                                      help='The transcription of the RNA sequence often '
                                      'often requires at least one G at the start of the sequence.')
        with col3:
            if st.button('Apply the sequence start'):

                seq_list = list(sequence_constraint)
                # pair_map = pf.dot_bracket_to_pair_map(structure) # unnecessary

                for i in range(len(new_start)):
                    seq_list[i] = new_start[i]
                    ### NOT NECESSARY
                    # paired = pair_map[i]
                    # if paired is not None:
                    #     paired_nucl = new_start[i].translate(pf.nucl_to_pair)
                    #     seq_list[paired] = paired_nucl

                st.session_state.generate_sequence = "".join(seq_list)
                st.rerun()
    
    col1, col2 = st.columns(2, vertical_alignment='bottom')
    with col1:
        filename =  st.text_input('Name of RNA origami', value='Origami')
    with col2: 
        generate = False
        if st.button("Generate RNA sequence"):
            generate = True
    if generate:
        generate_sequence()


    if not st.session_state.rna_origami_seq:
        st.stop()
        
    st.divider()
    sequence = st.session_state.rna_origami_seq
    folds = st.session_state.rna_origami_folds

    diversity = round(folds[6], 1)
    if diversity < 30:
        diversity_text = f":green[low {diversity}]"
    elif diversity < 50:
        diversity_text = f":orange[medium {diversity}]"
    else:
        diversity_text = f":red[high {diversity}]"

    cols = st.columns(4, vertical_alignment='bottom')
    st.markdown("### Last Optimized sequence "
                f"(ensemble diversity: {diversity_text})", 
                help='The ensemble diversity is the average distance between the '
                'structures in the ensemble (the set of all the possible structures '
                'that can be formed by the sequence). A lower value means that the '
                'structures are more similar to each other (and therefore more '
                'similar to the Minimum Free Energy structure). A higher value '
                'means that the structures are more diverse and there are less '
                'chances to obtain the minimum free energy structure. '
                )
    

    col1, col2 = st.columns(2)
    with col1:
        with st.columns(3)[1]:
            st.markdown("#### MFE Structure", 
                        help='The Minimum Free Energy (MFE) structure is the structure with '
                        'the lowest free energy. It is the most stable structure that can be '
                        'formed by the sequence. '
                        )
            
        subcol1, subcol2 = st.columns(2)
        with subcol1:
            st.markdown(f'Energy: {round(folds[1], 2)} Kcal/mol')
        with subcol2:
            st.markdown(f'Frequency in the ensemble: {round(folds[2] * 100, 4)} %',
                        help = 'The MFE frequency in the ensemble is the probability of '
                        'obtaining the MFE structure among all the possible structures '
                        'that can be formed by the sequence.'
                        )
        partial_forna(structure = folds[0],
                        sequence = sequence,
                        key='struct_mfe')
    with col2:
        with st.columns(3)[1]:
            st.markdown("#### Centroid",
                        help='The centroid structure is the structure that is the most similar to '
                        'the average structure of the ensemble. It is the closest structure to '
                        'represent the average of all the possible structures that can be '
                        'formed by the sequence. '
                        )
        subcol1, subcol2 = st.columns(2)
        with subcol1:
            st.markdown(f'Energy: {round(folds[4], 2)} Kcal/mol')
        with subcol2:
            st.markdown(f'Frequency in the ensemble: {round(folds[5] * 100, 4)} %',
                        help = 'The centroid frequency in the ensemble is the probability of '
                        'obtaining the centroid structure among all the possible structures '
                        'that can be formed by the sequence.'
                        )
        partial_forna(structure = folds[3],
                        sequence = sequence,
                        key='struct_centroid')
        
    st.divider()

    st.markdown(format_text(st.session_state.rna_origami_seq))
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.page_link("pages/3_Convert.py", 
                        label=":orange[Convert RNA to DNA]", 
                        icon=":material/genetics:")
    with col2:
        copy_to_clipboard(folds[0], 'Structure')
    with col3:
        copy_to_clipboard(sequence, 'Sequence')
    with col4:
        with open(st.session_state.zip_path, "rb") as fp:
            st.download_button("Download optimization files",
                                data=fp,
                                file_name = f"{filename}.zip",
                                mime = "application/zip",
                                on_click='ignore',
                                )

    if ('origami' in st.session_state 
            and st.session_state.origami
            and str(st.session_state.origami.sequence) == st.session_state.rna_origami_seq):

        st.session_state.prepare_ind = 1
        st.page_link("pages/4_Prepare.py",
                        label=":orange[Prepare MD simulations]", 
                        icon=":material/sync_alt:")

        save_origami(filename)
                
