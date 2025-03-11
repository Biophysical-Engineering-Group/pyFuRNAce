import os
import zipfile
import streamlit as st
from colour import Color
import warnings
import tempfile
################################OWN####LIBRARIES#######################################
from utils import check_import_pyfurnace, load_logo, save_origami
check_import_pyfurnace()

def format_text(text):
    return "```\n" + text + "\n```"

def forna_options(lenght):
    ### checkboxes for the options
    col1, col2, col3 = st.columns(3, vertical_alignment='bottom')
    with col1:
        st.session_state.zoomable = st.checkbox("Zoomable", value=True, help = """Enable zooming in and out the structure""")
    with col2:
        st.session_state.animation = st.checkbox("Interact", value=False, help = """Enable interaction with the structure""")
    with col3:
        st.session_state.node_label = st.checkbox("Label", value=True, help = """Show the nucleotide label""")
    
    ### second row of options
    col1, col2, col3 = st.columns(3, vertical_alignment='center')
    with col1:
        ### color scheme options
        scheme_options = ["sequence", "structure", "positions", "color range", "custom"]
        st.session_state.color_scheme = st.selectbox("Color scheme", scheme_options, index=1)
    with col2:
        st.session_state.height = st.slider("Frame height", 10, 800, 300, help = """Set the height of the frame""")
    with col3:
        st.session_state.label_interval = st.slider("Label every", 0,  max(lenght, 10), min(lenght, 10), help = """Show the nucleotide number every n nucleotides""")

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

    st.header('Generate', help='Generate the RNA sequence that matches the desired dot-bracket notation (and sequence constraints) for the nanostructure.')
    structure = st.text_input("RNA Strucutre (dot-bracket notation)", value=st.session_state.generate_structure)
    sequence_constraint = st.text_input("Sequence Constraints", value=st.session_state.generate_sequence)
    pseudoknot_info = st.text_input("Pseudoknot Constraints (semicolon-separated)", value=st.session_state.generate_pseudoknots)
    
    if not sequence_constraint:
        sequence_constraint = "N" * len(structure.replace("&", ""))

    # check if the dot-bracket notation contains multiple strands   
    if "&" in structure:
        st.warning("Experimental: the dot-bracket notation contains multiple strands")

    with st.columns([1.3, 3])[0]:
        with st.popover("RNA display options", use_container_width=True):
            forna_options(len(structure))

    # # Initialize the FORNA link for the target structure
    if structure:
    #     st.write("#### Target Structure")
    #     st.write(format_text(structure))
        edited = forna_component(structure = structure, 
                                sequence = sequence_constraint,
                                height = st.session_state.height,
                                animation = st.session_state.animation,
                                zoomable = st.session_state.zoomable,
                                label_interval = st.session_state.label_interval,
                                node_label = st.session_state.node_label,
                                editable = st.session_state.editable,
                                color_scheme = st.session_state.color_scheme,
                                colors = st.session_state.color_text,
                                key='target_forna')
        
    generate = True
    col1, col2 = st.columns(2)
    with col1:
        filename =  st.text_input('Name of RNA origami', value='Origami')
    with col2: 
        st.write("\n"); st.write("\n"); 
        if st.button("Generate RNA Sequence"):
            generate = False
    if not generate:
        # Initialize the UI elements
        progress_bar = st.empty()
        best_info = st.empty()
        stats = st.empty()
        # Initialize the progress bar
        progress_bar.progress(0, "Initializing...")
        # Generate the RNA origami
        tot_strands = structure.count("&") + 1

        
        with tempfile.TemporaryDirectory() as tmpdirname:
            try:
                pass
                ### TODO
            except Exception as e:
                st.error(f"Error: {e}")

    # TO FIX
    #     # Clear the progress bar
    #     progress_bar.progress(1.0, text = f"Generation Completed")
    #     st.session_state.rna_origami_seq = full_seq
    #     st.session_state.rna_origami_struct = full_struct
    #     # Clear progess bar
    #     update_ui(None, None, None, None)
    #     # progress_bar.empty()

    #     # Display the file to download
    
    # if st.session_state.rna_origami_seq:
    #     st.write("### Last Optimized sequence:")
    
    #     st.write(format_text(st.session_state.rna_origami_seq))
    #     col1, col2, col3 = st.columns(3)
    #     with col1:
    #         st.page_link("pages/3_Template.py", 
    #                      label=":orange[Prepare the DNA template]", 
    #                      icon=":material/genetics:")
    #     with col2:
    #         st_copy_to_clipboard(st.session_state.rna_origami_seq, 
    #                              before_copy_label='Copy Sequence📋', 
    #                              show_text=False)
    #     with col3:
    #         st.download_button("Download Results",
    #                 st.session_state.zip_results,
    #                 f"{filename}.zip")
            
    #     forna_component(structure = st.session_state.rna_origami_struct, 
    #                     sequence = st.session_state.rna_origami_seq,
    #                     height = st.session_state.height,
    #                     animation = st.session_state.animation,
    #                     zoomable = st.session_state.zoomable,
    #                     label_interval = st.session_state.label_interval,
    #                     node_label = st.session_state.node_label,
    #                     editable = st.session_state.editable,
    #                     color_scheme = st.session_state.color_scheme,
    #                     colors = st.session_state.color_text,
    #                     key='optimized_forna')

        if 'origami' in st.session_state and st.session_state.origami:
            if st.button("Load the sequence to the Origami"):
                try:
                    st.session_state.origami.sequence = st.session_state.rna_origami_seq
                    if 'code' in st.session_state:
                        st.session_state.code.append('origami.sequence = "' + st.session_state.rna_origami_seq + '"')
                    st.success("Sequence Loaded to the Origami")
                except Exception as e:
                    st.error(f"Error: {e}")
            if all(n in 'AUCG' for n in st.session_state.rna_origami_seq):
                save_origami(filename)
                    
