import inspect
import streamlit as st
from streamlit import session_state as st_state

### pyFuRNAce modules
from utils import load_logo, check_import_pyfurnace
from utils.template_functions import sanitize_input
check_import_pyfurnace()
import pyfurnace as pf
import pyfurnace.design.utils.origami_lib as origami_lib


def quickstart():
    with st.popover(":red[New to pyFuRNAce? Start here!]",
                    use_container_width=True,):
        st.write("Do you have a sequence/structure (dot-bracket)?")
        start = st.radio("Select an option:", 
                         index=None,
                         options=["Yes", 
                                  "No"],
                         horizontal=True,
                         label_visibility="collapsed",
                         key="design_choice")
        
        if not start:
            return
        
        elif start == "Yes":
            st.write("Great! You can start by adding either your sequence, "
                     "your dot-bracket structure or both.")
            col1, col2 = st.columns(2)
            with col1:
                sequence = st.text_input("Sequence:", 
                                        key="Sequence", 
                            help="To generate an RNA origami from a sequence," 
                                 " we fold it with viennaRNA and then convert " 
                                 "the dot bracket notation to an Origami, "
                                 "as an collection of stems and bulges.") 

            with col2:
                structure = st.text_input("Dot-bracket structure:", 
                                          key="Dot-bracket_structure", 
                            help="Add a dot-bracket structure to generate "
                                " a blueprint. If the secondary structure "
                                " of known aptamers is detected, the aptamer "
                                " motif will be used in the Origami. ")
                
            structure = sanitize_input(structure)
            sequence = sanitize_input(sequence) 

            if not (structure or sequence):
                return
            try:
                pf.Origami.from_structure(structure=structure, 
                                          sequence=sequence)
            except ValueError as e:
                st.error(f"Error: {e}", icon=":material/personal_injury:")
            else:
                st.success("Origami created successfully!", 
                           icon=":material/check_circle:")
                
                st_state.code = ['import pyfurnace as pf\n'
                                         f'origami = pf.Origami.from_structure('
                                         f'structure="{structure}", '
                                         f'sequence="{sequence}")']
                st_state.origami = pf.Origami.from_structure(structure=structure,
                                                             sequence=sequence)
                
                st.page_link("pages/1_Design.py", 
                            label=":orange[Switch to the Design page]", 
                            icon=":material/draw:")
                
        elif start == "No":
            st.write("No problem! You can design from scratch, "
                     "or loading a template.")

            templ_funcs = {name.replace('template_', ''): func 
                                for name, func in vars(origami_lib).items()
                                    if callable(func) and name.startswith("template_")}

            template_docs = 'Available templates:\n\n'
            for name, func in templ_funcs.items():
                template_docs += f'**--> {name}** \n\n ' 
                docs = inspect.getdoc(func)
                template_docs += f'{docs.split("Returns")[0].strip()}\n\n'

            template = st.selectbox("Select a template:",
                                    options=list(templ_funcs.keys()),
                                    index=None,
                                    help=template_docs)
            if template:
                func_code = inspect.getsource(templ_funcs[template])
                # remove one layer of indentation
                only_code = func_code.split('"""')[2].split('return')[0]
                func_code = '\n'.join([line[4:] for line in only_code.splitlines()])

                st_state.origami = templ_funcs[template]()
                st_state.code = [f'{func_code.strip()}\n']
            
                st.success(f"Template {template} loaded successfully!", 
                           icon=":material/check_circle:")
                
            st.page_link("pages/1_Design.py", 
                        label=":orange[Switch to the Design page]", 
                        icon=":material/draw:")



if __name__ == '__main__':
    load_logo() 
    check_import_pyfurnace()

    st.write("# Hello and Welcome to pyFuRNAce!")

    st.write('Design and generate RNA nanostructures in few simple steps.')

    st.page_link("pages/1_Design.py", 
                 label=":orange[Design:]", 
                 icon=":material/draw:")

    st.markdown("- Design your RNA nanostructure and download it as "
                "textfile/python script.")

    st.page_link("pages/2_Generate.py", 
                 label=":orange[Generate:]", 
                 icon=":material/network_node:")

    st.markdown("- Generate the RNA sequence that matches the desired dot-bracket"
                " notation for the nanostructure.")

    st.page_link("pages/3_Convert.py", 
                 label=":orange[Convert:]", 
                 icon=":material/genetics:")

    st.markdown("- Prepare the DNA template for you RNA Origami, search subsequences"
                " and search for dimers.")

    st.page_link("pages/4_Prepare.py", 
                 label=":orange[Prepare:]", 
                 icon=":material/sync_alt:")

    st.markdown("- Design primers for your DNA template or prepare the Origami for "
                "OxDNA simulation.")
    
    
    _, col, _ = st.columns(3)
    with col:
        quickstart()
    
    st.divider()

    st.write("### About pyFuRNAce")
    
    st.markdown("pyFuRNAce is an open-source Python package and web-based design "
                "engine for creating complex RNA nanostructures using the "
                " co-transcriptional RNA origami approach.")
    st.markdown(" - **GitHub**: [Biophysical-Engineering-Group/pyFuRNAce]"
                "(https://github.com/Biophysical-Engineering-Group/pyFuRNAce)")
    st.markdown(" - **PyPI**: [pyfurnace](https://pypi.org/project/pyfurnace/)")
    st.markdown(" - **Documentation**: [Read the Docs]"
                "(https://pyfurnace.readthedocs.io/en/latest/)")
    st.markdown(" - bug reports, feature requests or any other questions, "
                "please reach out to us via the " 
                "[GitHub Issues](https://github.com/Biophysical-Engineering"
                "-Group/pyFuRNAce/issues)"
                " or the "
                "[GitHub Discussions](https://github.com/Biophysical-Engineering"
                "-Group/pyFuRNAce/discussions)."
                )
    
    st.divider()
    
    st.write("#### Check out the 1-min demo video:")
    st.video("https://github.com/Biophysical-Engineering-Group/pyFuRNAce/blob"
             "/main/vid/demo_1min.mp4?raw=true", 
             format="video/mp4", 
             start_time=0, 
             subtitles=None, 
             loop=True, 
             autoplay=True, 
             muted=True)