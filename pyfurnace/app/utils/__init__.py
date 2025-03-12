from pathlib import Path
from copy import deepcopy
import tempfile
import sys
import streamlit as st
import importlib.util

app_path = Path(__file__).resolve().parent.parent
main_menu_style = { "container": {
                        "padding": "0!important", 
                        "border-radius": "10px",  # Rounded borders
                        "position": "relative"  # Ensures it doesn't overflow
                    },
                    "nav-link": {
                        "font-size": "inherit", 
                        "text-align": "center", 
                        "margin": "0px", 
                        "border-radius": "10px", # Rounded option links
                        "--hover-color": "#cccfcc"  # Hover color
                    },
                    "nav-link-selected": {
                        "border-radius": "10px"  # Rounded selected option
                    }
                }
second_menu_style = deepcopy(main_menu_style)
second_menu_style['nav-link-selected']['background-color'] = '#D00000'


def check_import_pyfurnace():
    ### If the module is already imported, skip
    if 'pyfurnace' in sys.modules: 
        return
    
    ### The module is installed but not imported, import it
    if importlib.util.find_spec('pyfurnace') is not None: # Check if the module is installed
        import pyfurnace
        return
    
    ### Last chance, try to import the module from the local path
    pyfurnace_path = app_path.parent.parent
    sys.path.insert(0, str(pyfurnace_path))
    import pyfurnace

def load_logo(page_title="pyFuRNAce", page_icon=str(app_path / "static" / "logo_lr.png")):
    # First instructions to run the app, set the layout and logo
    st.set_page_config(page_title=page_title, page_icon=page_icon, layout="wide", initial_sidebar_state='collapsed',)
    st.logo(str(app_path / "static" / "logo_text.png"), 
            icon_image=str(app_path / "static" / "logo_lr.png"),
            link='https://pyfurnace.streamlit.app',
            size='large')
    # st.html("""
    #     <style>
    #       [alt=Logo] {
    #         height: 3rem;
    #       }
    #     </style>
    #     """)

def copy_to_clipboard(text_to_copy, button_text=''):
    copied_text = 'Copied!'
    if not button_text:
        copied_text = ''
    st.components.v1.html(f"""
    <link href="https://fonts.googleapis.com/css2?family=Material+Symbols+Outlined" rel="stylesheet" />

        <script>
            function copyToClipboard() {{
                button = document.querySelector('.copy_button');
                navigator.clipboard.writeText("{text_to_copy}")
                button.innerHTML = '<span class="material-symbols-outlined">done</span>{copied_text}';
                setTimeout(() => {{
                     button.style.backgroundColor = 'white';  // Change label after copying
                     button.innerHTML = '<span class="material-symbols-outlined">content_copy</span>{button_text}';
                }}, 1000);
        }}
        </script>
        <style>
            .copy_container {{
                display: inline-block;  /* Ensures no extra space around */
                margin: 0;
                padding: 0;
            }}
            
            .copy_button {{
                display: inline-flex;
                background-color: inherit;
                margin: 0rem;
                margin-left: -0.5rem;
                margin-top: -0.5rem;
                align-items: center;
                border-radius: 0.5rem;
                font-size: inherit;
                font-weight: 350;
                cursor: pointer;
                border: none;
                }}
            .copy_button:hover {{
                color: #00856A;
                border-color: #00856A;
                }}

        </style>
        <div class="copy_container">
            <button class="stButton copy_button" onclick="copyToClipboard()">
                <span class="material-symbols-outlined">content_copy</span>
                {button_text}
            </button>
        </div>
        """, height=25)

def save_origami(origami_name='Origami'):
    if not st.session_state.origami:
        return
    origami = st.session_state.origami

    st.write('### Download RNA origami structure')
    col1, col2 = st.columns([1, 6])
    with col1:
        file_type = st.selectbox('Select file type', ['py', 'txt', 'txt with info', 'PDB', 'oxDNA'], key = 'file_type')
    with col2:
        ori_name = st.text_input('Name of RNA origami', value=origami_name, key = 'ori_name')
        if not ori_name:
            st.stop()

        if file_type == 'PDB':
            if any(nucl not in "AUCG&" for nucl in origami.sequence):
                st.error('The sequence contains non-standard nucleotides. PDB format only supports A, U, C, G.')
                st.stop()
            with tempfile.TemporaryDirectory() as tmpdirname:
                file_path = f"{tmpdirname}/origami"
                origami.save_3d_model(file_path, pdb=True)
                try:
                    with open(f"{file_path}.pdb", 'r') as f:
                        pdb_text = f.read()
                except:
                    st.warning('No PDB file found')
                    pdb_text = None
            if pdb_text:
                st.download_button('Download PDB', 
                                   pdb_text, 
                                   f"{ori_name}.pdb",
                                   on_click="ignore")
        
        elif file_type == 'oxDNA':
            with tempfile.TemporaryDirectory() as tmpdirname:
                file_path = f"{tmpdirname}/origami"
                origami.save_3d_model(file_path, forces=True)
                with open(f"{file_path}.dat", 'r') as f:
                    conf_text = f.read()
                with open(f"{file_path}.top", 'r') as f:
                    topo_text = f.read()
                try:
                    with open(f"{file_path}_forces.txt", 'r') as f:
                        forces = f.read()
                except:
                    forces = None
                    st.warning('No forces file found')
            col1, col2, col3 = st.columns(3)
            with col1:
                st.download_button('Download Configuration', 
                                   conf_text, 
                                   f"{ori_name}.dat",
                                   on_click="ignore")
            with col2:
                st.download_button('Download Topology', 
                                   topo_text, 
                                   f"{ori_name}.top",
                                   on_click="ignore")
            with col3:
                if forces:
                    st.download_button('Download Forces', 
                                       forces, 
                                       f"{ori_name}_forces.txt",
                                       on_click="ignore")

        else:
            if file_type == 'py' and 'code' in st.session_state:
                # create a text data with the structure of the RNA origami in python code
                text_data = '\n\n'.join(st.session_state.code)
            elif 'txt' in file_type:
                #create a text data with the structure of the RNA origami
                if 'to_road' in st.session_state and st.session_state.to_road:
                    origami_str = origami.to_road()
                else:
                    origami_str = str(origami)
                info = ''
                if file_type == 'txt with info':
                    info = f"\nSequence:{origami.sequence}\nStructure:\n{origami.structure}\nPseudoknots info:\n{origami.pseudoknots}\n"
                    file_type = 'txt'
                text_data = f'>{ori_name}\n{info}\n{origami_str}'
                
            st.download_button('Download', 
                               text_data, 
                               f"{ori_name}.{file_type}",
                               on_click="ignore")
    



