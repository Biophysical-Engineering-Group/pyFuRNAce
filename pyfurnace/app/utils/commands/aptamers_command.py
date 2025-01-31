import inspect
import streamlit as st
from streamlit_option_menu import option_menu
from pyfurnace.design.motifs import aptamers
from pyfurnace.design.motifs import Loop
from .motif_command import MotifCommand
from .. import second_menu_style
from ..motifs_icons import MOTIF_ICONS

# Get names of all aptamer functions defined in the `aptamers` module
aptamers_names = [
    func_name for func_name, member in inspect.getmembers(aptamers, inspect.isfunction)
    if member.__module__ == aptamers.__name__
]

aptamers_names.remove('create_aptamer')

common_aptamers = ['Broccoli',
                   'Pepper',
                   'Ispinach',
                   'MS2',
                   "Streptavidin",
                    ]


class AptamersCommand(MotifCommand):

    def execute(self):
        # override the theme
        col1, col2 = st.columns([1, 5])
        with col1:
            aptamers_box = st.selectbox(":green[Search an Aptamer or]", 
                                        ['No selection'] + aptamers_names, 
                                        key='aptamers_ddd',
                                        )
        with col2:
            st.markdown(
                """
                <div style="text-align: center;">
                    Select a common aptamer:
                </div>
                """,
                unsafe_allow_html=True,
            )
            aptamer_selction = option_menu(None,
                                           common_aptamers,
                                           icons=[MOTIF_ICONS[name] for name in common_aptamers],
                                           menu_icon="cast",
                                           orientation="horizontal",
                                           styles=second_menu_style,
                                           key='AptamerOption',
                                           )

        if aptamers_box != 'No selection':
            aptamer_selction = aptamers_box
        if aptamer_selction:
            motif = [func() for name, func in aptamers.__dict__.items() if callable(func) and name == aptamer_selction][0]
            flip_default = False
            if st.session_state.current_line_occupied and isinstance(motif, Loop):
                flip_default = True
            st.session_state.motif = motif
            st.session_state.motif_buffer = f"motif = pf.{aptamer_selction}()"
            if st.toggle("Flip the aptamer", value=flip_default, key='flip_aptamer'):
                motif.flip(1, 1)
                st.session_state.motif_buffer += ".flip(1, 1)"