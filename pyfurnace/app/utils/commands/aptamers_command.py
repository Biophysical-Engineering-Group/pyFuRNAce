import inspect
import streamlit as st
import hydralit_components as hc
from pyfurnace.design.motifs import aptamers
from pyfurnace.design.motifs import Loop
from .motif_command import MotifCommand
from .. import second_hc_theme
from ..motifs_icons import MOTIF_ICONS

# Get names of all aptamer functions defined in the `aptamers` module
aptamers_names = [
    func_name for func_name, member in inspect.getmembers(aptamers, inspect.isfunction)
    if member.__module__ == aptamers.__name__
]

common_aptamers = ['Broccoli',
                   'Pepper',
                   'Ispinach',
                   'MS2',
                   "Streptavidin",
                    ]


class AptamersCommand(MotifCommand):

    def execute(self):
        # override the theme
        col1, col2 = st.columns([5, 1])
        with col1:
            option_data = [{'icon': MOTIF_ICONS[name],
                            'label': name} for name in common_aptamers]
            aptamer_selction = hc.option_bar(option_definition=option_data, key='AptamerOption', override_theme=second_hc_theme, horizontal_orientation=True)
        with col2:
            aptamers_box = st.selectbox(":green[Or search an Aptamer]", 
                                        ['No selection'] + aptamers_names, 
                                        key='aptamers_ddd',
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