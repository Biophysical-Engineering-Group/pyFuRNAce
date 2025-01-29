import inspect
import streamlit as st
import hydralit_components as hc
from pyfurnace.design import motifs, Loop, structural
from .motif_command import MotifCommand
from .dovetail_command import DovetailCommand
from .stem_command import StemCommand
from .. import second_hc_theme
from ..motifs_icons import MOTIF_ICONS


struct_names = [ut_name for ut_name, obj in inspect.getmembers(structural) if inspect.isfunction(obj)]
          
struct_names = ['Stem', 'Dovetail'] + struct_names

option_data = [{'icon': MOTIF_ICONS[name], 
                'label': name.replace('_', ' ')} for name in struct_names]

class StructuralCommand(MotifCommand):

    def execute(self):
        # override the theme
        selected = hc.option_bar(option_definition=option_data, key='StructuralOption', override_theme=second_hc_theme, horizontal_orientation=True)
        match selected:
            case 'Stem':
                StemCommand().execute()
            case 'Dovetail':
                DovetailCommand().execute()
            case _:
                motif = [func() for name, func in motifs.__dict__.items() if callable(func) and name == selected][0]
                flip_default = False
                if st.session_state.current_line_occupied and isinstance(motif, Loop):
                    flip_default = True
                st.session_state.motif = motif
                st.session_state.motif_buffer = f"motif = pf.{selected}()"
                if st.toggle("Flip the aptamer", value=flip_default, key='flip_aptamer'):
                    motif.flip(1, 1)
                    st.session_state.motif_buffer += ".flip(1, 1)"