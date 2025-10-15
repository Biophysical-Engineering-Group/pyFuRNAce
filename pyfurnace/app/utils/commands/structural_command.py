import inspect
import streamlit as st
from streamlit_option_menu import option_menu
from pyfurnace.design import motifs, Loop, structural
from .motif_command import MotifCommand
from .dovetail_command import DovetailCommand
from .stem_command import StemCommand
from .. import second_menu_style
from ..motifs_icons import MOTIF_ICONS


struct_names = [
    ut_name
    for ut_name, obj in inspect.getmembers(structural)
    if inspect.isfunction(obj)
]

# ignore the three-way junction
del struct_names[struct_names.index("ThreeWayJunction")]

struct_names = ["Stem", "Dovetail"] + struct_names


class StructuralCommand(MotifCommand):

    def execute(self):
        selected = option_menu(
            None,
            struct_names,
            icons=[MOTIF_ICONS[name] for name in struct_names],
            menu_icon="cast",
            orientation="horizontal",
            styles=second_menu_style,
            key="StructuralOption",
        )
        if selected == "Stem":
            StemCommand().execute()
            return
        elif selected == "Dovetail":
            DovetailCommand().execute()
            return

        motif = [
            func()
            for name, func in motifs.__dict__.items()
            if callable(func) and name == selected
        ][0]
        flip_default = False

        if st.session_state.current_line_occupied and isinstance(motif, Loop):
            flip_default = True

        st.session_state.motif = motif
        st.session_state.motif_buffer = f"motif = pf.{selected}()"

        if st.toggle("Flip the motif", value=flip_default, key="flip_motif"):
            motif.flip(1, 1)
            st.session_state.motif_buffer += ".flip(1, 1)"
