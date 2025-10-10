import streamlit as st
from streamlit_option_menu import option_menu
from pyfurnace.design.motifs import aptamers, aptamers_list
from pyfurnace.design.motifs import Loop
from .motif_command import MotifCommand
from .. import second_menu_style, inactive_menu_style
from ..motifs_icons import MOTIF_ICONS

common_aptamers = [
    "Broccoli",
    "Pepper",
    "Ispinach",
    "MS2",
    "Streptavidin",
]


class AptamersCommand(MotifCommand):

    def execute(self):
        # override the theme
        col1, col2 = st.columns([1, 5], vertical_alignment="bottom")
        with col1:
            aptamers_box = st.selectbox(
                ":green[Search an Aptamer or]",
                ["No selection"] + aptamers_list,
                key="aptamers_box",
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
            menu_style = second_menu_style
            if aptamers_box != "No selection":
                menu_style = inactive_menu_style
            aptamer_selection = option_menu(
                None,
                common_aptamers,
                icons=[MOTIF_ICONS[name] for name in common_aptamers],
                menu_icon="cast",
                orientation="horizontal",
                styles=menu_style,
            )

        if aptamers_box != "No selection":
            aptamer_selection = aptamers_box
        if aptamer_selection:
            motif = aptamers.__dict__[aptamer_selection]()
            flip_default = False
            if st.session_state.current_line_occupied and isinstance(motif, Loop):
                flip_default = True
            st.session_state.motif = motif
            st.session_state.motif_buffer = f"motif = pf.{aptamer_selection}()"
            if st.toggle("Flip the aptamer", value=flip_default):
                motif.flip(1, 1)
                st.session_state.motif_buffer += ".flip(1, 1)"
