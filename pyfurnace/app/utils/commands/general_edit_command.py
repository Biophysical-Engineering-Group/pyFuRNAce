import streamlit as st
from streamlit import session_state as st_state
from .motif_command import MotifCommand


class GeneralEditCommand(MotifCommand):

    def execute(self, motif=None):
        ### Modify motif
        if motif:
            flip_vert, flip_hor, rotate = self.interface("mod")
            if flip_vert or flip_hor:
                st_state.modified_motif_text += (
                    f"\nmotif.flip(horizontally="
                    f"{flip_hor}, "
                    f"vertically={flip_vert})"
                )
                motif.flip(horizontally=flip_hor, vertically=flip_vert)
            elif rotate:
                st_state.modified_motif_text += f"\nmotif.rotate({rotate})"
                motif.rotate(rotate)

    @staticmethod
    def interface(key=""):
        col1, col2, col3 = st.columns(3, vertical_alignment="bottom")

        if key:
            with col1:
                subcol1, subcol2 = st.columns(2)
                with subcol1:
                    flip_vert = st.toggle(
                        "Flip vertically",
                        value=True,
                    )
                with subcol2:
                    flip_hor = st.toggle(
                        "Flip horizontally",
                        value=True,
                    )
            with col2:
                if flip_vert or flip_hor:
                    flip = st.button("Flip")
                    flip_vert &= flip
                    flip_hor &= flip
            with col3:
                rotate = st.button("Rotate 90° clockwise")

        else:
            with col1:
                flip_vert = st.toggle("Flip vertically")
            with col2:
                flip_hor = st.toggle("Flip horizontally")
            with col3:
                rotate = st.number_input(
                    "Rotate 90° clockwise:",
                    min_value=0,
                    max_value=4,
                )

        return flip_vert, flip_hor, rotate
