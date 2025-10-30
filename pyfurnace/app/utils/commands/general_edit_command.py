import streamlit as st
from streamlit import session_state as st_state
from .motif_command import MotifCommand
from .. import pyfurnace_layout_cols


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
        col1, col2 = pyfurnace_layout_cols([1.2, 2], vertical_alignment="bottom")

        if key:
            with col1:
                subcol1, subcol2 = st.columns(2)
                with subcol1:
                    flip_vert = st.toggle(
                        "Flip vertically",
                        value=False,
                    )
                with subcol2:
                    flip_hor = st.toggle(
                        "Flip horizontally",
                        value=False,
                    )
            with col2:
                subcol1, subcol2 = st.columns(2)
                if flip_vert or flip_hor:
                    with subcol1:
                        if flip_vert or flip_hor:
                            flip = st.button("Flip")
                            flip_vert &= flip
                            flip_hor &= flip
                with subcol2:
                    rotate = st.button("Rotate 90° clockwise")

        else:
            with col1:
                subcol1, subcol2 = st.columns(2)
                with subcol1:
                    flip_vert = st.toggle("Flip vertically")
                with subcol2:
                    flip_hor = st.toggle("Flip horizontally")
            with col2:
                subcol1, subcol2 = st.columns(2)
                with subcol2:
                    rotate = st.number_input(
                        "Rotate 90° clockwise:",
                        min_value=0,
                        max_value=4,
                    )

        return flip_vert, flip_hor, rotate
