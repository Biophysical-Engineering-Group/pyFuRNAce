import streamlit as st
from .motif_command import MotifCommand

class GeneralEditCommand(MotifCommand):
    
    def execute(self, motif=None):
        ### Modify motif
        if motif:
            flip_vert, flip_hor, rotate = self.interface('mod')
            if flip_vert or flip_hor:
                st.session_state.modified_motif_text += f"\nmotif.flip(horizontally={flip_hor}, vertically={flip_vert})"
                motif.flip(horizontally=flip_hor, vertically=flip_vert)
            elif rotate:
                st.session_state.modified_motif_text += f"\nmotif.rotate({rotate})"
                motif.rotate(rotate)

    @staticmethod
    def interface(key=''):
        col1, col2, col3 = st.columns(3)
        with col1:
            if key:
                flip_vert = st.button("Flip vertically", key=f"{key}flip_vert")
            else:
                st.write('\n'); st.write('\n')
                flip_vert = st.toggle("Flip vertically", key=f"{key}flip_vert")
        with col2:
            if key:
                flip_hor = st.button("Flip horizontally", key=f"{key}flip_hor")
            else:
                st.write('\n'); st.write('\n')
                flip_hor = st.toggle("Flip horizontally", key=f"{key}flip_hor")
        with col3:
            if key:
                rotate = st.button("Rotate 90° clockwise", key=f"{key}rotate")
            else:
                rotate = st.number_input("Rotate 90° clockwise:", min_value=0, max_value=4, key=f"{key}rotate")
        return flip_vert, flip_hor, rotate
