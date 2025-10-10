import streamlit as st
from streamlit import session_state as st_state
import pyfurnace as pf
from .motif_command import MotifCommand


class TetraLoopCommand(MotifCommand):

    def execute(self, motif=None):
        ### Modify the motif
        if motif:
            seq, flip = self.interface("mod", motif[0].sequence)
            if seq and motif[0].sequence != seq:
                st_state.modified_motif_text += f"\nmotif.set_sequence('{seq}')"
                motif.set_sequence(seq)
            elif flip:
                st_state.modified_motif_text += "\nmotif.flip()"
                motif.flip()

        ### Create a new motif
        else:
            seq, open_left = self.interface()
            if seq:
                st_state.motif_buffer = (
                    f"motif = pf.TetraLoop({open_left}, " f"sequence = '{seq}')"
                )
                motif = pf.TetraLoop(open_left, sequence=seq)

            else:
                st_state.motif_buffer = f"motif = pf.TetraLoop({open_left})"
                motif = pf.TetraLoop(open_left)
            # save the motif in the session state
            st_state.motif = motif

    def interface(self, key="", seq_default="UUCG"):
        seq = None
        open_left = False
        col1, col2, col3 = st.columns([1, 1, 4], vertical_alignment="bottom")
        with col1:
            custom_seq = st.toggle("Custom Sequence")
        with col2:
            if key == "mod":
                open_left = st.button("Flip")
            else:
                open_left = st.toggle(
                    "Flip",
                    value=st_state.current_line_occupied,
                )
        if custom_seq:
            with col3:
                seq = st.text_input("Sequence:", value=seq_default)
        return seq, open_left
