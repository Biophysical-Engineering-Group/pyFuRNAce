import streamlit as st
import pyfurnace as pf
from .motif_command import MotifCommand


class StemCommand(MotifCommand):

    def execute(self, motif=None):
        ### Modify the motif
        if motif:
            top_seq, seq_length, wobble_interval, wobble_tolerance = self.interface('mod', 
                                                                                    motif[0].sequence, 
                                                                                    motif.length, 
                                                                                    motif.wobble_interval, 
                                                                                    motif.wobble_tolerance)
            if top_seq and motif[0].sequence != top_seq:
                st.session_state.modified_motif_text += f"\nmotif.set_top_sequence('{top_seq}')"
                motif.set_top_sequence(top_seq)
            elif seq_length and motif.length != seq_length:
                st.session_state.modified_motif_text += f"\nmotif.length = {seq_length}"
                motif.length = seq_length
            elif motif.wobble_interval != wobble_interval:
                st.session_state.modified_motif_text += f"\nmotif.wobble_interval = {wobble_interval}"
                motif.wobble_interval = wobble_interval
            elif motif.wobble_tolerance != wobble_tolerance:
                st.session_state.modified_motif_text += f"\nmotif.wobble_tolerance = {wobble_tolerance}"
                motif.wobble_tolerance = wobble_tolerance

        ### Create a new motif
        else:
            top_seq, seq_length, wobble_interval, wobble_tolerance = self.interface()
            if top_seq:
                st.session_state.motif_buffer = f"motif = pf.Stem(sequence = '{top_seq}')"
                motif = pf.Stem(sequence = top_seq)
            else:
                st.session_state.motif_buffer = f"motif = pf.Stem(length = {seq_length}, wobble_interval = {wobble_interval}, wobble_tolerance = {wobble_tolerance}, wobble_insert = 'middle', strong_bases = True)"
                motif = pf.Stem(length = seq_length, wobble_interval = wobble_interval, wobble_tolerance=wobble_tolerance)
            # save the motif in the session state
            st.session_state.motif = motif


    def interface(self, key='', top_seq=None, len_default=7, wobble_interval=7, wobble_tolerance=3):
        ### initialize the variables
        seq_length = 0

        ### create the interface
        col1, col2 = st.columns([1, 5], vertical_alignment='bottom')
        with col1:
            specific_seq = st.toggle("Custom Sequence", key=f'seq_stem{key}')
        with col2:
            if specific_seq:   
                col1, col2 = st.columns([5, 1])
                with col1:
                    top_seq = st.text_input('Sequence:', key=f'txt_top_seq_stem{key}', value=top_seq)
            else:
                subcol1, subcol2, subcol3= st.columns([3, 1, 1])
                with subcol1:
                    seq_length = st.number_input('Length:', key=f'stem_length{key}', min_value=1, value=len_default)
                with subcol2:
                    wobble_interval = st.number_input('Wobble interval:', key=f'stem_wobble_interval{key}', min_value=0, value=wobble_interval, help="Add a wobble every n° nucleotides")
                with subcol3:
                    wobble_tolerance = st.number_input('Wobble tolerance:', key=f'stem_wobble_tolerance{key}', min_value=0, value=wobble_tolerance)
        return top_seq, seq_length, wobble_interval, wobble_tolerance