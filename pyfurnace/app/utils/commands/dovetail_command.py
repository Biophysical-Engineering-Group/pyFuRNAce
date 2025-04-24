import streamlit as st
import pyfurnace as pf
from .motif_command import MotifCommand


class DovetailCommand(MotifCommand):

    def execute(self, motif=None):
        ### Modify the motif
        if motif:
            sign, top_seq, seq_length, wobble_interval, wobble_tolerance, up_cross, down_cross = self.interface('mod', 
                                                                                                                motif[0].sequence, 
                                                                                                                motif.length, 
                                                                                                                motif.wobble_interval, 
                                                                                                                motif.wobble_tolerance,
                                                                                                                motif.up_cross, 
                                                                                                                motif.down_cross)
            if top_seq and str(motif[0].sequence) != top_seq and str(motif[1].sequence) != top_seq or top_seq and sign and sign != motif._sign:
                st.session_state.modified_motif_text += f"\nmotif.set_top_sequence('{top_seq}', sign = str({sign})"
                motif.set_top_sequence(top_seq, sign)
            elif type(seq_length) == int and motif.length != seq_length:
                st.session_state.modified_motif_text += f"\nmotif.length = {seq_length}"
                motif.length = seq_length
            elif motif.wobble_interval != wobble_interval:
                st.session_state.modified_motif_text += f"\nmotif.wobble_interval = {wobble_interval}"
                motif.wobble_interval = wobble_interval
            elif motif.wobble_tolerance != wobble_tolerance:
                st.session_state.modified_motif_text += f"\nmotif.wobble_tolerance = {wobble_tolerance}"
                motif.wobble_tolerance = wobble_tolerance
            elif motif.up_cross != up_cross:
                st.session_state.modified_motif_text += f"\nmotif.up_cross = {up_cross}"
                motif.up_cross = up_cross
            elif motif.down_cross != down_cross:
                st.session_state.modified_motif_text += f"\nmotif.down_cross = {down_cross}"
                motif.down_cross = down_cross
        ### Create a new motif
        else:
            # take the values from the interface
            sign, top_seq, seq_length, wobble_interval, wobble_tolerance, up_cross, down_cross = self.interface()

            # Eventually remove the top and bottom cross from the motif
            if not up_cross:
                st.session_state.motif.up_cross = False
            if not down_cross:
                st.session_state.motif.down_cross = False
            if top_seq:
                st.session_state.motif_buffer = f"motif = pf.Dovetail(sequence = '{top_seq}', sign = str({sign}), up_cross = {up_cross}, down_cross = {down_cross})"
                motif = pf.Dovetail(sequence = top_seq, sign=sign, up_cross=up_cross, down_cross=down_cross)
            else:
                st.session_state.motif_buffer = f"motif = pf.Dovetail(length = {seq_length}, wobble_interval = {wobble_interval}, wobble_tolerance = {wobble_tolerance}, up_cross = {up_cross}, down_cross = {down_cross}, wobble_insert = 'middle')"
                motif = pf.Dovetail(length = seq_length, sign=sign, wobble_interval = wobble_interval, up_cross=up_cross, down_cross=down_cross)
            # save the motif in the session state
            st.session_state.motif = motif

    def interface(self, key='', top_seq=None, len_default=-2, wobble_interval=8, wobble_tolerance=2, up_cross=True, down_cross=True):
        ### initialize the variables
        seq_length = 0
        sign = 0

        ### create the interface
        col1, col2 = st.columns([1, 5], vertical_alignment='bottom')
        with col1:
            specific_seq = st.toggle("Custom Sequence", key=f'seq_dt{key}')
        with col2:
            if specific_seq:   
                col1, col2, col3, col4 = st.columns([4, 1, 1, 1], vertical_alignment='bottom')
                with col1:
                    top_seq = st.text_input('Sequence:', key=f'txt_top_seq_dt{key}', value=top_seq)
                with col2:
                    sign = st.selectbox('Sign:', [-1, +1], key=f'txt_sign_dt{key}')
                with col3:
                    up_cross = st.toggle('Up cross', key=f'up_cross_dt{key}', value=up_cross)
                with col4:
                    down_cross = st.toggle('Down cross', key=f'down_cross_dt{key}', value=down_cross)
            else:
                subcol1, subcol2, subcol3, subcol4, subcol5 = st.columns([1, 1, 1, 1, 1], vertical_alignment='bottom')
                with subcol1:
                    seq_length = st.number_input('Length:', key=f'dt_length{key}', min_value=-100, value=len_default)
                with subcol2:
                    up_cross = st.toggle('Up cross', key=f'up_cross_dt{key}', value=up_cross)
                with subcol3:
                    down_cross = st.toggle('Down cross', key=f'down_cross_dt{key}', value=down_cross)
                with subcol4:
                    wobble_interval = st.number_input('Wobble interval:', key=f'dt_wobble_interval{key}', min_value=0, value=wobble_interval, help="Add a wobble every n° nucleotides")
                with subcol5:
                    wobble_tolerance = st.number_input('Wobble tolerance:', key=f'dt_wobble_tolerance{key}', min_value=0, value=wobble_tolerance)


        return sign, top_seq, seq_length, wobble_interval, wobble_tolerance, up_cross, down_cross
