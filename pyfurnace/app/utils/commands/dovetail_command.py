import streamlit as st
import pyfurnace as pf
from .motif_command import MotifCommand


class DovetailCommand(MotifCommand):

    def execute(self, motif=None):
        ### Modify the motif
        if motif:
            sign, top_seq, seq_length, wobble_interval, wobble_tolerance, wobble_insert, top_cross, bot_cross, strong_bases = self.interface('mod', motif[0].sequence, motif.length, motif.wobble_interval, motif.wobble_tolerance, motif.wobble_insert, motif.top_cross, motif.bot_cross, 'S' in motif.sequence)
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
            elif motif.wobble_insert != wobble_insert:
                st.session_state.modified_motif_text += f"\nmotif.wobble_insert = '{wobble_insert}'"
                motif.wobble_insert = wobble_insert
            elif motif.top_cross != top_cross:
                st.session_state.modified_motif_text += f"\nmotif.top_cross = {top_cross}"
                motif.top_cross = top_cross
            elif motif.bot_cross != bot_cross:
                st.session_state.modified_motif_text += f"\nmotif.bot_cross = {bot_cross}"
                motif.bot_cross = bot_cross
            elif strong_bases and 'N' in motif.sequence or not strong_bases and 'S' in motif.sequence:
                st.session_state.modified_motif_text += f"\nmotif.set_strong_bases({strong_bases})"
                motif.set_strong_bases(strong_bases)
        ### Create a new motif
        else:
            # take the values from the interface
            sign, top_seq, seq_length, wobble_interval, wobble_tolerance, wobble_insert, top_cross, bot_cross, strong_bases = self.interface()

            # Eventually remove the top and bottom cross from the motif
            if not top_cross:
                st.session_state.motif.top_cross = False
            if not bot_cross:
                st.session_state.motif.bot_cross = False
            if top_seq:
                st.session_state.motif_buffer = f"motif = pf.Dovetail(sequence = '{top_seq}', sign = str({sign}), top_cross = {top_cross}, bot_cross = {bot_cross})"
                motif = pf.Dovetail(sequence = top_seq, sign=sign, top_cross=top_cross, bot_cross=bot_cross)
            else:
                st.session_state.motif_buffer = f"motif = pf.Dovetail(length = {seq_length}, wobble_interval = {wobble_interval}, wobble_tolerance = {wobble_tolerance}, wobble_insert = '{wobble_insert}', top_cross = {top_cross}, bot_cross = {bot_cross}, strong_bases = {strong_bases})"
                motif = pf.Dovetail(length = seq_length, sign=sign, wobble_interval = wobble_interval, wobble_insert = wobble_insert, top_cross=top_cross, bot_cross=bot_cross, strong_bases=strong_bases)
            # save the motif in the session state
            st.session_state.motif = motif

    def interface(self, key='', top_seq=None, len_default=-2, wobble_interval=8, wobble_tolerance=2, wobble_insert='middle', top_cross=True, bot_cross=True, strong_bases=True):
        ### initialize the variables
        seq_length = None
        sign = 0

        ### create the interface
        col1, col2 = st.columns([1, 5], vertical_alignment='bottom')
        with col1:
            specific_seq = st.toggle("Custom Sequence", key=f'seq_dt{key}')
        with col2:
            if specific_seq:   
                col1, col2, col3, col4 = st.columns([4, 1, 1, 1])
                with col1:
                    top_seq = st.text_input('Sequence:', key=f'txt_top_seq_dt{key}', value=top_seq)
                with col2:
                    sign = st.selectbox('Sign:', [-1, +1], key=f'txt_sign_dt{key}')
                with col3:
                    top_cross = st.toggle('Top cross', key=f'top_cross_dt{key}', value=top_cross)
                with col4:
                    bot_cross = st.toggle('Bottom cross', key=f'bot_cross_dt{key}', value=bot_cross)
            else:
                subcol1, subcol2, subcol3, subcol4, subcol5, subcol6, subcol7, = st.columns(7, 
                                                                                            vertical_alignment='bottom')
                with subcol1:
                    seq_length = st.number_input('Length:', key=f'dt_length{key}', min_value=-100, value=len_default)
                with subcol2:
                    wobble_interval = st.number_input('Wobble interval:', key=f'dt_wobble_interval{key}', min_value=0, value=wobble_interval, help="Add a wobble every n° nucleotides")
                with subcol3:
                    wobble_tolerance = st.number_input('Wobble tolerance:', key=f'dt_wobble_tolerance{key}', min_value=0, value=wobble_tolerance)
                with subcol4:
                    wobble_insert = st.selectbox('Wobble insert:', ['middle', 'start', 'end'], key=f'dt_wobble_ins{key}', index=['middle', 'start', 'end'].index(wobble_insert))
                with subcol5:
                    strong_bases = st.toggle('Strong bases', key=f'strong_bases_dt{key}', value=strong_bases, help='Use strong bases for dovetails shoter than 4 nt')
                with subcol6:
                    top_cross = st.toggle('Top cross', key=f'top_cross_dt{key}', value=top_cross)
                with subcol7:
                    bot_cross = st.toggle('Bottom cross', key=f'bot_cross_dt{key}', value=bot_cross)

        return sign, top_seq, seq_length, wobble_interval, wobble_tolerance, wobble_insert, top_cross, bot_cross, strong_bases
