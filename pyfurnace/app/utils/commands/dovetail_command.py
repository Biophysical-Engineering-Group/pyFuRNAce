import streamlit as st
from streamlit import session_state as st_state
import pyfurnace as pf
from .motif_command import MotifCommand


class DovetailCommand(MotifCommand):

    def execute(self, motif=None):
        ### Modify the motif
        if motif:
            # interface data:
            # sign,top_seq,length,w_interval,w_tolerance, up_cross, down_cross
            data = self.interface(
                motif[0].sequence,
                motif.length,
                motif.wobble_interval,
                motif.wobble_tolerance,
                motif.up_cross,
                motif.down_cross,
            )
            sign, top_seq, seq_length = data[:3]
            wobble_interval, wobble_tolerance = data[3:5]
            up_cross, down_cross = data[5:7]

            if (
                top_seq
                and str(motif[0].sequence) != top_seq
                and str(motif[1].sequence) != top_seq
                or top_seq
                and sign
                and sign != motif._sign
            ):
                st_state.modified_motif_text += (
                    f"\nmotif.set_up_sequence('" f"{top_seq}', sign = str({sign})"
                )
                motif.set_up_sequence(top_seq, sign)
            elif isinstance(seq_length, int) and motif.length != seq_length:
                st_state.modified_motif_text += f"\nmotif.length = " f"{seq_length}"
                motif.length = seq_length

            elif motif.wobble_interval != wobble_interval:
                st_state.modified_motif_text += (
                    f"\nmotif.wobble_interval = " f"{wobble_interval}"
                )
                motif.wobble_interval = wobble_interval

            elif motif.wobble_tolerance != wobble_tolerance:
                st_state.modified_motif_text += (
                    f"\nmotif.wobble_tolerance = " f"{wobble_tolerance}"
                )
                motif.wobble_tolerance = wobble_tolerance

            elif motif.up_cross != up_cross:
                st_state.modified_motif_text += f"\nmotif.up_cross = " f"{up_cross}"
                motif.up_cross = up_cross

            elif motif.down_cross != down_cross:
                st_state.modified_motif_text += f"\nmotif.down_cross = " f"{down_cross}"
                motif.down_cross = down_cross

        ### Create a new motif
        else:
            # take the values from the interface
            data = self.interface()
            sign, top_seq, seq_length = data[:3]
            wobble_interval, wobble_tolerance = data[3:5]
            up_cross, down_cross = data[5:7]

            # Eventually remove the top and bottom cross from the motif
            if not up_cross:
                st_state.motif.up_cross = False
            if not down_cross:
                st_state.motif.down_cross = False
            if top_seq:
                st_state.motif_buffer = (
                    f"motif = pf.Dovetail(sequence = "
                    f"'{top_seq}', sign = str({sign}), "
                    f"up_cross = {up_cross}, "
                    f"down_cross = {down_cross})"
                )
                motif = pf.Dovetail(
                    sequence=top_seq,
                    sign=sign,
                    up_cross=up_cross,
                    down_cross=down_cross,
                )
            else:
                st_state.motif_buffer = (
                    f"motif = pf.Dovetail(length = {seq_length}, "
                    f"wobble_interval = {wobble_interval}, "
                    f"wobble_tolerance = {wobble_tolerance}, "
                    f"up_cross = {up_cross}, "
                    f"down_cross = {down_cross}, "
                    f"wobble_insert = 'middle')"
                )
                motif = pf.Dovetail(
                    length=seq_length,
                    sign=sign,
                    wobble_interval=wobble_interval,
                    up_cross=up_cross,
                    down_cross=down_cross,
                )
            # save the motif in the session state
            st_state.motif = motif

    def interface(
        self,
        top_seq=None,
        len_default=-2,
        wobble_interval=5,
        wobble_tolerance=2,
        up_cross=True,
        down_cross=True,
    ):
        ### initialize the variables
        seq_length = 0
        sign = 0

        ### create the interface
        col1, col2 = st.columns([1, 5], vertical_alignment="bottom")
        with col1:
            specific_seq = st.toggle("Custom Sequence")
        with col2:
            if specific_seq:
                col1, col2, col3, col4 = st.columns(
                    [4, 1, 1, 1], vertical_alignment="bottom"
                )
                with col1:
                    top_seq = st.text_input("Sequence:", value=top_seq)
                with col2:
                    sign = st.selectbox("Sign:", [-1, +1])
                with col3:
                    up_cross = st.toggle("Up cross", value=up_cross)
                with col4:
                    down_cross = st.toggle("Down cross", value=down_cross)
            else:
                subcols = st.columns([1, 1, 1, 1, 1], vertical_alignment="bottom")
                with subcols[0]:
                    seq_length = st.number_input(
                        "Length:",
                        min_value=-100,
                        value=len_default,
                    )
                with subcols[1]:
                    up_cross = st.toggle("Up cross", value=up_cross)
                with subcols[2]:
                    down_cross = st.toggle("Down cross", value=down_cross)
                with subcols[3]:
                    wobble_interval = st.number_input(
                        "Wobble interval:",
                        min_value=0,
                        value=wobble_interval,
                        help="Add a wobble every n° nucleotides, to avoid secondary "
                        "structures in the DNA",
                    )
                with subcols[4]:
                    wobble_tolerance = st.number_input(
                        "Wobble tolerance:",
                        min_value=0,
                        value=wobble_tolerance,
                        help="Randomize the number of wobbles by ± n, to speed up "
                        "sequence optimization",
                    )

        return (
            sign,
            top_seq,
            seq_length,
            wobble_interval,
            wobble_tolerance,
            up_cross,
            down_cross,
        )
