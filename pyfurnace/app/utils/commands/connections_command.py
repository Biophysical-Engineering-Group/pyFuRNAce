import inspect
from pyfurnace.design import utils, single_strand, start_end_stem

import streamlit as st
from streamlit_option_menu import option_menu
from .motif_command import MotifCommand
from .general_edit_command import GeneralEditCommand
from .tetraloop_command import TetraLoopCommand
from .. import second_menu_style
from ..motifs_icons import MOTIF_ICONS

# Filter and collect the motif utils
util_names = [
    ut_name
    for ut_name, obj in inspect.getmembers(utils.motif_lib)
    if inspect.isfunction(obj)
]

# ignore the simple vertical connection
del util_names[util_names.index("vertical_link")]

# Add the tetraloop option
util_names = ["Tetraloop"] + util_names


class ConnectionsCommand(MotifCommand):

    def execute(self, motif=None):
        ### Modify motif
        if motif:
            GeneralEditCommand().execute(motif)

        ### Create new motif
        else:
            util_option = option_menu(
                None,
                util_names,
                icons=[MOTIF_ICONS[name] for name in util_names],
                menu_icon="cast",
                orientation="horizontal",
                styles=second_menu_style,
                key="UtilsOption",
            )

            if util_option == "Tetraloop":
                TetraLoopCommand().execute()
                return
            if util_option == "single_strand":
                single_strandCommand().execute()
                return
            if util_option == "start_end_stem":
                start_end_stemCommand().execute()
                return

            name = util_option.replace(" ", "_")
            motif_util = getattr(utils.motif_lib, name)
            motif_text = f"motif = pf.{name}"
            flip_vert, flip_hor, rotate = self.interface()
            if flip_vert or flip_hor:
                motif_text += f".flip({flip_vert}, {flip_hor})"
            if rotate:
                motif_text += f".rotate({rotate})"
            # save the motif in the session state
            st.session_state.motif_buffer = (
                f"motif = pf.{name}(hflip={flip_hor}, "
                f"vflip={flip_vert}, rotate={rotate})"
            )
            st.session_state.motif = motif_util(
                hflip=flip_hor, vflip=flip_vert, rotate=rotate
            )

    def interface(self, key=""):
        return GeneralEditCommand.interface(key)


class single_strandCommand(MotifCommand):
    def execute(self):
        seq, loop, flip = self.interface()
        st.session_state.motif_buffer = (
            f"motif = pf.single_strand(sequence='{seq}', "
            "loop={loop}, vflip={flip}, hflip={flip})"
        )
        st.session_state.motif = single_strand(
            sequence=seq, loop=loop, vflip=flip, hflip=flip
        )

    def interface(
        self,
        seq_def="",
        loop_def=False,
    ):
        col1, col2, col3 = st.columns([3, 1, 1], vertical_alignment="bottom")
        with col1:
            seq = st.text_input(
                "Sequence:",
                value=seq_def,
            )
        with col2:
            loop = st.toggle(
                "Make a loop",
                value=loop_def,
                help="If checked, the strand will form a loop.",
            )
        with col3:
            flip = st.toggle("Flip", value=False, help="Flip the strand diagonally.")
        return seq, loop, flip


class start_end_stemCommand(MotifCommand):

    def execute(self):
        t_l, t_r, b_l, b_r = self.interface()
        st.session_state.motif_buffer = (
            f"motif = pf.start_end_stem(up_left='{t_l}', "
            f"up_right='{t_r}', down_left='{b_l}', "
            f"down_right='{b_r}')"
        )
        st.session_state.motif = start_end_stem(
            up_left=t_l, up_right=t_r, down_left=b_l, down_right=b_r
        )

    def interface(
        self,
        up_l_def="3",
        up_r_def="5",
        down_l_def="5",
        down_r_def="3",
        up_ind=0,
        down_ind=1,
    ):
        _, _, col1, col2, _, _ = st.columns([1, 1, 1, 1, 1, 1])
        with col1:
            t_l = st.selectbox(
                "Top left:",
                [up_l_def] + ["─", None],
                index=up_ind,
            )
            b_l = st.selectbox(
                "Bottom left:",
                [down_l_def] + ["─", None],
                index=down_ind,
            )
        with col2:
            t_r = st.selectbox(
                "Top right:",
                [up_r_def] + ["─", None],
                index=up_ind,
            )
            b_r = st.selectbox(
                "Bottom right:",
                [down_r_def] + ["─", None],
                index=down_ind,
            )
        if not t_l:
            t_l = ""
        if not t_r:
            t_r = ""
        if not b_l:
            b_l = ""
        if not b_r:
            b_r = ""
        return t_l, t_r, b_l, b_r
