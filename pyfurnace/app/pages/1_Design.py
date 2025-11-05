import streamlit as st
from streamlit import session_state as st_state
import matplotlib.pyplot as plt

### pyFuRNAce modules
from utils import load_logo, save_origami
import utils.design_functions as des_func
from utils.st_fixed_container import sticky_container


# motif menu function
def motif_menu_expander():
    with st.expander("**Add motifs to the origami:**", expanded=True):
        des_func.make_motif_menu(st.session_state.origami)


if __name__ == "__main__":

    load_logo()

    ### initiate the session state
    des_func.initiate_session_state()

    st.header(
        "Design",
        help="Design your RNA nanostructure and "
        "download it as textfile/python script.",
    )

    ### make the general options for the RNA origami
    des_func.origami_general_options(st.session_state.origami, expanded=False)

    with st.container(
        horizontal=True,
        vertical_alignment="center",
        horizontal_alignment="distribute",
        gap="large",
    ):

        # simple origami popover
        with st.popover(
            "Make a simple origami",
            use_container_width=False,
            icon=":material/mist:",
            help="Start by creating a simple origami rather than "
            "starting from scratch",
        ):
            des_func.simple_origami()

        # motif menu sidebar toggle
        with st.container(horizontal_alignment="center"):
            motif_menu_sidebar = st.toggle(
                "Motif sidebar",
                value=st_state.sidebar_motif_menu,
                help="Show the motif menu in the sidebar "
                "instead of below the general options.",
            )
            if motif_menu_sidebar != st_state.sidebar_motif_menu:
                st_state.sidebar_motif_menu = motif_menu_sidebar
                st.rerun()

            motif_width = 0.01
            if motif_menu_sidebar:
                # make the slider to adjust the sidebar width
                motif_width = st.slider(
                    "Split sidebar width (%)",
                    min_value=10,
                    max_value=90,
                    value=42,
                    format="%d%% sidebar",
                    label_visibility="collapsed",
                    help="Adjust the width of the sidebar "
                    "to better fit the motif menu.",
                )

        # colormap
        cmap = st.selectbox(
            "Colormap:",
            ["Reds", None] + plt.colormaps(),
            key="colormap",
            # label_visibility="collapsed",
            help="Change the color of the OxView visualization.",
        )
        st.session_state.oxview_colormap = cmap

        # gradient toggle
        grad = st.toggle(
            "Color gradient path",
            key="grad",
            help="Toggle the gradient color scheme for the nucleotides",
        )
        st.session_state.gradient = grad

    motif_sidebar, main = st.columns(
        [motif_width / 100, 1 - motif_width / 100],
        border=st_state.get("sidebar_motif_menu", False),
    )

    ### display the motif menu in the sidebar
    if st_state.sidebar_motif_menu:
        with motif_sidebar:
            with sticky_container(mode="top", border=False):
                des_func.make_motif_menu(st.session_state.origami)

    with main:
        # motif menu in the sidebar, so stick the origami view options
        if st.session_state.sidebar_motif_menu:
            with sticky_container(mode="top", border=False):
                view_opt = des_func.origami_select_display()

        # if the motif menu is not in the sidebar, display it here
        else:
            # sticky the motif menu
            if st.session_state.motif_menu_sticky:
                with sticky_container(mode="top", border=False):
                    motif_menu_expander()
                    view_opt = des_func.origami_select_display()
            # normal motif menu
            else:
                motif_menu_expander()
                view_opt = des_func.origami_select_display()

        ### select the render mode
        if not st.session_state.origami:
            st.success("The origami is empty, add a motif!")
            st.stop()
        else:
            des_func.origami_build_view(view_opt)

        ### display dot-bracket structure, sequence constraints
        # and link to the Generate page
        des_func.display_structure_sequence()

        ### Download the RNA origami structure
        save_origami()
