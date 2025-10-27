import streamlit as st
from streamlit import session_state as st_state
import matplotlib.pyplot as plt

### pyFuRNAce modules
from utils import load_logo, save_origami
import utils.design_functions as des_func
from utils.st_fixed_container import sticky_container


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

    ### make 3 common options for the RNA origami
    cols = st.columns([1, 1, 1, 1], vertical_alignment="center")

    # simple origami popover
    with cols[0]:
        with st.popover(
            "Make a simple origami",
            use_container_width=False,
            help="Start by creating a simple origami rather than "
            "starting from scratch",
        ):
            des_func.simple_origami()

    with cols[1]:
        motif_menu_sidebar = st.toggle(
            "Motif menu in the sidebar",
            value=st_state.sidebar_motif_menu,
            help="Show the motif menu in the sidebar "
            "instead of below the general options.",
        )
        if motif_menu_sidebar != st_state.sidebar_motif_menu:
            st_state.sidebar_motif_menu = motif_menu_sidebar
            st.rerun()

    with cols[2]:
        cmap = st.selectbox(
            "Colormap:",
            ["Reds", None] + plt.colormaps(),
            key="colormap",
            # label_visibility="collapsed",
            help="Change the color of the OxView visualization.",
        )
        st.session_state.oxview_colormap = cmap

    # gradient toggle
    with cols[3]:
        grad = st.toggle(
            "Color gradient path",
            key="grad",
            help="Toggle the gradient color scheme for the nucleotides",
        )
        st.session_state.gradient = grad

    # motif menu function
    def motif_menu_expander():
        with st.expander("**Add motifs to the origami:**", expanded=True):
            des_func.make_motif_menu(st.session_state.origami)

    ### display the motif menu in the sidebar or as sticky container
    if st.session_state.sidebar_motif_menu:
        st.markdown(
            """
            <style>
            [data-testid="stSidebar"][aria-expanded="true"]{
                min-width: 30%;
                max-width: 50%;
            }
            """,
            unsafe_allow_html=True,
        )
        with st.sidebar:
            st.markdown("#### :green[Drag the border of the sidebar to resize it ->]")
            motif_menu_expander()
        # if the motif menu is in the sidebar,
        # keep the origami view menu sticky, so no need to scroll back up
        with sticky_container(mode="top", border=False):
            view_opt = des_func.origami_select_display()
    else:
        if st.session_state.motif_menu_sticky:
            with sticky_container(mode="top", border=False):
                motif_menu_expander()
                view_opt = des_func.origami_select_display()
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
