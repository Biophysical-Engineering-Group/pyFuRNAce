import streamlit as st
import matplotlib.pyplot as plt
### import the design functions
from utils import check_import_pyfurnace, load_logo, save_origami
check_import_pyfurnace()
import utils.design_functions as des_func
from utils.st_fixed_container import sticky_container

### set the logo of the app

if __name__ == "__main__":
    load_logo()
    ### initiate the session state
    des_func.initiate_session_state()
    origami = st.session_state.origami
    st.write("## Design RNA origami")

    ### make the general options for the RNA origami
    des_func.origami_general_options(origami, expanded=False)

    cols = st.columns([1] * 7)   
    with cols[0]:    
        with st.popover("Make a simple Origami",
                        help='Start by creating a simple Origami rather than starting from scratch'):
            des_func.simple_origami()
    with cols[3]:
        st.selectbox('OxView 3D ColorMap:', ['Reds', None] + plt.colormaps() , key='oxview_colormap', help='Change the color of the OxView visualization.')
    with cols[-1]:
        st.toggle('Gradient Stran Path', key='gradient', help='Toggle the gradient color scheme for the nucleotides')


    ### option menu to manage the motifs in the origami
    st.write("#### Add Motifs to the Origami:")
    # with st.popover(":green[Add Motifs to the Origami:]",
    #                 use_container_width=True,):
    #     make_motif_menu(origami)
    with sticky_container(mode="top", border=False):
        des_func.make_motif_menu(origami)
        # st.markdown("""<hr style="height:1px;border:none;color:#DDDDDD;background-color:#DDDDDD;" /> """, unsafe_allow_html=True)
        # st.markdown("<hr style='margin-top:+0em;border:none;margin-bottom:-1em;color:#FFFFFFw;background-color:#FFFFFF;' />", unsafe_allow_html=True)
        st.markdown("<hr style='margin-top:+0.5em;margin-bottom:+1.0em;' />", unsafe_allow_html=True)
        view_opt = des_func.origami_select_display()
        # st.divider()


    ### add separator
    # st.markdown("""<hr style="height:1px;border:none;color:#DDDDDD;background-color:#DDDDDD;" /> """, unsafe_allow_html=True)
 
    ### display the RNA origami structure
    # st.write('#### Structure of your RNA origami')
    # st.divider()

    ### select the render mode
    if not st.session_state.origami:
        st.success('The Origami is empty, add a Motif!')
        st.stop()
    else:
        des_func.origami_build_view(view_opt)

    ### display the dot-bracket notation and sequence constraints and link to the Generate page
    des_func.display_structure_sequence()
    ### Download the RNA origami structure
    st.markdown("""<hr style="height:5px;border:none;color:#CCCCCC;background-color:#CCCCCC;" /> """, unsafe_allow_html=True)
    save_origami()

    
    