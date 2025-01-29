import streamlit as st
### import the design functions
from utils import check_import_pyfurnace, load_logo, save_origami
check_import_pyfurnace()
from utils.design_functions import initiate_session_state, simple_origami, origami_general_options, make_motif_menu, origami_display_menu, display_structure_sequence

### set the logo of the app

if __name__ == "__main__":
    load_logo()
    ### initiate the session state
    initiate_session_state()
    origami = st.session_state.origami
    st.write("## Design RNA origami")

    ### make the general options for the RNA origami
    origami_general_options(origami)
    simple_ori =  st.toggle("Make a simple Origami", 
                            key='simple_origami', 
                            help='Start by creating a simple Origami rather than starting from scratch')
    if simple_ori:
        simple_origami()

    ### option menu to manage the motifs in the origami
    if not simple_ori:
        st.write("#### Add Motifs to the Origami:")
        make_motif_menu(origami)
    elif st.session_state.origami:
        st.markdown(':orange[Close] the "Make a simple Origami" toggle :orange[to edit] the origami.')
    ### add separator
    # st.markdown("""<hr style="height:1px;border:none;color:#DDDDDD;background-color:#DDDDDD;" /> """, unsafe_allow_html=True)
 
    ### display the RNA origami structure
    # st.write('#### Structure of your RNA origami')
    st.divider()

    ### select the render mode
    if not st.session_state.origami:
        st.success('The Origami is empty, add a Motif!')
        st.stop()
    else:
        origami_display_menu()
    ### display the dot-bracket notation and sequence constraints and link to the Generate page
    display_structure_sequence()
    ### Download the RNA origami structure
    st.markdown("""<hr style="height:5px;border:none;color:#CCCCCC;background-color:#CCCCCC;" /> """, unsafe_allow_html=True)
    save_origami()

    
    