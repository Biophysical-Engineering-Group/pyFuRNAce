import streamlit as st
from utils import load_logo, app_path

if __name__ == '__main__':
    load_logo() 

    st.write("# Hello and Welcome to pyFuRNAce!")

    st.write('Design and generate RNA nanostructures in few simple steps.')

    st.page_link("pages/1_Design.py", 
                 label=":orange[Design:]", 
                 icon=":material/draw:")

    st.markdown("- Design your RNA nanostructure and download it as "
                "textfile/python script.")

    st.page_link("pages/2_Generate.py", 
                 label=":orange[Generate:]", 
                 icon=":material/network_node:")

    st.markdown("- Generate the RNA sequence that matches the desired dot-bracket" 
                " notation for the nanostructure.")

    st.page_link("pages/3_Convert.py", 
                 label=":orange[Template:]", 
                 icon=":material/genetics:")

    st.markdown("- Prepare the DNA template for you RNA Origami, search subsequences "
                "and search for dimers.")

    st.page_link("pages/4_Prepare.py", 
                 label=":orange[Prepare:]", 
                 icon=":material/sync_alt:")

    st.markdown("- Design primers for your DNA template or prepare the Origami for "
                "OxDNA simulation.")
    
    st.write("### Check out the 1-min walkthrough video:")
    st.video(str(app_path / "static" / "walkthrough_1min.mp4"), 
             format="video/mp4", 
             start_time=0, 
             subtitles=None, 
             loop=True, 
             autoplay=True, 
             muted=True)