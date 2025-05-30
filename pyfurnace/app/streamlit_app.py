import streamlit as st

if __name__ == '__main__':

    pg = st.navigation([
                        st.Page("pages/0_Home.py", title="Home", icon=":material/home:", url_path="Home"),
                        st.Page("pages/1_Design.py", title="Design", icon=":material/draw:", url_path="/design"),
                        st.Page("pages/2_Generate.py", title="Generate", icon=":material/network_node:", url_path="/generate"),
                        st.Page("pages/3_Convert.py", title="Template", icon=":material/genetics:", url_path="/template"),
                        st.Page("pages/4_Prepare.py", title="Prepare", icon=":material/sync_alt:", url_path="/prepare")])

    pg.run()