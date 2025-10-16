import inspect
import streamlit as st
from streamlit import session_state as st_state

### pyFuRNAce modules
from utils import load_logo
from utils.template_functions import sanitize_input
from utils.design_functions import update_code, initiate_session_state
import pyfurnace as pf
import pyfurnace.design.utils.origami_lib as origami_lib


def load_origami_button():
    st_state.new_comer = False
    if st.button("Load the Origami Design", icon=":material/draw:", type="primary"):

        st_state.code.append(st_state.qs_text)
        st_state.origami = st_state.qs_origami
        st.switch_page("pages/1_Design.py")


@st.dialog("Load a sequence/dot-bracket")
def load_seq_dot():
    st.write(
        "You can start by adding either your sequence, "
        "your dot-bracket structure or both."
    )
    col1, col2 = st.columns(2)
    with col1:
        sequence = st.text_input(
            "Sequence:",
            key="Sequence",
            help="To generate an RNA origami from a sequence,"
            " we fold it with viennaRNA and then convert "
            "the dot bracket notation to an Origami, "
            "as an collection of stems and bulges.",
        )

    with col2:
        structure = st.text_input(
            "Dot-bracket:",
            key="Dot-bracket_structure",
            help="Add a dot-bracket structure to generate "
            " a blueprint. If the secondary structure "
            " of known aptamers is detected, the aptamer "
            " motif will be used in the Origami. ",
        )

    structure = sanitize_input(structure)
    sequence = sanitize_input(sequence)

    if not (structure or sequence):
        return

    if structure and sequence:
        col = st.columns(1)[0]
    else:
        col = col2 if structure else col1

    with col:
        try:
            # try on the fly to see if it would work
            st_state.qs_origami = pf.Origami.from_structure(
                structure=structure, sequence=sequence
            )
            st_state.qs_text = (
                f"origami = pf.Origami.from_structure("
                f'structure="{structure}", '
                f'sequence="{sequence}")'
            )
        except ValueError as e:
            st.error(f"Error: {e}", icon=":material/personal_injury:")
        else:
            st.success("Origami created successfully!", icon=":material/check_circle:")
    load_origami_button()


@st.dialog("Load a template")
def load_template():
    templ_funcs = {
        name.replace("template_", ""): func
        for name, func in vars(origami_lib).items()
        if callable(func) and name.startswith("template_")
    }

    template_docs = "Available templates:\n\n"
    for name, func in templ_funcs.items():
        template_docs += f"**--> {name}** \n\n "
        docs = inspect.getdoc(func)
        template_docs += f'{docs.split("Returns")[0].strip()}\n\n'

    template = st.selectbox(
        "Select a template:",
        options=list(templ_funcs.keys()),
        index=None,
        help=template_docs,
    )
    if not template:
        return

    func_code = inspect.getsource(templ_funcs[template])
    # remove one layer of indentation
    only_code = func_code.split('"""')[2].split("return")[0]
    func_code = "\n".join([line[4:] for line in only_code.splitlines()])
    # remove the pyfurnace import if present
    func_code = func_code.replace("import pyfurnace as pf", "").strip()

    st_state.qs_origami = templ_funcs[template]()
    st_state.qs_text = f"{func_code.strip()}\n"

    st.success(f"{template} loaded successfully!", icon=":material/check_circle:")

    load_origami_button()


@st.dialog("Load a file")
def load_file():
    uploaded_file = st.file_uploader(
        "Add a fasta file (.fasta) with sequence and "
        " structure or a python script (.py)"
        " with the pyfurnace origami code.",
        type=["py", "fasta"],
    )

    if uploaded_file is None:
        return

    if ".fasta" in uploaded_file.name:
        file_contents = uploaded_file.read().decode("utf-8")
        dot_b_symbols = set(pf.all_pk_symbols).union({"(", ")", "."})
        try:
            lines = file_contents.strip().split("\n")
            structure = ""
            for i, line in enumerate(lines):
                if line.startswith(">"):
                    st_state.qs_text = f"# Loading sequence: {line[1:].strip()}\n"
                    sequence = lines[i + 1].strip()
                    if (
                        i + 2 < len(lines)
                        and not lines[i + 2].startswith(">")
                        and all(sym in dot_b_symbols for sym in lines[i + 2])
                    ):
                        structure = lines[i + 2].strip()
                    break

            st_state.qs_origami = pf.Origami.from_structure(
                structure=structure, sequence=sequence
            )
            st_state.qs_text = (
                f"origami = pf.Origami.from_structure("
                f'structure="{structure}", '
                f'sequence="{sequence}")'
            )

        except Exception as e:
            st.error(
                f"Error loading FASTA file: {e}", icon=":material/personal_injury:"
            )

    elif ".py" in uploaded_file.name:
        file_contents = uploaded_file.read().decode("utf-8")
        st_state.qs_text = file_contents.strip()

        with st.spinner("Loading the Origami..."):
            origami = update_code(st_state.qs_text, return_origami=True)
        if not origami:
            return
        st_state.qs_origami = origami
        st.success("File loaded successfully!", icon=":material/check_circle:")

    load_origami_button()


def newcomer_flow():
    cols_ratio = [1.7, 2.1, 2, 4]
    col1, col2, col3, _ = st.columns(cols_ratio, vertical_alignment="bottom")
    sequence_type = "secondary"
    template_type = "secondary"
    file_type = "secondary"

    with col1:
        if st.button("New to pyFuRNAce?", icon=":material/new_window:", type="primary"):
            st_state.new_comer = True

    if st_state.new_comer:
        question_col = st.columns(3)
        with question_col[0]:
            start = st.radio(
                "Do you have a sequence or structure (dot-bracket)?",
                index=None,
                options=["No", "Yes"],
                horizontal=True,
                key="design_choice",
            )

        if start == "Yes":
            with question_col[1]:
                start2 = st.radio(
                    "-> Do you have a fasta or python file?",
                    index=None,
                    options=["No", "Yes"],
                    horizontal=True,
                    key="design_choice2",
                )

            col1, col2, col3, _ = st.columns(cols_ratio, vertical_alignment="bottom")

            if start2 == "Yes":
                file_type = "primary"
                with col3:
                    st.write("Try this one ↓")

            elif start2 == "No":
                sequence_type = "primary"
                with col2:
                    st.write("Add your sequence or structure: ↓")

        elif start == "No":
            col1, col2, col3, _ = st.columns(cols_ratio, vertical_alignment="bottom")
            template_type = "primary"
            with col1:
                st.write("Try this one ↓")

    col1, col2, col3, _ = st.columns(cols_ratio, vertical_alignment="bottom")

    with col1:
        if st.button(
            "Load a template", icon=":material/wallpaper:", type=template_type
        ):
            load_template()

    with col2:
        if st.button(
            "Load sequence/structure", icon=":material/stylus_note:", type=sequence_type
        ):
            load_seq_dot()

    with col3:
        if st.button("Load a file", icon=":material/upload:", type=file_type):
            load_file()


if __name__ == "__main__":
    load_logo()
    st_state.qs_text = None
    st_state.qs_origami = None
    st_state.setdefault("new_comer", False)
    initiate_session_state()

    st.write("# Hello and Welcome to pyFuRNAce!")

    st.write("Design and generate RNA nanostructures in few simple steps.")

    st.page_link("pages/1_Design.py", label=":orange[Design:]", icon=":material/draw:")

    st.markdown(
        "- Design your RNA nanostructure and download it as " "textfile/python script."
    )

    st.page_link(
        "pages/2_Generate.py",
        label=":orange[Generate:]",
        icon=":material/network_node:",
    )

    st.markdown(
        "- Generate the RNA sequence that matches the desired dot-bracket"
        " notation for the nanostructure."
    )

    st.page_link(
        "pages/3_Convert.py", label=":orange[Convert:]", icon=":material/genetics:"
    )

    st.markdown(
        "- Prepare the DNA template for you RNA Origami, search subsequences"
        " and search for dimers."
    )

    st.page_link(
        "pages/4_Prepare.py", label=":orange[Prepare:]", icon=":material/sync_alt:"
    )

    st.markdown(
        "- Design primers for your DNA template or prepare the Origami for "
        "OxDNA simulation."
    )

    newcomer_flow()

    st.divider()

    st.write("### About pyFuRNAce")

    st.markdown(
        "pyFuRNAce is an open-source Python package and web-based design "
        "engine for creating complex RNA nanostructures using the "
        " co-transcriptional RNA origami approach."
    )
    st.markdown(
        " - **GitHub**: [Biophysical-Engineering-Group/pyFuRNAce]"
        "(https://github.com/Biophysical-Engineering-Group/pyFuRNAce)"
    )
    st.markdown(" - **PyPI**: [pyfurnace](https://pypi.org/project/pyfurnace/)")
    st.markdown(
        " - **Documentation**: [Read the Docs]"
        "(https://pyfurnace.readthedocs.io/en/latest/)"
    )
    st.markdown(
        " - bug reports, feature requests or any other questions, "
        "please reach out to us via the "
        "[GitHub Issues](https://github.com/Biophysical-Engineering"
        "-Group/pyFuRNAce/issues)"
        " or the "
        "[GitHub Discussions](https://github.com/Biophysical-Engineering"
        "-Group/pyFuRNAce/discussions)."
    )

    st.markdown("If you use pyFuRNAce in your research, please cite:")
    st.markdown(
        "Monari, L., Braun, I., Poppleton, E. & Göpfrich, K. PyFuRNAce: "
        "An integrated design engine for RNA origami (2025) "
        "[doi:10.1101/2025.04.17.647389]"
        "(https://doi.org/10.1101/2025.04.17.647389)."
    )

    st.divider()

    st.write("#### Check out the 1-min demo video:")
    st.video(
        "https://github.com/Biophysical-Engineering-Group/pyFuRNAce/blob"
        "/main/vid/demo_1min.mp4?raw=true",
        format="video/mp4",
        start_time=0,
        subtitles=None,
        loop=True,
        autoplay=True,
        muted=True,
    )
