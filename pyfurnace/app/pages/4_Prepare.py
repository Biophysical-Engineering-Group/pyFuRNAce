import streamlit as st
from streamlit_option_menu import option_menu
import json
import os
import warnings

# App utilities
from utils import (
    load_logo,
    main_menu_style,
    copy_to_clipboard,
    write_format_text,
)
from utils.template_functions import symbols, reference, sanitize_input
from utils.design_functions import origami_build_view
from pyfurnace import prepare as prep

# Melting and primer logic
from pyfurnace.prepare.utils import (
    tm_methods,
    tm_models,
    default_values,
    calculate_gc,
    make_tm_calculator,
    check_dimer,
    auto_design_primers,
    annealing_temp,
)


def upload_setting_button():
    """Allow to upload setting"""
    st.session_state["upload_setting"] = True
    return


def calculate_annealing(seq, mts, c_primer, nc_primer, tm_kwargs):

    col1, col2 = st.columns(2, gap="large")
    with col1:
        subcol1, subcol2 = st.columns([11, 1], vertical_alignment="center")
        with subcol1:
            write_format_text(c_primer)
        with subcol2:
            copy_to_clipboard(c_primer, "")

        subcol1, subcol2, subcol3 = st.columns(3)
        with subcol1:
            st.markdown(f"GC content: {calculate_gc(c_primer)}%")
        with subcol2:
            st.markdown(f"**:green[Tm: {round(mts[0], 1)}°C]**")
        with subcol3:
            st.markdown(f"Length: {len(c_primer)}")

    with col2:

        subcol1, subcol2 = st.columns([11, 1], vertical_alignment="center")
        with subcol1:
            write_format_text(nc_primer)
        with subcol2:
            copy_to_clipboard(nc_primer, "")

        subcol1, subcol2, subcol3 = st.columns(3)
        with subcol1:
            st.markdown(f"GC content: {calculate_gc(nc_primer)}%")
        with subcol2:
            st.markdown(f"**:green[Tm: {round(mts[1], 1)}°C]**")
        with subcol3:
            st.markdown(f"Length: {len(nc_primer)}")

    # show the primers preview, check the self-dimerization of the primer and
    # check the dimer between the primers and the sequence
    c_bases = len(c_primer)
    nc_bases = len(nc_primer)
    with st.expander("Dimerization preview"):
        st.markdown(
            f'<div style="text-align: center;"><span style="color: #D00000">'
            f"{c_primer}</span>{seq[c_bases]}[...]{seq[-nc_bases - 1:]}</div>",
            unsafe_allow_html=True,
        )
        st.markdown(
            f'<div style="text-align: center;">'
            f"{prep.Seq(seq[:c_bases + 1]).complement()}"
            f"[...]"
            f"{prep.Seq(seq[-nc_bases - 1]).complement()}"
            f'<span style="color: #D00000">'
            f"{nc_primer[::-1]}</span></div>",
            unsafe_allow_html=True,
        )
        st.divider()
        col1, col2 = st.columns(2, gap="large")
        with col1:
            self_dimer1 = check_dimer(
                c_primer, c_primer, basepair={"A": "T", "T": "A", "C": "G", "G": "C"}
            )
            write_format_text("Self-Dimer forward primer\n" + self_dimer1)

        with col2:
            self_dimer2 = check_dimer(
                nc_primer, nc_primer, basepair={"A": "T", "T": "A", "C": "G", "G": "C"}
            )
            write_format_text("Self-Dimer reverse primer\n" + self_dimer2)

        self_dimer12 = check_dimer(
            c_primer, nc_primer, basepair={"A": "T", "T": "A", "C": "G", "G": "C"}
        )
        write_format_text("Dimer between primers\n" + self_dimer12)

    # add a warning if the melting temperature of the primers is too different
    anneal_color = "green"
    if abs(mts[0] - mts[1]) >= 5:
        st.warning("The difference of Tm should be below 5°C", icon="⚠️")
        anneal_color = "red"

    # take into account the two main method to calculate the PCR annealing
    # temperature: IDT method and Phusion Buffer method
    t_anneal = annealing_temp(
        mts, seq, tm_kwargs, method=st.session_state.annealing_model.split()[0]
    )
    if t_anneal < 50:
        st.warning(
            "Is suggested a temperature of annealing higher than 50°C.", icon="⚠️"
        )
        anneal_color = "red"

    # write the annealing temperature in big
    st.markdown(f"### Anneal at: :{anneal_color}[{t_anneal}°C]")


def primers_tab(seq, mt_correct):
    """Calculate the melting temperature of the primers and check
    the dimer between the primers and the sequence"""

    # show the settings for the two primers:
    # choose the number of bases and show the gc content and melting temperature
    st.write("\n")
    col1, col2 = st.columns(2, gap="large")
    mts = [0, 0]
    # settings for the coding primer
    with col1:
        with st.columns(5)[2]:
            st.markdown("###### Forward")

        c_bases = st.slider(
            "Primer length:",
            min_value=5,
            max_value=50,
            value=st.session_state.prim_fwd_len,
            key=f"coding_primer_{st.session_state.auto_count}",
        )
        c_primer = seq[:c_bases]

        mts[0] = round(mt_correct(c_primer), 2)

    # settings for the non-coding primer
    with col2:
        with st.columns(5)[2]:
            st.markdown("###### Reverse")

        nc_bases = st.slider(
            "Primer length:",
            min_value=5,
            max_value=50,
            value=st.session_state.prim_rev_len,
            key=f"non_coding_primer_{st.session_state.auto_count}",
        )

        nc_primer = str(prep.Seq(seq[-nc_bases:]).reverse_complement())

        mts[1] = round(mt_correct(nc_primer), 2)

    return mts, c_primer, nc_primer


def tm_parameters():
    # load settings if you want to upload custom settings
    load_settings = st.checkbox(
        "Upload previous energy parameter:",
        help="If you have saved a json file with the "
        "energy parameters, you can upload it and"
        " use the same settings.",
    )
    if load_settings:
        saved_setting = st.file_uploader(
            "Upload previous energy parameter " "(optional):",
            on_change=upload_setting_button,
            type="json",
            help="If you have saved a json file "
            "with the energy parameters, you "
            "can upload it and use the same "
            "settings.",
        )

        # if load the settings, upload the valuse from the json file
        if saved_setting and st.session_state["upload_setting"]:
            # To read file as bytes:
            session_old = json.load(saved_setting)
            for key, value in session_old.items():
                st.session_state[key] = value
            st.session_state["upload_setting"] = False
            st.rerun()

    # initialize session state values for the energy model
    for key, val in default_values.items():
        if key not in st.session_state:
            st.session_state[key] = val

    # select the energy model
    st.write("Energy model:")
    col1, col2, col3 = st.columns([2, 2, 1])

    # select the melting temperature method
    with col1:
        mt_type = st.selectbox(
            "Melting Temperature Method",
            list(tm_methods),
            help=(prep.mt.__doc__),
            key="mt_method",
            label_visibility="collapsed",
        )
        method = tm_methods[mt_type]

    # select the melting temperature model
    with col2:
        mt_model = st.selectbox(
            "Energy Correction",
            list(tm_models[mt_type]),
            key="mt_model",
            label_visibility="collapsed",
        )

    # add a button to show the help of the selected method
    with col3:
        help_model = st.button("Show info", key="energy_model_info")
    if help_model:
        st.help(method)

    # buffer correction parameters
    st.divider()
    if method != list(tm_methods)[2]:

        st.write("Buffer corrections (mM):")

        cols = st.columns(8)
        with cols[0]:
            na = st.number_input("Na", key="Na")
        with cols[1]:
            k = st.number_input("K", key="K")
        with cols[2]:
            tris = st.number_input("Tris", key="Tris")
        with cols[3]:
            mg = st.number_input("Mg", value=st.session_state["Mg"], key="Mg")
        with cols[4]:
            dntps = st.number_input(
                "dNTPs", value=st.session_state["dNTPs"], key="dNTPs"
            )
        with cols[5]:
            dmso = st.number_input(
                "DMSO (%)", value=st.session_state["DMSO (%)"], key="DMSO (%)"
            )
        with cols[6]:
            correction_met = st.selectbox("Method", [0, 1, 2, 3, 4, 6, 7], key="Method")
        with cols[7]:
            help_correction = st.button("Show info", key="energy_correction_info")

        tm_kwargs = {
            "Na": na,
            "K": k,
            "Tris": tris,
            "Mg": mg,
            "saltcorr": correction_met,
            "dNTPs": dntps,
            "DMSO": dmso,
        }

    if help_correction:
        # add a button to show the help of the selected method
        st.help(prep.mt.salt_correction)
        st.help(prep.mt.chem_correction)

    primer_conc = st.number_input("Primer conc. (nM)", key="Primer")

    # create the function to calculate the TM
    calculate_mt = make_tm_calculator(
        method_name=mt_type,
        model_name=mt_model,
        primer_conc=primer_conc,
        tm_kwargs=tm_kwargs,
    )

    # save settings and allow download
    with open("energy_parameters.json", "w") as settings:
        # cannot save session state as it is,
        # I have to convert it to a dictionary
        session_dict = {key: st.session_state[key] for key in default_values}
        json.dump(session_dict, settings)
    with open("energy_parameters.json", "rb") as file:
        st.download_button(
            label="Download energy parameters",
            data=file,
            file_name="energy_parameters.json",
            mime="application/json",
            help="Save the current settings (e.g. ions concentration, "
            "Tm model), so you can easily reload them if you refresh "
            "the page!",
        )

    return calculate_mt, tm_kwargs


def primers_setup():
    if "dna_template" not in st.session_state:
        st.session_state["dna_template"] = ""

    # take the input sequence and sanitize it
    seq = sanitize_input(
        st.text_area("Input sequence:", value=st.session_state["dna_template"])
    )

    # check the symbols in the sequence
    if set(seq) - symbols:
        st.warning(
            "The sequence contains symbols not included in the "
            "[IUPAC alphabet](https://www.bioinformatics.org/sms/iupac.html).",
            icon="⚠️",
        )
    if "U" in seq:
        st.error("The DNA template contains U", icon=":material/personal_injury:")

    if not seq:
        st.stop()
    elif len(seq) < 40:
        st.error("The sequence is too short", icon=":material/personal_injury:")
        st.stop()

    ###
    # Various for PCR settings
    ###

    mcol1, mcol2, mcol3, mcol4 = st.columns(4, vertical_alignment="center", gap="large")
    with mcol1:
        with st.popover("Melting temperature parameters", use_container_width=True):

            calculate_mt, tm_kwargs = tm_parameters()

    with mcol2:
        with st.popover(
            "Add primers overhangs",
            use_container_width=True,
            help="Add restriction sites or other sequences to the primers"
            " without affecting the gene sequence. The overhangs will be "
            "added to the 5' end of the primers and taken into account "
            "for the melting temperature calculation.",
        ):
            restric_site_1 = st.text_input(
                "Coding overhang sequence (5'->3')" " before the promoter:"
            )
            restric_site_2 = st.text_input(
                "Non-coding overhang sequence (5'->3')" " after the fragment:"
            )

    with mcol3:
        with st.popover("Annealing calculation", use_container_width=True):
            st.selectbox(
                "Calculate annealing:",
                ["IDT method [2]", "Phusion method [3]"],
                key="annealing_model",
                label_visibility="collapsed",
                help="Select the method to calculate the annealing "
                "temperature. Check the references for more information.",
            )
    with mcol4:
        # default primers length
        with st.popover(
            "**Auto design primers**",
            icon=":material/precision_manufacturing:",
            use_container_width=True,
        ):

            target_temp = st.number_input(
                "Target melting temperature (°C)",
                value=65,
                min_value=0,
                max_value=100,
                step=1,
            )
            if st.button(
                "Auto design primers",
                key="auto_design_primers",
                help="Automatically design the primers for the sequence. "
                "The primers are designed to have a melting temperature "
                "and follow best practices (see documentation).",
            ):

                status = st.empty()
                status.info(
                    "Designing primers...", icon=":material/precision_manufacturing:"
                )

                primers, final_mts = auto_design_primers(seq, target_temp, calculate_mt)

                st.session_state.auto_count += 1

                # set the new primers length in the session state
                if st.session_state.prim_fwd_len != len(
                    primers[0]
                ) or st.session_state.prim_rev_len != len(primers[1]):
                    st.session_state.prim_fwd_len = len(primers[0])
                    st.session_state.prim_rev_len = len(primers[1])

                status.success(
                    "Primer designed!", icon=":material/precision_manufacturing:"
                )

    seq = restric_site_1 + seq + str(prep.Seq(restric_site_2).reverse_complement())

    mts, c_primer, nc_primer = primers_tab(seq, calculate_mt)

    calculate_annealing(seq, mts, c_primer, nc_primer, tm_kwargs)

    # add bibliography
    reference(True)


def md_setup():

    # check if the origami blueprint is good
    if not st.session_state.get("origami") or not all(
        n in "AUCG&" for n in st.session_state.origami.sequence
    ):

        st.error(
            "The origami blueprint is empty or it doesn't contain a sequence",
            icon=":material/personal_injury:",
        )

        col1, col2 = st.columns(2)
        with col1:
            st.page_link(
                "pages/1_Design.py",
                label="**:orange[Have a look to the blueprint]**",
                icon=":material/draw:",
            )

        with col2:
            st.page_link(
                "pages/2_Generate.py",
                label="**:orange[Generate the sequence]**",
                icon=":material/network_node:",
            )
        st.stop()

    # show the origami blueprint
    with st.expander("**Origami preview:**", expanded=True):
        origami_build_view("Origami split view")

    st.divider()
    col1, col2 = st.columns(2)
    with col1:
        oxdna_dir = st.text_input(
            "OxDNA directory " '(e.g. "~/Documents/software/oxDNA")'
        )
    with col2:
        temp = st.number_input(
            "Temperature (°C)", value=37, min_value=0, max_value=100, step=1
        )

    help_mc_rel = (
        "Use Monte Carlo relaxation to relax the origami before the MD "
        "simulation, using external forces to keep the basepairing."
    )
    help_md_rel = (
        "Use Molecular Dynamics to relax the structure, using forces to"
        "maintain the basepairing"
    )
    help_md_equ = (
        "An extra Molecular Dynamics simulation to make sure the Origami"
        "is fully relaxed and not bias the starting configuration of the "
        "production. **External forces for pseudoknots only** are used "
        "in this simulation, to make sure the pseudoknots are paired "
        "at the start of the production"
    )
    help_md_run = (
        "The final Molecular Dynamics simulation, to simulate the "
        "structure without any external force."
    )

    num_input = st.toggle("Use numerical input", key="num_input")

    col1, col2, col3, col4 = st.columns(4)
    with col1:
        if num_input:
            step_mc_rel = st.number_input(
                "MC relaxation steps",
                min_value=0,
                value=int(5e3),
                help=help_mc_rel,
            )
        else:
            step_mc_rel = st.slider(
                "MC relaxation steps",
                min_value=0,
                value=int(5e3),
                max_value=int(10e3),
                step=int(1e3),
                help=help_mc_rel,
                format="%0.0e",
            )
    with col2:
        if num_input:
            step_md_rel = st.number_input(
                "MD relaxation steps",
                min_value=0,
                value=int(1e7),
                help=help_md_rel,
            )
        else:
            step_md_rel = st.slider(
                "MD relaxation steps",
                min_value=0,
                value=int(1e7),
                max_value=int(1e8),
                step=int(1e6),
                help=help_md_rel,
                format="%0.0e",
            )
    with col3:
        if num_input:
            step_md_equ = st.number_input(
                "MD equilibration steps",
                min_value=0,
                value=int(1e8),
                help=help_md_equ,
            )
        else:
            step_md_equ = st.slider(
                "MD equilibration steps",
                min_value=0,
                value=int(1e8),
                max_value=int(1e9),
                step=int(1e7),
                help=help_md_equ,
                format="%0.0e",
            )
    with col4:
        if num_input:
            step_md_run = st.number_input(
                "MD production steps",
                min_value=0,
                value=int(1e9),
                help=help_md_run,
            )
        else:
            step_md_run = st.slider(
                "MD production steps",
                min_value=0,
                value=int(1e9),
                max_value=int(1e10),
                step=int(1e8),
                help=help_md_run,
                format="%0.0e",
            )

    if not oxdna_dir or not temp:
        st.stop()

    zip_path = prep.oxdna_simulations(
        origami=st.session_state.origami,
        oxdna_directory=oxdna_dir,
        temperature=temp,
        mc_relax_steps=step_mc_rel,
        md_relax_steps=step_md_rel,
        md_equil_steps=step_md_equ,
        md_prod_steps=step_md_run,
    )
    # add a button to download the zip file
    st.divider()

    col1, col2 = st.columns(2, vertical_alignment="bottom")
    with col1:
        name = st.text_input("Name of the simulation:", value="origami_simulation")
    with col2:
        with open(zip_path, "rb") as file:
            st.download_button(
                label="Download MD simulation files",
                data=file,
                file_name=f"{name}.zip",
                mime="application/zip",
                help="Download the MD simulation files. "
                "The zip file contains the input files for the MD simulation.",
            )

    # delete the zip file
    if os.path.exists(zip_path):
        os.remove(zip_path)
    return


if __name__ == "__main__":
    ### set the logo of the app
    load_logo()
    warnings.filterwarnings("ignore")  # ignore warnings

    if "prepare_ind" not in st.session_state:
        st.session_state.prepare_ind = 0
    if "prim_fwd_len" not in st.session_state:
        st.session_state.prim_fwd_len = 21
    if "prim_rev_len" not in st.session_state:
        st.session_state.prim_rev_len = 21
    if "auto_count" not in st.session_state:
        st.session_state.auto_count = 0

    # create the tabs with the functions
    st.header(
        "Prepare",
        help="Design primers for your DNA template "
        "or prepare the Origami for OxDNA simulation.",
    )
    option_data = {"Primers": "bi bi-arrow-left-right", "MD simulations": "bi bi-cpu"}

    selected_operation = option_menu(
        None,
        list(option_data.keys()),
        icons=list(option_data.values()),
        menu_icon="cast",
        orientation="horizontal",
        default_index=st.session_state.prepare_ind,
        styles=main_menu_style,
    )

    if selected_operation == "MD simulations":
        if st.session_state.prepare_ind == 0:
            st.session_state.prepare_ind = 1
            st.rerun()
        md_setup()
    elif selected_operation == "Primers":
        if st.session_state.prepare_ind == 1:
            st.session_state.prepare_ind = 0
            st.rerun()
        primers_setup()
