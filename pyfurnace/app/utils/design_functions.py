import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import html
from time import sleep
import warnings

# Streamlit imports
import streamlit as st
from streamlit import session_state as st_state
from streamlit_option_menu import option_menu
from st_click_detector import click_detector
from st_oxview import oxview_from_text
from code_editor import code_editor

### pyFuRNAce modules
from . import main_menu_style, second_menu_style, copy_to_clipboard
import pyfurnace as pf
from utils.commands import (  # noqa: F401
    AptamersCommand,
    ConnectionsCommand,
    start_end_stemCommand,
    DovetailCommand,
    GeneralEditCommand,
    KissingLoopsCommand,
    StemCommand,
    StructuralCommand,
    TetraLoopCommand,
)

code_editor_buttons = (
    {
        "name": "Copy",
        "feather": "Copy",
        "hasText": True,
        "alwaysOn": True,
        "commands": ["copyAll"],
        "style": {"top": "0.46rem", "right": "0.4rem"},
    },
    {
        "name": "Run",
        "feather": "Play",
        "primary": True,
        "hasText": True,
        "showWithIcon": True,
        "commands": ["submit"],
        "style": {"bottom": "0.44rem", "right": "0.4rem"},
    },
)

funny_bootstrap_icons = [
    "robot",
    "trash",
    "umbrella",
    "camera",
    "cart",
    "cpu",
    "cup-straw",
    "trophy",
    "palette",
    "cup-straw",
    "camera-reels",
    "puzzle",
    "hourglass-split",
    "mortarboard",
]

direction_list = pf.Direction.names()


def origami_general_options(origami, expanded=True):
    with st.expander("General settings", expanded=expanded):
        cols = st.columns(3)

        ### select alignment
        with cols[0]:
            align = st.radio(
                "Alignment type:",
                ["To the left", "Junctions alignment: first", "Line center"],
                index=0,
                key="alignment",
            )

            ### add the alignment to the code
            if st.button("Set alignment", key="set_alignment"):
                align = align.split()[-1]
                origami.align = align
                st_state.code.append(f"origami.align = '{align}' # Align the lines")

        with cols[1]:
            st.toggle(
                "Optimize the blueprint for ROAD",
                value=False,
                key="to_road",
                help="Optimize the blueprint for the ROAD software. "
                'This substitutes the Kissing Loops base pairings with "*";'
                ' and the short stem base pairings with "!".',
            )

        with cols[2]:
            new_sticky = st.toggle(
                "Stick the motif menu at the top",
                value=st_state.motif_menu_sticky,
                key="sticky_menu",
                help="Keep the motif menu and origami visualization"
                " menu to stick to the top of the page.",
            )
            if new_sticky != st_state.motif_menu_sticky:
                st_state.motif_menu_sticky = new_sticky
                st.rerun()

        col1, col2 = st.columns(2)
        with col1:
            font_size = st.slider(
                "Origami font size",
                min_value=2,
                max_value=50,
                value=14,
                key="ori_font_size",
            )
            st_state.origami_font_size = font_size

        with col2:
            frame_size = st.slider(
                "Oxview frame size (disable and renable "
                "the OxView to apply changes)",
                min_value=0,
                max_value=2000,
                value=500,
                key="oxFrame_size",
            )
            st_state.oxview_frame_size = frame_size


def simple_origami():
    with st.form(key="simple_origami_form"):
        dt_text = st.text_input(
            "Enter a list of the angles between helices, " "separated by commas",
            "120, ",
            key="dt_text_list",
            help=f"The program calculates the best connections"
            f" bewteen the helices to fit the given angles."
            f"The connection between helices (Dovetails) are "
            f"obtained roughly with this lookup table "
            f"(angle --> dt): {pf.ANGLES_DT_DICT}",
        )

        angle_list = [int(x) for x in dt_text.split(",") if x and x.strip()]
        dt_list = pf.convert_angles_to_dt(angle_list)
        main_stem_default = 11 * (
            (max([abs(dt) for dt in dt_list], default=0) + 17) // 11 + 1
        )
        col1, col2 = st.columns(2, vertical_alignment="bottom")
        with col1:
            kl_columns = st.number_input(
                "Kissing loop columns:",
                min_value=1,
                value=1,
                help="number of KL repeats in the helix",
            )
        with col2:
            main_stem = st.number_input(
                "Spacing between crossovers (bp):",
                min_value=22,
                value=main_stem_default,
                step=11,
                help="The length of the consecutive stems " "in the helix",
            )

        submitted = st.form_submit_button("Submit")
        if submitted:
            st_state.origami = pf.simple_origami(
                dt_list=angle_list,
                kl_columns=kl_columns,
                main_stem=main_stem,
                add_terminal_helix=True,
                align=st_state.origami.align,
                use_angles=True,
            )
            st_state.code.append(
                f"origami = pf.simple_origami(dt_list={angle_list}, "
                f"kl_columns={kl_columns}, "
                f"main_stem={main_stem}, "
                f"add_terminal_helix=True, "
                f'align="{st_state.origami.align}", use_angles=True) '
                f"# Create a simple origami"
            )

            # select the end of the origami
            st_state.line_index = len(st_state.origami) - 1
            st_state.motif_index = len(st_state.origami[-1])


def motif_text_format(motif: pf.Motif) -> str:
    # take a motif string and expand it on right, left, top and bottom with spaces
    motif_list = (
        [[" "] * (motif.num_char + 2)]
        + [[" "] + [char for char in line] + [" "] for line in str(motif).split("\n")]
        + [[" "] * (motif.num_char + 2)]
    )

    ### Add the 5' and 3' to the motif text
    for s in motif:

        if not s.sequence:  # skip strand without sequence
            continue

        y0, x0 = s.prec_pos[1] + 1, s.prec_pos[0] + 1
        if s[0] not in "35" and motif_list[y0][x0] == " ":
            motif_list[y0][x0] = "<" + s.directionality[0] + ">"

        y1, x1 = s.next_pos[1] + 1, s.next_pos[0] + 1
        if s[-1] not in "35" and motif_list[y1][x1] == " ":
            motif_list[y1][x1] = "<" + s.directionality[1] + ">"

    # remove the top and bottom lines if they are empty
    for i in (0, -1):
        if all([char == " " for char in motif_list[i]]):
            motif_list.pop(i)

    motif_str = "\n".join(["".join(line) for line in motif_list])

    preview_txt = (
        motif_str.replace(" ", "&nbsp;")
        .replace("\n", "<br />")
        .replace("<5>", '<span style="color: #D52919;">5</span>')
        .replace("<3>", '<span style="color: #D52919;">3</span>')
    )
    return preview_txt


def initiate_session_state():
    default_state = {
        "origami": pf.Origami(),
        "code": [
            "import pyfurnace as pf",
            "origami = pf.Origami()",
            "RENDER_TARGET = origami",
        ],
        "motif_buffer": "",
        "mod_motif_buffer": "",
        "motif": pf.Motif(),
        #   "redo": [],
        "motif_menu_sticky": False,
        "copied_motif": None,
        "copied_motif_text": "",
        "modified_motif_text": "",
        "upload_key": 0,
        # Custom motif parameters:
        "custom_strands": [],
        "custom_motifs": [("Custom Motif", "dpad", pf.Motif(lock_coords=False))],
        "custom_key": 0,
        "custom_edit": True,
        "last_clicked_position": "",
        "line_index": 0,
        "motif_index": 0,
        "current_line_occupied": False,
        "ori_click_count": 0,
        "max_pk_index": 1,
        "oxview_selected": (),
        "selected_motif": "Connections",
        "gradient": False,
        "oxview_colormap": "Reds",
        "origami_font_size": 14,
        "oxview_frame_size": 500,
    }
    for key, value in default_state.items():
        if key not in st_state:
            st_state[key] = value


def update_file_uploader():
    st_state.upload_key += 1


@st.fragment
def make_motif_menu(origami):
    # update max pk index
    pk_motifs = origami[lambda m: hasattr(m, "pk_index")]
    st_state.max_pk_index = (
        max(
            [abs(int(x.pk_index.replace("'", ""))) for line in pk_motifs for x in line],
            default=0,
        )
        + 1
    )

    option_data = {
        "Connections": "bi-sliders",
        "Structural": "bi-bricks",
        "Kissing Loops": "bi-heart-half",
        "Aptamers": "bi-palette",
        "Custom": "bi-joystick",
        "Edit": "bi-bandaid",
        "Code": "bi-code",
    }

    select_ind = 0
    if st_state.selected_motif in option_data:
        select_ind = list(option_data).index(st_state.selected_motif)

    selected_motif = option_menu(
        None,
        list(option_data.keys()),
        icons=list(option_data.values()),
        menu_icon="cast",
        orientation="horizontal",
        manual_select=select_ind,
        styles=main_menu_style,
        key="motif_menu",
    )

    if selected_motif != st_state.selected_motif:
        st_state.selected_motif = selected_motif
        try:
            st.rerun(scope="fragment")
        except Exception as e:
            print(f"Error in rerun: {e}")
            st.rerun()

    motif_add = True

    if selected_motif == "Connections":
        ConnectionsCommand().execute()
    elif selected_motif == "Structural":
        StructuralCommand().execute()
    elif selected_motif == "Kissing Loops":
        KissingLoopsCommand().execute()
    elif selected_motif == "Aptamers":
        AptamersCommand().execute()
    elif selected_motif == "Custom":

        ### Adding the menu here, otherwise there are issues choosing the
        # custom motif with st.fragment when the edit mode is off
        col1, col2 = st.columns([5, 1.5], vertical_alignment="bottom")
        with col1:
            motif_selected = option_menu(
                None,
                [nam_ico[0] for nam_ico in st_state.custom_motifs],
                icons=[f"bi bi-{nam_ico[1]}" for nam_ico in st_state.custom_motifs],
                menu_icon="cast",
                orientation="horizontal",
                styles=second_menu_style,
            )
        with col2:
            new_name = st.text_input(
                ":green[Add custom motif with name:]",
                key=f"new_custom_motif{st_state.custom_key}",
            )
            if new_name:
                icon = funny_bootstrap_icons[0]
                # rotate the icons
                funny_bootstrap_icons.append(funny_bootstrap_icons.pop(0))
                st_state.custom_motifs.append(
                    (new_name, icon, pf.Motif(lock_coords=False))
                )
                st_state.custom_key += 1
                st.rerun()

        current_custom_motif = [
            m for m in st_state.custom_motifs if m[0] == motif_selected
        ][0][2]
        custom(current_custom_motif)

    elif selected_motif == "Edit":
        motif_add = False
        edit(st_state.motif_index, st_state.line_index)

    elif selected_motif == "Code":
        code()
        motif_add = False
    # elif selected_motif == 'Undo/Redo':
    #     _, col1, _, col2, _ = st.columns([1] * 5)
    #     with col1:
    #         undo()
    #     with col2:
    #         redo()
    #     motif_add = False

    if motif_add:
        add_motif(origami)


def select_line(f_col1=None, f_subcol2=None, f_subcol3=None):
    # add divider
    st.markdown(
        "<hr style='margin-top:-0em;margin-bottom:-1em' />", unsafe_allow_html=True
    )

    warnings.filterwarnings("error")  # raise warnings as errors
    origami = st_state.origami
    origami_len = len(origami)
    line_index = st_state.line_index
    motif_index = st_state.motif_index
    clicked_indexes = st_state.get(f"origami_click{st_state.ori_click_count}")

    if line_index < 0:
        line_index = 0
    elif motif_index < 0:
        motif_index = 0

    ### Check if the origami is empty and add an empty line
    if origami_len == 0:
        origami.append([])  # add empty line if the Origami is empty
        st_state.code.append("origami.append([]) # Add empty line")
        origami_len = 1

    col1, col2, col3, col4, col5 = st.columns(
        [1, 1.4, 1.4, 1, 0.8],
        vertical_alignment="top",
    )

    with col1:
        if f_col1:
            f_col1()

    ### Select the line in which to add the motif
    with col2:
        line_index = min(line_index, origami_len - 1)
        line_index = st.number_input(
            "**Select line index:**",
            min_value=-1,
            max_value=origami_len - 1,
            value=line_index,
        )
        st_state.current_line_occupied = False

        if line_index < origami_len and origami[line_index]:
            st_state.current_line_occupied = True
        if line_index != st_state.line_index:
            st_state.line_index = line_index
            if clicked_indexes:
                st_state.ori_click_count += 1
            st.rerun()

        subcol1, subcol2 = st.columns(2)

        with subcol1:
            if f_subcol2:
                f_subcol2(line_index + 1)

        ### delete line
        with subcol2:
            # check that the origami has more than one line, or the line is not empty
            if len(st_state.origami) > 1 or st_state.origami[0]:
                # If we want to add a new line, we stop the function here
                if st_state.line_index >= len(st_state.origami):
                    return
                if len(st_state.origami) > 0:
                    del_line = st.button(
                        "Delete line", key="del_line", help="Delete the choosen line"
                    )
                    if del_line:
                        st_state.origami.pop(
                            st_state.line_index
                        )  # remove choosen helix
                        st_state.code.append(f"origami.pop({st_state.line_index})")
                        st_state.line_index -= 1
                        st.rerun()  # rerun app
            # if there is the button to add a motif, add a empty button to align
            elif f_subcol3:
                st.button("", type="tertiary")

    if st_state.line_index >= len(origami):
        return

    ### Motif index selection
    with col3:
        max_val = len(origami[line_index])
        if motif_index > max_val:
            motif_index = max_val
        motif_index = st.number_input(
            "**Select motif index:**",
            min_value=0,
            max_value=max_val,
            value=motif_index,
        )

        if motif_index != st_state.motif_index:
            st_state.motif_index = motif_index
            if clicked_indexes:
                st_state.ori_click_count += 1
            st.rerun()
        subcol1, subcol2 = st.columns(2)

        with subcol1:
            if f_subcol3:
                f_subcol3()
        with subcol2:
            ### Delete motif
            if len(origami[st_state.line_index]) > 0:
                delete_button = st.button(":red[Delete motif]", key="delete_motif")
                if delete_button:
                    if st_state.motif_index == len(origami[st_state.line_index]):
                        st_state.motif_index -= 1
                    origami.pop((st_state.line_index, st_state.motif_index))
                    st_state.code.append(
                        f"origami.pop(({st_state.line_index}, "
                        f"{st_state.motif_index})) # Delete motif"
                    )
                    # st_state.redo = [] # clear the redo buffer
                    st.rerun()  # rerun app

    ### copy/paste motif
    with col4:
        if f_col1:
            copy_motif("preview")
        slice = (line_index, motif_index)
        if 0 <= line_index < len(origami) and 0 <= motif_index < len(
            origami[line_index]
        ):

            motif = origami[slice]
            copy_motif("selected", motif=motif, motif_slice=slice)

        # Paste motif
        if st_state.copied_motif:
            paste_button = st.button("Paste motif", key="paste_motif")
            if paste_button:
                origami.insert((line_index, motif_index), st_state.copied_motif)
                st_state.code.append(
                    st_state.copied_motif_text + f"\norigami.insert(({line_index}, "
                    f"{motif_index}), motif) # Paste motif"
                )
                # st_state.redo = []
                st.rerun()

    # Duplicate line and undo
    with col5:
        undo(key="motif_undo")

        if (
            0 <= line_index < len(origami)
            and origami[line_index]
            and st.button("Duplicate line", key="duplicate_line")
        ):

            # warnings.filterwarnings("ignore") # ignore kl energy warning
            origami.duplicate_line(line_index, insert_idx=len(origami))
            st_state.code.append(
                f"origami.duplicate_line({line_index}, "
                f"insert_idx={len(origami)}) # Duplicate line"
            )
            st_state.line_index = len(origami) - 1
            st.rerun()  # rerun app


def add_motif(origami):
    ### Load the motif and prepare the motif buffer for custom motifs
    motif = st_state.motif
    if not motif:
        return

    # If the motif is a generic motif, generate the motif text buffer
    elif st_state.selected_motif == "Custom":
        ### generate the motif code
        st_state.motif_buffer = "strands = []\n"
        for s in st_state["custom_strands"]:
            coord_txt = ""
            if s._coords:
                coord_txt = (
                    f", coords=pf.Coords({s.coords.array.tolist()}, "
                    f"dummy_ends=({s.coords.dummy_ends[0].tolist()}, "
                    f"{s.coords.dummy_ends[1].tolist()}))"
                )
            st_state.motif_buffer += (
                f"strands.append(pf.Strand('{s}', "
                f"directionality='{s.directionality}', "
                f"start={tuple(s.start)}, "
                f"direction={tuple(s.direction)}{coord_txt}"
                f"))\n"
            )

        st_state.motif_buffer += (
            f"motif = pf.Motif(strands, "
            f"structure = '{motif.structure}', "
            f"lock_coords={motif.lock_coords})"
        )

    ### Show the preview of the motif
    def f_col1():
        # write Preview in small and in orange
        st.write(":orange[Preview:]")
        # display the Motif to be added next
        scrollable_text(motif_text_format(motif))

    ### add line
    def f_subcol2(linex_index):
        # check that the origami has more than one line, or the line is not empty
        if len(st_state.origami) > 1 or st_state.origami[0]:
            add_line = st.button(
                "Add line", key="Insert_line", help=f"Add a new line n° {linex_index}."
            )
            if add_line:
                # add empty line at the choosen position
                origami.insert(linex_index, [])
                st_state.code.append(f"origami.insert({linex_index}, [])")
                st_state.line_index += 1
                st.rerun()

    ### Add motif
    def f_subcol3():
        button_label = "Insert motif"
        if st_state.motif_index == len(origami[st_state.line_index]):
            button_label = "Add motif"
        add_button = st.button(
            f"**:green[{button_label}]**",
            key="insert_motif",
            help="Insert the motif at the choosen position",
        )

        if add_button:
            origami.insert((st_state.line_index, st_state.motif_index), motif)
            st_state.code.append(
                st_state.motif_buffer + f"\norigami.insert(({st_state.line_index}, "
                f"{st_state.motif_index}), motif) # Add motif"
            )
            st_state.motif_buffer = ""  # clear the buffer
            # st_state.redo = [] # clear the redo buffer
            st_state.motif_index += 1
            st.rerun()

    select_line(f_col1, f_subcol2, f_subcol3)


def copy_motif(key="", motif=None, motif_slice=None):
    copy_button = st.button(f"Copy {key} motif", key=f"copy_motif{key}")
    if copy_button:
        if not motif:
            st_state.copied_motif = st_state.motif
            st_state.copied_motif_text = st_state.motif_buffer
        else:
            st_state.copied_motif = motif
            st_state.copied_motif_text = (
                f"motif = origami[{motif_slice}].copy() "
                f"# Copy motif at line {motif_slice[0]}, "
                f"index {motif_slice[1]}"
            )


###
# TAB FUNCTIONS
###


@st.fragment
def generate_custom_motif_text(strand, x_size=50, y_size=10):
    ### update the content
    content = (
        f"<div style='font-family: monospace; "
        f"font-size: {st_state['origami_font_size'] + 8}px;'>"
    )
    current_motif = st_state.motif
    m_pos_to_sym = current_motif.get_position_map()
    s_pos_to_sym = strand.get_position_map()

    for y in range(y_size):
        for x in range(x_size):
            if (x, y) in s_pos_to_sym:
                content += (
                    f'<a href="javascript:void(0);" id="{x},{y}" '
                    f'style="color: #D00000;">{s_pos_to_sym[(x, y)]}</a>'
                )

            elif (x, y) in m_pos_to_sym:  # strand not currently selected
                symbol = m_pos_to_sym[(x, y)]
                content += (
                    f'<a href="javascript:void(0);" id="{x},{y}" '
                    f'style="color: #00856A; opacity: 0.5;">{symbol}</a>'
                )

            else:
                content += (
                    f'<a href="javascript:void(0);" id="{x},{y}" '
                    f'style="color: grey; opacity: 0.5;">•</a>'
                )
        content += "<br />"
    content += "</div>"
    return content


def custom_text_input(current_custom_motif):
    # make a motif list adding a space around the motif
    motif_lines = str(current_custom_motif).split("\n")
    add_line = False
    add_char = False
    if current_custom_motif:
        if any(line[0] != " " for line in motif_lines):
            add_char = True

        motif_list = [
            [" "] * add_char + [char for char in line] + [" "] for line in motif_lines
        ]

        if motif_lines[0].strip():
            add_line = True
            motif_list = [[" "] * (current_custom_motif.num_char + 2)] + motif_list
        if motif_lines[-1].strip():
            motif_list += [[" "] * (current_custom_motif.num_char + 2)]

    else:
        motif_list = [[""]]

    # Add the  5 before the start of each motif.
    for s in current_custom_motif:

        ys, xs = s.prec_pos[1] + int(add_line), s.prec_pos[0] + int(add_char)
        ye, xe = s.next_pos[1] + int(add_line), s.next_pos[0] + int(add_char)
        # add 5' to the start of the strand
        if s.directionality == "53" and motif_list[ys][xs] == " ":
            motif_list[ys][xs] = "5"

        # add 5' to the end of the strand
        elif s.directionality == "35" and motif_list[ye][xe] == " ":
            motif_list[ye][xe] = "5"

    current_custom_motif_str = "\n".join("".join(line) for line in motif_list)

    strand_text = st.text_area(
        "Motif text: draw a motif where each Strand starts " "with a 5",
        value=current_custom_motif_str,
        key="Motif_text",
        help="Draw a motif, where each strand has to start "
        'with "5". If you want to start a strand with 5, '
        "add an additional 5 at the beginning of the "
        "strand.",
    )

    if strand_text != current_custom_motif_str:
        if "5" not in strand_text:
            st.warning('Don\'t forget to start the strand with "5"')
        else:
            if st.button("Convert", key="start_convert", type="primary"):
                new_motif = pf.Motif.from_text(strand_text).strip()
                current_custom_motif.replace_all_strands(new_motif._strands, copy=False)
                current_custom_motif.basepair = new_motif.basepair
                st.rerun()


def structure_converter(current_custom_motif):
    current_struct = current_custom_motif.structure
    current_seq = str(current_custom_motif.sequence)
    if len(current_custom_motif.strands) == current_struct.count("&") + 2:
        current_struct += "&"

    ### Add the structure converter
    col1, col2, col3 = st.columns([4, 4, 1], vertical_alignment="bottom")
    with col1:
        structure = st.text_input(
            "Dot-bracket structure:",
            value=current_struct,
            help="Add a dot-bracket structure to convert it into " "a 3D motif.",
        )
    with col2:
        def_seq = st_state.get("Sequence_converter", current_seq)
        if len(structure) > len(def_seq):
            def_seq += "".join(
                "N" if sym != "&" else "&"
                for sym in structure[len(def_seq) : len(structure) + 1]
            )
        sequence = st.text_input(
            "Sequence:",
            value=def_seq,
            help="Add the sequence of the motif.",
        )
    # check that there is an input
    structure = structure.replace(" ", "")
    sequence = sequence.replace(" ", "")
    optimize = True
    if not structure and not sequence:
        optimize = False
    elif structure != current_struct or sequence != current_seq:
        optimize = True

    if current_custom_motif:
        with col3:
            convert = st.button("Convert", key="start_convert", type="primary")
            optimize = optimize and convert

    if not optimize:
        return

    try:
        new_motif = pf.Motif.from_structure(structure, sequence=sequence)
    except Exception as e:
        st.error(f"Error: {e}")
        return
    current_custom_motif.replace_all_strands(new_motif._strands, copy=False)
    current_custom_motif.basepair = new_motif.basepair
    st.rerun()


def upload_3d_interface(strand, strand_num, current_custom_motif):
    file_3d = st.file_uploader(
        f"3D coordinates (OxDNA format) of " f"**strand {strand_num}**",
        type=["dat", "pdb"],
        key=f"custom_{st_state.upload_key}",
        help='Upload an Oxview configuration ".dat" '
        "file with the 3D coordinates of one strand.",
    )
    dummy_cols = st.columns(2)

    dummy_help = (
        "The dummy base is a base in the coordinates that is not "
        "part of the sequence, but it is used to connect other strands "
        "to the beginning or end of the strand."
    )

    with dummy_cols[0]:
        dummy_start = st.toggle(
            "Start dummy nucleotide", key="dummy_start_custom", help=dummy_help
        )
    with dummy_cols[1]:
        dummy_end = st.toggle(
            "End dummy nucleotide", key="dummy_end_custom", help=dummy_help
        )

    if file_3d:
        pdb_format = file_3d.name.endswith(".pdb")
        strand.coords = pf.Coords.load_from_text(
            file_3d.getvalue().decode("utf-8"),
            dummy_ends=(dummy_start, dummy_end),
            pdb_format=pdb_format,
        )
        update_strand(strand, strand_num, current_custom_motif, coords=True)


def update_strand(strand, strand_num, current_custom_motif, coords=None):
    st_state["custom_strands"][strand_num] = strand
    current_custom_motif.replace_all_strands(st_state["custom_strands"], copy=False)
    if coords:
        st.success(
            "3D coordinates uploaded successfully!", icon=":material/view_in_ar:"
        )
        sleep(2)
        update_file_uploader()
    st.rerun()


def strand_number_button(current_custom_motif, delete=True):
    strands = [s.copy() for s in st_state["custom_strands"]]
    strand_num = st.radio(
        "Select strand:", list(range(len(strands))), index=len(strands) - 1
    )
    strand = strands[strand_num]

    if (
        delete
        and st.button("Delete strand", key="delete_strand")
        and current_custom_motif
    ):

        st_state["custom_strands"].pop(strand_num)
        current_custom_motif.pop(strand_num)
        if not st_state["custom_strands"]:
            st_state["custom_strands"].append(pf.Strand(""))
        st.rerun()

    return strand, strand_num


@st.fragment
def custom(current_custom_motif):
    ### Silence the warnings
    warnings.filterwarnings("ignore")
    st_state.motif = current_custom_motif

    if not st_state["custom_edit"]:
        if st.button(":orange[Edit the motif]", key="edit_strand"):
            st_state["custom_edit"] = True
            st.rerun()
        return

    GeneralEditCommand().execute(current_custom_motif)

    st_state["custom_strands"] = [s.copy() for s in current_custom_motif]
    if not current_custom_motif:
        st_state["custom_strands"] = [pf.Strand("")]

    col1, col2 = st.columns([4, 1.5], vertical_alignment="bottom")
    with col1:
        method = st.segmented_control(
            "Build the motif with:",
            [
                "Structure converter",
                "Full text input",
                "Drawing tool",
            ],
            default="Structure converter",
            help="Structure converter: Convert a dot-bracket"
            " structure into a 3D motif automatically.\n\n"
            "Drawing tool: Draw each strand on a grid by "
            "clicking on the positions where you want the "
            "strand to go; curves and crossings are "
            "calculated automatically.\n\n"
            "Full text input: Create a motif on a blank text"
            " area by typing the motif strand "
            "and basepairing.\n",
        )

    with col2:
        ### add instructions and common symbols
        with st.popover("Common symbols to copy"):
            common_pyroad_sym = ["─", "│", "╭", "╮", "╰", "╯", "^", "┼", "┊", "&"]
            cols = st.columns(len(common_pyroad_sym))
            for i, sym in enumerate(common_pyroad_sym):
                with cols[i]:
                    copy_to_clipboard(sym, sym)

    if method == "Structure converter":
        structure_converter(current_custom_motif)
    elif method == "Full text input":
        custom_text_input(current_custom_motif)
    if method != "Drawing tool":

        with st.popover("Upload 3D coordinates"):
            subcol1, subcol2 = st.columns([1, 4])
            with subcol1:
                strand, strand_nr = strand_number_button(
                    current_custom_motif, delete=False
                )
            with subcol2:
                upload_3d_interface(strand, strand_nr, current_custom_motif)

        if current_custom_motif:
            st.write("Current motif preview:")
            scrollable_text(motif_text_format(current_custom_motif))
            finish_editing(current_custom_motif)
        return

    cols = st.columns(4, vertical_alignment="bottom")
    with cols[0]:
        if st.button("Add strand", key="add_strand"):
            st_state["custom_strands"].append(pf.Strand(""))
            current_custom_motif.append(pf.Strand(""), copy=False, join=False)
    with cols[1]:
        x_dots = st.number_input("Canvas x size", min_value=1, value=64, key="x_size")
    with cols[2]:
        y_dots = st.number_input("Canvas y size", min_value=1, value=8, key="y_size")
    with cols[3]:
        if st.button("Clear", key="clear_strand"):
            st_state["custom_strands"] = [pf.Strand("")]
            current_custom_motif.replace_all_strands([pf.Strand("")], copy=False)

    # Update and display the motif text
    col1, col2 = st.columns([2, 7])
    with col1:
        strand, strand_num = strand_number_button(current_custom_motif)
    with col2:
        try:
            clicked = click_detector(
                generate_custom_motif_text(strand, x_size=x_dots, y_size=y_dots),
                key="custom_motif_click_detector",
            )
        except Exception as e:
            st.error(str(e))

    def move_strand(clicked, add_nucl=None):
        # nonlocal strand
        pos = tuple([int(i) for i in clicked.split(",")])

        ### if the strand is empty, add the first position
        if not str(strand):
            strand._start = pos
            strand.strand = "─"

        ### if the strand has only one position, add the strand direction
        elif len(str(strand)) == 1:
            start = strand.start
            direction = (pos[0] - start[0], pos[1] - start[1])
            if direction[0] and direction[1]:
                st.error("Invalid strand, the strand can't form diagonal lines.")
                return
            # the Start horizontally
            if direction[0]:
                sym = "-"
                distance = abs(direction[0]) + 1
                if direction[0] > 0:
                    strand.direction = (1, 0)
                else:
                    strand.direction = (-1, 0)
            # the Start vertically
            else:
                sym = "|"
                distance = abs(direction[1]) + 1
                if direction[1] > 0:
                    strand.direction = (0, 1)
                else:
                    strand.direction = (0, -1)
            if add_nucl:
                sym = add_nucl
            strand.strand = sym * distance
        else:
            ### if the clicked position is the start, remove it and update the start
            if pos == strand.start:
                strand._start = strand.positions[0]
                strand._direction = strand.directions[0]
                strand.strand = strand.strand[1:]

            ### check the direction and in case add the new position
            else:
                # calculate the direction
                direction = (pos[0] - strand.end[0], pos[1] - strand.end[1])
                # calculate the direction symbol
                if add_nucl:
                    sym = add_nucl
                elif direction[0]:
                    sym = "─"
                else:
                    sym = "│"
                # calulcate the difference between the new direction and
                # the end direction (use to determine the symbol to add)
                normal_direction = tuple(
                    [i // abs(i) if i != 0 else 0 for i in direction]
                )
                direction_difference = (
                    normal_direction[0] - strand.end_direction[0],
                    normal_direction[1] - strand.end_direction[1],
                )
                distance = max(abs(pos[0] - strand.end[0]), abs(pos[1] - strand.end[1]))
                # if the direction isn't changed, add the last symbol
                # to match the strand
                if normal_direction == strand.end_direction:
                    strand.strand = strand.strand + sym * distance
                elif normal_direction not in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
                    st.error("Invalid strand. The strand can't form diagonal lines.")
                elif direction_difference == (-1, 1):
                    strand.strand = strand.strand[:-1] + "╮" + sym * distance
                elif direction_difference == (1, 1):
                    strand.strand = strand.strand[:-1] + "╭" + sym * distance
                elif direction_difference == (-1, -1):
                    strand.strand = strand.strand[:-1] + "╯" + sym * distance
                elif direction_difference == (1, -1):
                    strand.strand = strand.strand[:-1] + "╰" + sym * distance
                else:
                    st.error("Invalid strand")

    ### Update the strand on click
    if clicked and clicked != st_state["last_clicked_position"]:
        try:
            move_strand(clicked)
            st_state["last_clicked_position"] = clicked
        except pf.MotifStructureError as e:
            st.error(str(e))

    ### Show the strand options
    col1, col2, col3, col4, col5 = st.columns([1, 1, 1, 1, 5])
    with col1:
        start_x = st.number_input(
            "Start x:", min_value=0, value=strand.start[0], key="start_x_custom"
        )
    with col2:
        start_y = st.number_input(
            "Start y:", min_value=0, value=strand.start[1], key="start_y_custom"
        )
    with col3:
        strand_direction_ind = [d for d in pf.Direction].index(strand.direction)
        new_dir = st.selectbox(
            "Start direction:",
            direction_list,
            index=strand_direction_ind,
            key="dir_custom",
        )
        new_dir_tuple = pf.Direction[new_dir]
    with col4:
        seq_dir = st.selectbox(
            "Directionality:",
            ["35", "53"],
            index=["35", "53"].index(strand.directionality),
            key="seq_dir_custom",
        )
    with col5:
        new_strand = st.text_input(
            f"New strand (strand directionality: " f"{strand.directionality}) ",
            value=str(strand),
            key="strand_custom",
        )

    ### update the strand
    strand.start = (start_x, start_y)
    strand.direction = new_dir_tuple
    strand.strand = new_strand
    strand.directionality = seq_dir
    with st.popover("Add strand 3D coordinates"):
        upload_3d_interface(strand, strand_num, current_custom_motif)

    # compare the strand with the version strand: if the strand is different,
    #  update it and rerun
    old_strand = st_state["custom_strands"][strand_num]
    if strand != old_strand:
        update_strand(strand, strand_num, current_custom_motif)

    ### check the base pair symbols of the motif
    current_structure = current_custom_motif.structure
    new_db = st.text_input(
        "Add dot-bracket notation:",
        value=current_structure,
        key="structure_custom",
        help="Add the dot-bracket notation of the motif for each "
        'strand, separated by a "&". If the paired bases are '
        'more than one position apart, the pairing symbol "┊"'
        " is not shown.",
    )
    if new_db != current_structure:
        current_custom_motif.structure = new_db

    finish_editing(current_custom_motif)
    st.stop()


def finish_editing(current_custom_motif):
    if st.button(":green[Finish editing to add motif to the origami]", key="end_edit"):
        st_state["custom_edit"] = False
        current_custom_motif.replace_all_strands(st_state["custom_strands"], copy=False)
        st.rerun()


def update_code(code_text, return_origami=False):
    # Check for forbidden keywords
    forbidden_keywords = ["exec", "eval", "compile", " open", " os", " sys", " __"]
    for kw in forbidden_keywords:
        if kw in code_text:
            st.error(
                f"Forbidden keyword detected in the code: '{kw}'",
                icon=":material/sentiment_extremely_dissatisfied:",
            )
            return False

    # Define the local environment in which the user's code will be executed
    local_context = {"origami": st_state.origami}
    # Attempt to execute the code safely
    try:
        exec(code_text, {"__builtins__": __builtins__, "pf": pf}, local_context)
        # Retrieve the modified origami variable
        origami = [
            v
            for k, v in local_context.items()
            if isinstance(v, pf.Origami) and k == "origami"
        ][0]
        render_target = local_context.get("RENDER_TARGET", None)
        if render_target and isinstance(render_target, pf.Origami):
            origami = render_target

        if origami:
            # select the end of the origami
            st_state.line_index = len(origami) - 1
            st_state.motif_index = len(origami[-1])

    except Exception as e:
        st.error(f"Error in executing the code: {e}")
        return False

    if return_origami:
        return origami

    st.success("Nanostructure updated successfully!")

    st_state.origami = origami
    st_state.code = code_text.split("\n\n")
    st.rerun()


def code():
    if "last_code_id" not in st_state:
        st_state["last_code_id"] = ""
    code_text = "\n\n".join(st_state.code)

    col1, col2, col3 = st.columns([1, 1, 1], gap="large", vertical_alignment="center")
    with col1:
        st.link_button(
            "Check the documentation!",
            "https://pyfurnace.readthedocs.io/en/latest/api.html",
            icon=":material/document_search:",
            help="The code editor run python code, so you can program your "
            "RNA origami. The GUI will render the variable named 'origami'"
            " as the RNA origami structure. If you want to render a different"
            " structure, assign it to the variable RENDER_TARGET. If you need"
            " more info, click the button to open the documentation in a new tab.",
        )
    with col2:
        render_lines = st.slider(
            "Number of lines to render:",
            min_value=1,
            max_value=200,
            value=10,
            key="render_lines",
        )
    with col3:
        wrap = st.toggle(
            "Wrap lines",
            value=False,
            key="wrap_lines",
            help="Wrap lines in the code editor to fit the screen.",
        )

    response_dict = code_editor(
        code_text,
        buttons=code_editor_buttons,
        height=render_lines,
        options={"showLineNumbers": True, "wrap": wrap},
        allow_reset=True,
        key="pyroad_code_editor",
    )

    if response_dict["id"] != st_state["last_code_id"] and (
        response_dict["type"] == "submit" or response_dict["type"] == "selection"
    ):
        st_state["last_code_id"] = response_dict["id"]
        update_code(response_dict["text"])

    else:
        st.success("Structure updated!")


def undo(key=""):
    if not st_state.origami:
        return
    if len(st_state.code) <= 1:
        st.warning("Nothing to undo.")
        return
    undo_button = st.button(
        ":red[Undo]",
        key=f"undo{key}",
        use_container_width=False,
        help="Undo the last action.",
    )
    if not undo_button:
        return
    # st_state.redo.append(st_state.code[-1])
    update_code("\n\n".join(st_state.code[:-1]))


# def redo():
#     if not st_state.redo:
#         # st.warning("Nothing to redo.")
#         return
#     redo_button = st.button("Redo", help="Redo the last action.")
#     if not redo_button:
#         return
#     last_action = st_state.redo.pop()
#     update_code('\n\n'.join(st_state.code + [last_action]))
#     st.rerun()


def origami_select_display():
    option_data = {
        "Origami 2D view": "bi bi-square",
        "Origami 3D view": "bi bi-box",
        "Origami split view": "bi bi-window-split",
        "Folding barrier": "bi bi-graph-up",
    }

    selected_display = option_menu(
        None,
        list(option_data.keys()),
        icons=list(option_data.values()),
        orientation="horizontal",
        key="DisplayMenu",
        styles=main_menu_style,
    )
    return selected_display


@st.fragment
def origami_build_view(selected_display):
    ### Display the RNA origami structure with clickable elements
    #  and modify them in case
    warnings.filterwarnings("ignore")  # ignore numpy warnings

    if selected_display == "Origami 2D view":
        clicked = display_origami()
        clicked_options(clicked)

    elif selected_display == "Origami 3D view":
        display3d()

    elif selected_display == "Origami split view":
        col1, col2 = st.columns(2)
        with col1:
            clicked = display_origami()
        with col2:
            display3d()
        clicked_options(clicked)

    elif selected_display == "Folding barrier":
        barriers, score = st_state.origami.folding_barriers()
        clicked = display_origami(barriers=barriers)
        clicked_options(clicked)

        help_message = (
            "The folding barriers are a simplified calculation to check "
            "if a structure is feasible to form co-transcriptionally. "
            "They are calculated considering that pseudoknots lock the "
            "structure in a specific conformation; preventing the "
            "formation of successive stems. A barrier happens when "
            "a stem is open before the formation of a pseudoknots, "
            "and is closed after the formation of the pseudoknot "
            "(after a 150-bases delay). The penalty is calculated as W"
            " (weak barrier): 1 point, S (strong barrier): 2 points, "
            "and the total penalty is the sum of the points."
        )

        cols = st.columns(4, vertical_alignment="center")
        with cols[0]:
            if score < 30:
                st.success(f"Low folding barrier penalty: {score}.")
            elif score < 100:
                st.warning(f"Medium folding barrier penalty: {score}.")
            else:
                st.error(f"High folding barrier penalty: {score}.")
        with cols[1]:
            st.markdown("", help=help_message)
        with cols[2]:
            if score > 30:
                if st.button(
                    "Try to fix the folding pathway",
                    key="fix_pathway",
                    help="Try to change the starting position of the sequence, "
                    "keeping the same structure, to reduce the folding "
                    "pathway penalty. This function is designed for "
                    "canonical RNA Origami structures, it may not work "
                    "for other complex structures.",
                ):

                    # calculate the simplfied dot bracket
                    with st.spinner("Optimizing the folding pathway..."):
                        code_lines = st_state.code
                        code_lines += [
                            "\norigami = origami.improve_"
                            "folding_pathway"
                            "(kl_delay=150)\n"
                        ]
                        update_code("\n".join(code_lines))


def display3d():
    origami = st_state.origami

    if st_state.gradient:
        seq = origami.sequence.replace("&", "")
        index_colors = list(range(len(seq)))

    else:
        m_slice = (st_state.line_index, st_state.motif_index)
        index_colors = ()
        if m_slice[0] < len(origami) and m_slice[1] < len(origami[m_slice[0]]):
            motif = origami[m_slice[0]][m_slice[1]]
            motif_shift = origami.index_shift_map[m_slice]
            index_colors = [0] * len(origami.sequence.replace("&", ""))
            for pos in motif.seq_positions:
                shifted = tuple([pos[0] + motif_shift[0], pos[1] + motif_shift[1]])
                index_colors[origami.seq_positions.index(shifted)] = 1

    if index_colors:
        for s in origami.strands:
            for protein in s.coords.proteins:
                index_colors += [0] * len(protein)

    conf, topo = origami.save_3d_model("origami", return_text=True)
    oxview_from_text(
        configuration=conf,  # path to the configuration file
        topology=topo,  # path to the topology file
        width="99%",  # width of the viewer frame
        colormap=st_state.oxview_colormap,  # colormap for the viewer
        height=st_state.oxview_frame_size,  # height of the viewer frame
        index_colors=index_colors,  # color the bases in the viewer
        frame_id=1,
        key="display_nano",
    )


def build_origami_content(barriers=None):
    barriers_colors = {"▂": "#FFBA08", "▄": "#FFBA08", "█": "#D00000"}
    highlight_color = "#D00000"
    normal_color = "#333333"  # 80% black

    origami = st_state.origami
    motif = origami.assembled

    if barriers:
        origami_lines = origami.barrier_repr(barriers=barriers, return_list=True)
    else:
        if st_state.get("to_road"):
            origami_str = origami.to_road()

        else:
            origami_str = str(origami)
        origami_lines = origami_str.split("\n")

    # create color gradient
    if st_state.gradient:
        # create a dictionary from positions to index
        pos_to_index = {pos: ind for ind, pos in enumerate(motif.seq_positions)}
        tot_len = 0
        for s in origami.strands:
            tot_len += len(s.sequence)
            for protein in s.coords.proteins:
                tot_len += len(protein)
        cmap = plt.get_cmap(st_state.oxview_colormap)
        oxview_offset = tot_len // 5
        c_map = [
            mcolors.to_hex(cmap(i)) for i in np.linspace(0, 1, tot_len + oxview_offset)
        ]
        c_map = c_map[oxview_offset:]

    # Prepare the string to add 5' and 3' symbols for the strands
    motif_list = (
        [[" "] * (motif.num_char + 2)]
        + [[" "] + [char for char in line] + [" "] for line in origami_lines]
        + [[" "] * (motif.num_char + 2)]
    )

    for s in motif:  # Add the 5' and 3' symbols to the motif as 1 and 2
        if not s.sequence:
            continue
        ys, xs = s.prec_pos[1] + 1, s.prec_pos[0] + 1
        ye, xe = s.next_pos[1] + 1, s.next_pos[0] + 1
        # add 5' to the start of the strand
        if s.sequence and s[0] not in "35" and motif_list[ys][xs] == " ":
            motif_list[ys][xs] = "1" if s.directionality == "53" else "2"
        # add 5' to the end of the strand
        if s.sequence and s[-1] not in "35" and motif_list[ye][xe] == " ":
            motif_list[ye][xe] = "2" if s.directionality == "53" else "1"
    origami_list = ["".join(line) for line in motif_list]

    content = (
        "<div style='font-family: monospace; white-space: nowrap; "
        "overflow-x: auto;'>"
    )

    # add a column with the line number
    content += (
        f'<div style="display:inline-block; font-family: monospace;'
        f' font-size: {st_state.origami_font_size}px;">Line:<br />'
    )

    line_nr = -2
    origami_list_len = len(origami_list)
    span_text = '<span style="font-family: monospace; '

    for y, line in enumerate(origami_list):
        line_color = normal_color
        motif_slice = None
        current_line_nr = [
            origami.pos_index_map[(x - 1, y - 1)]
            for x, _ in enumerate(line)
            if (x - 1, y - 1) in origami.pos_index_map
        ]

        next_line = origami_list[y + 1] if (y + 1) < origami_list_len else []

        # check the next line index
        next_line_nr = [
            origami.pos_index_map[(x - 1, y)]
            for x, _ in enumerate(next_line)
            if (x - 1, y) in origami.pos_index_map
        ]

        if current_line_nr:
            current_line_nr = current_line_nr[0][0]  # get the first line number
        if next_line_nr:
            next_line_nr = next_line_nr[0][0]

        # is a new origami line
        if line_nr != current_line_nr and isinstance(current_line_nr, int):
            if current_line_nr == st_state.line_index:
                line_color = highlight_color
            content += (
                span_text
                + (
                    f'color: {line_color}; line-height:1;">'
                    f"{current_line_nr})</span>"
                )
                + (4 - len(str(current_line_nr))) * "&nbsp;"
            )
            line_nr = current_line_nr

        # the origami line is empty
        elif (
            isinstance(next_line_nr, int)
            and line_nr < st_state.line_index < next_line_nr
        ):
            content += span_text + 'color: #D52919; line-height:1;">_____</span>'

        # is not a new origami line
        else:
            content += (
                span_text + 'line-height:1;">&nbsp;&nbsp;&nbsp;' "&nbsp;&nbsp;</span>"
            )

        hit = False  # to highlight the first symbol
        for x, char in enumerate(line):
            ori_pos = (x - 1, y - 1)
            color = normal_color
            if barriers_colors and char in barriers_colors:
                color = barriers_colors[char]

            if char == " ":
                content += span_text + 'line-height:1;">&nbsp;</span>'

            elif char == "1":
                content += (
                    span_text + f'color: {highlight_color}; line-height:1;">5</span>'
                )
            elif char == "2":
                content += (
                    span_text + f'color: {highlight_color}; line-height:1;">3</span>'
                )

            # do not highlight the base pair in red
            elif char in pf.bp_symbols:
                content += span_text + f'color: {color}; line-height:1;">{char}</span>'

            elif ori_pos in origami.pos_index_map:  # a motif symbol
                motif_slice = origami.pos_index_map[ori_pos]
                if st_state.gradient:
                    index = pos_to_index.get(ori_pos)
                    if index is not None:
                        color = c_map[index]

                # This is the selected motif
                if (
                    motif_slice
                    and motif_slice[0] == st_state.line_index
                    and motif_slice[1] == st_state.motif_index
                ):

                    # Don't add color if the gradient is not active
                    if st_state.gradient:
                        color = normal_color

                    # the first symbol is yellow
                    elif not hit:
                        color = "#FF8800"
                        hit = True

                    else:
                        color = highlight_color

                content += (
                    f'<a style="font-family: monospace; color: {color}; '
                    'line-height:1;" href="javascript:void(0);" '
                    f'id="{motif_slice[0]},{motif_slice[1]},{x - 1},'
                    f'{y - 1}">{char}</a>'
                )

            # is a junction symbol
            else:
                content += span_text + f'color: {color}; line-height:1;">{char}</span>'

        # Check if you wanna add the cursor
        if motif_slice and motif_slice[0] == st_state.line_index:
            if len(origami[motif_slice[0]]) == st_state.motif_index:
                content += (
                    span_text + f'color: {highlight_color}; line-height:1;">│</span>'
                )
        content += "<br />"

    if line_nr < st_state.line_index:
        # the origami line is empty
        content += span_text + 'color: #D52919; line-height:1;">_____</span>'

    content += "</div>"
    content += "</div>"
    return content


def display_origami(barriers=None):
    ### SHORTCUTS NOT READY YET!!!
    # if not split_view:
    #     def move_selection(x = 0, y = 0):
    #         st_state.line_index += y
    #         st_state.motif_index += x
    #         # st.rerun()

    #     with st.expander("Shortcuts to move over the Origami", expanded=False):
    #         st.write("Use the arrow keys to move the selection over the Origami "
    #                   "lines and motifs. Not available in split view.")
    #         col1, col2, col3, col4 = st.columns(4)
    #         with col1:
    #             streamlit_shortcuts.button("Move Up", on_click=move_selection,
    #                                           shortcut="Shift+ArrowUp",
    #                                            hint=True, kwargs={"y" : -1})
    #         with col2:
    #             streamlit_shortcuts.button("Move Down",
    #                                         on_click=move_selection,
    #                                         shortcut="Shift+ArrowDown",
    #                                          hint=True, kwargs={"y" : 1})
    #         with col3:
    #             streamlit_shortcuts.button("Move Left",
    #                                        on_click=move_selection,
    #                                         shortcut="Shift+ArrowLeft",
    #                                         hint=True, kwargs={"x" : -1})
    #         with col4:
    #             streamlit_shortcuts.button("Move Right",
    #                                         on_click=move_selection,
    #                                         shortcut="Shift+ArrowRight",
    #                                         hint=True, kwargs={"x" : 1})

    try:
        content = build_origami_content(barriers=barriers)
        clicked = click_detector(
            content, key=f"origami_click{st_state.ori_click_count}"
        )
    except pf.MotifStructureError as e:
        st.error(f"Structure error: \n {e}", icon=":material/personal_injury:")
        undo(key="error")
        st.stop()
    except pf.AmbiguosStructure as e:
        st.error(f"Ambigouse structure: \n {e}", icon=":material/theater_comedy:")
        undo(key="warning")
        st.write("You can try flipping the motif or changing the sequence direction.")
        st.stop()
    return clicked


def clicked_options(clicked):
    # undo()
    # col1, _ , _, col2, _ = st.columns([1] * 5)
    # with col1:
    #     undo()
    # with col2:
    #     redo()

    if clicked:
        y_slice, x_slice, x_pos, y_pos = [int(i) for i in clicked.split(",")]
        motif_slice = (y_slice, x_slice)
        pos = (x_pos, y_pos)

        # Write the selected motif on screen
        try:
            motif = st_state.origami[motif_slice]
        except IndexError:
            return  # this is a problem with updating the keys

        motif_class_name = motif.__class__.__name__
        if motif_class_name == "Motif":  # is custom Motif
            motif_class_name = "Custom " + motif_class_name

        nucl_text = ""
        # Indicate the nucleotide index:
        if pos in st_state.origami.seq_positions:
            ind = st_state.origami.seq_positions.index(pos)
            nucl_text += f":orange[nucleotide index {ind}] in"

        ### highlight the selected motif
        if (
            st_state.line_index != motif_slice[0]
            or st_state.motif_index != motif_slice[1]
        ):
            st_state.line_index = motif_slice[0]
            st_state.motif_index = motif_slice[1]
            st.rerun()

        st.markdown(
            f"Selected {nucl_text} :orange[{motif_class_name}]: "
            f":orange[line {motif_slice[0]}], "
            f":orange[motif {motif_slice[1]}]"
        )


def display_structure_sequence():
    st.write("#### Origami structure and sequence constraints:")
    origami = st_state.origami
    if origami.structure:
        structure_list = [char for char in origami.structure]
        sequence_list = [char for char in origami.sequence]
        motif_slice = (st_state.line_index, st_state.motif_index)
        indexes = []
        if motif_slice in origami.index_shift_map:
            motif = origami[motif_slice]
            shifts_x, shift_y = origami.index_shift_map[motif_slice]
            for pos in motif.seq_positions:
                indexes.append(
                    origami.seq_positions.index((pos[0] + shifts_x, pos[1] + shift_y))
                )

        for i in range(len(structure_list)):
            if i in indexes:  # highlight the selected motif
                structure_list[i] = (
                    f'<span style="color: #D52919; line-height:1;">'
                    f"{html.escape(structure_list[i])}</span>"
                )
                sequence_list[i] = (
                    f'<span style="color: #D52919; line-height:1;">'
                    f"{html.escape(sequence_list[i])}</span>"
                )
            else:
                structure_list[i] = html.escape(structure_list[i])
                sequence_list[i] = html.escape(sequence_list[i])
            # structure_list[i] = f":orange[{structure_list[i]}]"
            # sequence_list[i] = f":orange[{sequence_list[i]}]"

        ### Pseudoknots info
        pseudoknot_text = (
            "; ".join(
                [f"id: {k}, " + str(val) for k, val in origami.pseudoknots.items()]
            )
            + ";"
        )
        remove_syms = str.maketrans("", "", "{}'\"")
        pseudoknot_text = html.escape(pseudoknot_text.translate(remove_syms))

        st.markdown(
            f"**Structure length: "
            f':green[{str(len(origami.structure.replace("&", "")))}]**'
        )

        scrollable_text(
            f">Origami</br>{''.join(sequence_list)}</br>"
            f"{''.join(structure_list)}</br></br>"
            f"Pseudoknots info:</br>{pseudoknot_text}"
        )

        st.markdown(
            "<hr style='margin-top:-0em;margin-bottom:-1em;"
            "color:white;border-color:white' />",
            unsafe_allow_html=True,
        )

        st_state.generate_structure = origami.structure
        st_state.generate_sequence = str(origami.sequence)
        st_state.generate_pseudoknots = pseudoknot_text
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.page_link(
                "pages/2_Generate.py",
                label="**:orange[Generate the sequence]**",
                icon=":material/network_node:",
            )
        with col2:
            copy_to_clipboard(origami.structure, "Structure")
        with col3:
            copy_to_clipboard(origami.sequence, "Sequence")
        with col4:
            copy_to_clipboard(pseudoknot_text, "Pseudoknots")


def edit(x, y):
    origami = st_state.origami
    ### select the motif to change
    st_state.modified_motif_text = ""  # initialize the modification
    try:
        motif_slice = (y, x)
        if y >= len(origami) or x >= len(origami[y]):
            select_line()
            return
        motif = origami[motif_slice]
    except KeyError:
        return  # this is a problem with updating the keys
    motif_class_name = motif.__class__.__name__
    if motif_class_name == "Motif":  # is custom Motif
        motif_class_name = "Custom " + motif_class_name
    elif isinstance(motif, pf.KissingLoop):
        motif_class_name = "KissingLoops"

    ### Directly modify the motif and make a copy for the advanced feature
    command = f"{motif_class_name}Command"
    if command in globals():
        globals()[command]().execute(motif)
    else:
        GeneralEditCommand().execute(motif)

    if st_state.modified_motif_text:
        st_state.code.append(
            f"motif = origami[({motif_slice[0]}, "
            f"{motif_slice[1]})] "
            f"# motif slice: line {motif_slice[0]}, "
            f"index {motif_slice[1]}" + st_state.modified_motif_text
        )
        st_state.modified_motif_text = ""
        st.rerun()

    select_line()
    st_state["mot_adv_edit"] = motif.copy()
    advaced_edit(motif_slice)


@st.fragment
def advaced_edit(motif_slice):
    ### try to change each strand
    with st.popover(
        "Advanced feature, modify the strands of the selected motif:",
        use_container_width=True,
    ):

        cols = st.columns(
            [1.4, 1.4, 0.7, 0.9, 1.4, 1.4, 0.8], vertical_alignment="bottom"
        )
        with cols[0]:
            flip_vert = st.toggle("Flip vertically", value=False)
        with cols[1]:
            flip_h = st.toggle("Flip horizontally", value=False)
        with cols[2]:
            if flip_vert or flip_h:
                flip = st.button("Flip")
                flip_vert &= flip
                flip_h &= flip

        with cols[3]:
            if st.button("Rotate clockwise"):
                st_state.mot_adv_edit.rotate()
                st_state.modified_motif_text += "\nmotif.rotate()"
                st.rerun(scope="fragment")

        if flip_h or flip_vert:
            st_state.mot_adv_edit.flip(horizontally=flip_h, vertically=flip_vert)
            st_state.modified_motif_text += (
                f"\nmotif.flip(" f"horizontally={flip_h}, " f"vertically={flip_vert})"
            )
            st.rerun(scope="fragment")

        with cols[4]:
            x_shift = st.number_input("Horizontal shift amount:", value=0)
        with cols[5]:
            y_shift = st.number_input("Vertical shift amount:", value=0)

        with cols[6]:
            if x_shift != 0 or y_shift != 0:
                if st.button(
                    "Shift motif",
                    help="Shift the motif by the specified x and y values. "
                    "The shift is applied to all strands in the motif.",
                ):
                    st_state.mot_adv_edit.shift((x_shift, y_shift))
                    st_state.modified_motif_text += (
                        f"\nmotif.shift((" f"{x_shift}, {y_shift}))"
                    )
                    st.rerun(scope="fragment")

        for i, s in enumerate(st.session_state.mot_adv_edit):
            col1, col2, col3, col4, col5 = st.columns([1, 1, 1, 1, 3])
            with col1:
                start_x = st.number_input(
                    f"{i}) Start x:",
                    min_value=0,
                    value=s.start[0],
                )
            with col2:
                start_y = st.number_input(
                    f"{i}) Start y:",
                    min_value=0,
                    value=s.start[1],
                )

            with col3:
                strand_direction_ind = [d for d in pf.Direction].index(s.direction)
                new_dir = st.selectbox(
                    f"{i}) Start direction:",
                    direction_list,
                    index=strand_direction_ind,
                )
                new_dir_tuple = pf.Direction[new_dir]

            with col4:
                seq_dir = st.selectbox(
                    f"{i}) Directionality:",
                    ["35", "53"],
                    index=["35", "53"].index(s.directionality),
                )

            with col5:
                new_strand = st.text_input(
                    f"{i}) New strand (strand directionality: " f"{s.directionality}) ",
                    value=str(s),
                    help=" The strand contains all "
                    "the symbols of the structure "
                    "(nucleotides and direction "
                    "symbols). ",
                )

            ### check if the strand is modified and update the motif
            strand_stripped = new_strand.strip()
            if s.start != (start_x, start_y):
                s.start = (start_x, start_y)
                st_state.modified_motif_text += (
                    f"\nmotif[{i}].start = " f"({start_x}, {start_y})"
                )
                st.rerun(scope="fragment")

            if s.directionality != seq_dir:
                s.directionality = seq_dir
                st_state.modified_motif_text += (
                    f"\nmotif[{i}].directionality " f'= "{seq_dir}"'
                )
                st.rerun(scope="fragment")

            if s.direction != new_dir_tuple:
                s.direction = new_dir_tuple
                st_state.modified_motif_text += (
                    f"\nmotif[{i}].direction = " f"{new_dir_tuple}"
                )
                st.rerun(scope="fragment")

            if s.strand != strand_stripped:
                s.strand = strand_stripped
                st_state.modified_motif_text += (
                    f"\nmotif[{i}].strand = " f'"{strand_stripped}"'
                )
                st.rerun(scope="fragment")

        ### check the base pair symbols of the motif
        current_structure = st_state.mot_adv_edit.structure
        new_db = st.text_input(
            "Add dot-bracket notation:",
            value=current_structure,
            help="Add the dot-bracket notation of the motif for"
            ' each strand, separated by a "&". If the '
            "paired bases are more than one position "
            'apart, the pairing symbol "┊" is not shown.',
        )

        if new_db != current_structure:
            if not new_db:
                st_state.mot_adv_edit.autopairing = True
                st_state.modified_motif_text += "\nmotif.autopairing = True"
            else:
                st_state.mot_adv_edit.structure = new_db
                st_state.modified_motif_text += f'\nmotif.structure = "{new_db}"'
            st.rerun(scope="fragment")

        ### update the motif
        try:
            st.write(":orange[Preview:]")

            scrollable_text(motif_text_format(st_state.mot_adv_edit))
            if not st_state.modified_motif_text:
                st.write(":green[Structure is updated]")
            else:
                st.warning(
                    "Be careful before updating, there is a high risk of "
                    "breaking the structure."
                )
                update = st.button("Update strands")

                if update:
                    st_state.code.append(
                        f"motif = origami[({motif_slice[0]}, "
                        f"{motif_slice[1]})] # select the motif "
                        f"at index, line: {motif_slice}" + st_state.modified_motif_text
                    )

                    st_state.origami[motif_slice] = st_state.mot_adv_edit.copy()
                    st.rerun()

        except pf.MotifStructureError as e:
            st.error(e)


def scrollable_text(text: str):
    st.markdown(
        f"""
        <style>
        .scroll-box {{
            overflow-x: auto;
            white-space: nowrap;
            background-color: #fafafa;
            padding: 10px;
            border-radius: 20px; /* Rounded corners */
            font-size: {st_state.origami_font_size + 1}px;
        }}
        </style>
        <div class="scroll-box">{text}</div>
        """,
        unsafe_allow_html=True,
    )
