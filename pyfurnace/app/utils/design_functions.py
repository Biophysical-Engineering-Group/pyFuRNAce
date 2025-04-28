from time import sleep
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import streamlit as st
from streamlit_option_menu import option_menu
from st_click_detector import click_detector
from st_oxview import oxview_from_text
from code_editor import code_editor
import warnings
# from streamlit_shortcuts import button
import matplotlib.pyplot as plt
### Custom libraries
from . import main_menu_style, second_menu_style, copy_to_clipboard
import pyfurnace as pf
from utils.commands import *

code_editor_buttons = [
                        {
                        "name": "Copy",
                        "feather": "Copy",
                        "hasText": True,
                        "alwaysOn": True,
                        "commands": ["copyAll"],
                        "style": {"top": "0.46rem", "right": "0.4rem"}
                        },
                        {
                        "name": "Run",
                        "feather": "Play",
                        "primary": True,
                        "hasText": True,
                        "showWithIcon": True,
                        "commands": ["submit"],
                        "style": {"bottom": "0.44rem", "right": "0.4rem"}
                        },
                        ]

funny_bootstrap_icons = ['robot', 'trash', 'umbrella', 'camera', 'cart', 'cpu', 'cup-straw', 'trophy', 'palette', 'cup-straw', 'camera-reels', 'puzzle', 'hourglass-split', 'mortarboard']

direction_list = pf.Direction.names()

def origami_general_options(origami, expanded=True):
    with st.expander('General settings', expanded=expanded):
        cols = st.columns(3)

        ### select alignment
        with cols[0]:
            # with st.expander('Align Lines:'):
                alignment = st.radio("Alignment type:", ['To the left', 'Junctions alignment: first', 'Line center'], index=0, key='alignment')
                
                ### add the alignment to the code
                if st.button('Set alignment', key='set_alignment'):
                    alignment = alignment.split()[-1]
                    origami.align = alignment
                    st.session_state.code.append(f"origami.align = '{alignment}' # Align the lines")
        
        # not really useful
        # ### select the sequence direction
        # with cols[1]:
        #     top_strand_dir = st.radio("Sequence direction in the helix top strand:", ['53', '35'], key='top_strand',  
        #                               help=""" PyFuRNAce uses the 5' to 3' direction as the top strand of the helix.
        #                                         To allow compatibility with previus ROAD blueprint, you can change the direction of the top strand.
        #                                         Some functionalialities may not work properly with the 3' to 5' direction.""")
        #     if top_strand_dir == '53':
        #         st.session_state.flip = False
        #     else:
        #         st.session_state.flip = True

        with cols[1]:
            st.toggle('Optimize the blueprint for ROAD', 
                      value=False,
                      key='to_road', 
                      help='Optimize the blueprint for the ROAD software. This substitutes the Kissing Loops base pairings with "*"; and the short stem base pairings with "!".')

        # Not really useful
        # with cols[2]:
        #     def submit_ss_rna():
        #         st.session_state.code.append(f"origami.ss_assembly = {st.session_state.single_stranded}")
        #         st.session_state.origami.ss_assembly = st.session_state.single_stranded

        #     st.toggle("Single stranded 3D assembly", key='single_stranded', on_change=submit_ss_rna, 
        #               help='The Origami 3d assembly is created concatenating each RNA strand rather than concatenating the motifs')

        with cols[2]:
            new_sticky = st.toggle('Stick the motif menu at the top', 
                                   value=st.session_state.motif_menu_sticky, 
                                   key='sticky_menu', 
                                   help='Keep the motif menu and origami visualization menu to stick to the top of the page.')
            if new_sticky != st.session_state.motif_menu_sticky:
                st.session_state.motif_menu_sticky = new_sticky
                st.rerun()
            
        # font_size = st.slider('Font size:', min_value=4, max_value=40, value=14, key='general_font_size', help='Change the font size of the motif preview.')
        # st.markdown(
        #     f"""
        #     <style>
        #     /* Apply global font size */
        #     html, body, [class*="css"]  {{
        #         font-size: {font_size}px;
        #     }}
        #     </style>
        #     """,
        #     unsafe_allow_html=True,
        # )

        col1, col2 = st.columns(2)
        with col1:
            font_size = st.slider('Origami font size', 
                                  min_value=2, 
                                  max_value=50, 
                                  value=14, 
                                  key='ori_font_size')
            st.session_state.origami_font_size = font_size
        with col2:
            frame_size = st.slider('Oxview frame size (disable and renable '
                                   'the OxView to apply changes)',
                                    min_value=0, 
                                    max_value=2000, 
                                    value=500, 
                                    key='oxFrame_size')
            st.session_state.oxview_frame_size = frame_size

def simple_origami():
    with st.form(key='simple_origami_form'):
        dt_text = st.text_input("Enter a list of the angles between helices, separated by commas", 
                                "120, ", 
                                key='dt_text_list',
                                help= f"""
The program calculates the best connections bewteen the helices to fit the given angles.
The connection between helices (Dovetails) are obtained roughly with this lookup table (angle --> dt):
{pf.ANGLES_DT_DICT}
""")

        if not dt_text:
            st.stop()
        angle_list = [int(x) for x in dt_text.split(",") if x and x.strip()]
        dt_list = pf.convert_angles_to_dt(angle_list)
        main_stem_default = 11 * ((max([abs(dt) for dt in dt_list], default=0) + 17) // 11 + 1)
        col1, col2 = st.columns(2, vertical_alignment='bottom')
        with col1:
            kl_columns = st.number_input('Kissing loop columns:', min_value=1, value=1, help='number of KL repeats in the helix')
        with col2:
            main_stem = st.number_input('Spacing between crossovers (bp):', min_value=22, value=main_stem_default, step=11, help='The length of the consecutive stems in the helix')

        submitted = st.form_submit_button("Submit")
        if submitted:
            st.session_state.origami = pf.simple_origami(dt_list=angle_list, 
                                                         kl_columns=kl_columns, 
                                                         main_stem=main_stem, 
                                                         add_terminal_helix=True, 
                                                         align=st.session_state.origami.align, 
                                                         use_angles=True)
            st.session_state.code.append(f'origami = pf.simple_origami(dt_list={angle_list}, kl_columns={kl_columns}, main_stem={main_stem}, add_terminal_helix=True, align="{st.session_state.origami.align}", use_angles=True) # Create a simple origami')
            # select the end of the origami
            st.session_state.line_index = len(st.session_state.origami) - 1
            st.session_state.motif_index = len(st.session_state.origami[-1])

def motif_text_format(motif):
    if isinstance(motif, pf.Motif): ### Add the 5' and 3' to the motif text
        # take a motif string and expand it on right, left, top and bottom with spaces
        motif_list = [[' '] * (motif.num_char + 2)] + [[' '] + [char for char in line] + [' '] for line in str(motif).split('\n')] + [[' '] * (motif.num_char + 2)]
        for s in motif:
            if not s.sequence: # skip strand without sequence
                continue
            if s[0] not in '35' and motif_list[s.prec_pos[1] + 1][s.prec_pos[0] + 1] == ' ':
                motif_list[s.prec_pos[1] + 1][s.prec_pos[0] + 1] = '<' + s.directionality[0] + '>'
            if s[-1] not in '35' and motif_list[s.next_pos[1] + 1][s.next_pos[0] + 1] == ' ':
                motif_list[s.next_pos[1] + 1][s.next_pos[0] + 1] = '<' + s.directionality[1] + '>'

        # remove the top and bottom lines if they are empty
        for i in (0, -1):
            if all([char == ' ' for char in motif_list[i]]): 
                motif_list.pop(i)

        motif_str = '\n'.join([''.join(line) for line in motif_list])
    else:
        motif_str = str(motif)
    preview_txt = motif_str.replace(' ', '&nbsp;').replace('\n', '<br />').replace('<5>', '<span style="color: #D52919;">5</span>').replace('<3>', '<span style="color: #D52919;">3</span>')
    return preview_txt

def initiate_session_state():
    if 'origami' not in st.session_state:
        st.session_state.origami = pf.Origami()
    if "code" not in st.session_state:
        st.session_state.code = ["import pyfurnace as pf", "origami = pf.Origami()"]
        st.session_state.motif_buffer = ""
        st.session_state.mod_motif_buffer = ""
    if "motif" not in st.session_state:
        st.session_state.motif = pf.Motif()
    # if "redo" not in st.session_state:
    #     st.session_state.redo = []
    if 'motif_menu_sticky' not in st.session_state:
        st.session_state.motif_menu_sticky = False
    if "copied_motif" not in st.session_state:
        st.session_state.copied_motif = None
    if "copied_motif_text" not in st.session_state:
        st.session_state.copied_motif_text = ""
    if 'modified_motif_text' not in st.session_state:
        st.session_state.modified_motif_text = ""
    if 'upload_key' not in st.session_state:
        st.session_state.upload_key = 0
    ### Initialize the session state variables for the custom motif
    if 'custom_strands' not in st.session_state:
        st.session_state.custom_strands = []
    if 'custom_motifs' not in st.session_state:
        st.session_state.custom_motifs = [('Custom Motif', 'dpad', pf.Motif())]
    if 'custom_key' not in st.session_state:
        st.session_state.custom_key = 0
    if "custom_edit" not in st.session_state:
        st.session_state["custom_edit"] = True
    if "last_clicked_position" not in st.session_state:
        st.session_state["last_clicked_position"] = ''
    if 'line_index' not in st.session_state:
        st.session_state.line_index = 0
    if 'motif_index' not in st.session_state:
        st.session_state.motif_index = 0
    if 'current_line_occupied' not in st.session_state:
        st.session_state.current_line_occupied = False
    if 'origami_click_detector_counter' not in st.session_state:
        st.session_state.origami_click_detector_counter = 0
    if 'max_pk_index' not in st.session_state:
        st.session_state.max_pk_index = 1
    if 'oxview_selected' not in st.session_state:
        st.session_state.oxview_selected = ()
    if 'flip' not in st.session_state:
        st.session_state.flip = False
    if 'selected_motif' not in st.session_state:
        st.session_state.selected_motif = 'Connections'
    if 'gradient' not in st.session_state:
        st.session_state.gradient = False
    if 'oxview_colormap' not in st.session_state:
        st.session_state.oxview_colormap = 'Reds'
    if 'origami_font_size' not in st.session_state:
        st.session_state.origami_font_size = 14
    if 'oxview_frame_size' not in st.session_state:
        st.session_state.oxview_frame_size = 500

def update_file_uploader():
    st.session_state.upload_key += 1

@st.fragment
def make_motif_menu(origami):
    #update Pk index
    st.session_state.max_pk_index = max([abs(int(x.pk_index.replace("'", ""))) for line in origami[lambda m: hasattr(m, 'pk_index')] for x in line], default=0) + 1

    option_data = {'Connections': 'bi-sliders',
                   'Structural': 'bi-bricks',
                   'Kissing Loops': 'bi-heart-half',
                   'Aptamers': 'bi-palette',
                   'Custom': 'bi-joystick',
                   'Edit': 'bi-bandaid',
                   'Code': 'bi-code'}
    
    select_ind = 0
    if st.session_state.selected_motif in option_data:
        select_ind = list(option_data).index(st.session_state.selected_motif)        
    
    selected_motif = option_menu(None, 
                                 list(option_data.keys()),
                                 icons=list(option_data.values()),
                                 menu_icon="cast", 
                                 orientation="horizontal",
                                 manual_select=select_ind,
                                 styles=main_menu_style,
                                 key='motif_menu')
    
    if selected_motif != st.session_state.selected_motif:
        st.session_state.selected_motif = selected_motif
        try:
            st.rerun(scope='fragment')
        except Exception as e:
            print(f"Error in rerun: {e}")
            st.rerun()

    motif_add = True

    if selected_motif == 'Connections':
        ConnectionsCommand().execute()
    elif selected_motif == 'Structural':
        StructuralCommand().execute()
    elif selected_motif == 'Kissing Loops':
        KissingLoopsCommand().execute()
    elif selected_motif == 'Aptamers':
        AptamersCommand().execute()
    elif selected_motif == 'Custom':
        ### Adding the menu here, otherwise there are issues choosing the custom motif with st.fragment when the edit mode is off
        col1, col2 = st.columns([5, 1], vertical_alignment='bottom')
        with col1:
            motif_selected = option_menu(None,
                                        [l[0] for l in st.session_state.custom_motifs],
                                        icons=[f"bi bi-{l[1]}" for l in st.session_state.custom_motifs],
                                        menu_icon="cast", 
                                        orientation="horizontal",
                                        styles=second_menu_style)
        with col2:
            new_name = st.text_input(":green[Add custom motif with name:]", key=f'new_custom_motif{st.session_state.custom_key}')
            if new_name:
                icon = funny_bootstrap_icons[0]
                funny_bootstrap_icons[:] = funny_bootstrap_icons[1:] + [funny_bootstrap_icons[0]]
                st.session_state.custom_motifs.append((new_name, icon, pf.Motif()))
                st.session_state.custom_key += 1
                st.rerun()
        current_custom_motif = [m for m in st.session_state.custom_motifs if m[0] == motif_selected][0][2]
        custom(current_custom_motif)
    elif selected_motif == 'Edit':
        motif_add = False
        edit(st.session_state.motif_index, st.session_state.line_index)
    elif selected_motif == 'Code':
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
        # st.divider()
        # st.markdown("""<hr style="height:1px;border:none;color:#DDDDDD;background-color:#DDDDDD;" /> """, unsafe_allow_html=True)
        st.markdown("<hr style='margin-top:-0em;margin-bottom:-1em' />", unsafe_allow_html=True)
        add_motif(origami)

def select_line(f_col1=None, f_subcol2=None, f_subcol3=None):
    warnings.filterwarnings("error") # raise warnings as errors
    origami = st.session_state.origami
    origami_len = len(origami)
    line_index = st.session_state.line_index
    motif_index = st.session_state.motif_index
    clicked_indexes = st.session_state.get(f"origami_click_detector{st.session_state.origami_click_detector_counter}")

    if line_index < 0:
        line_index = 0
    elif motif_index < 0:
        motif_index = 0

    ### Check if the origami is empty and add an empty line
    if origami_len == 0:
        origami.append([]) # add empty line if the Origami is empty
        st.session_state.code.append('origami.append([]) # Add empty line')
        origami_len = 1
    
    col1, col2, col3, col4, col5, col6 = st.columns([1, 1.5, 1.5, 1, 1, 1], vertical_alignment='bottom')

    with col1:
        if f_col1: f_col1()

    ### Select the line in which to add the motif
    with col2:
        line_index = min(line_index, origami_len - 1)
        line_index = st.number_input("**Select line index:**", min_value=-1, max_value=origami_len - 1, value=line_index, key=f'select_line_index{st.session_state.origami_click_detector_counter}')
        st.session_state.current_line_occupied = False
        if line_index < origami_len and origami[line_index]:
            st.session_state.current_line_occupied = True
        if line_index != st.session_state.line_index:
            st.session_state.line_index = line_index
            if clicked_indexes:
                st.session_state.origami_click_detector_counter += 1
            st.rerun()
        subcol1, subcol2 = st.columns(2)
        
        with subcol1:
            if f_subcol2:
                f_subcol2(line_index + 1)

        ### delete line
        with subcol2:
            if len(st.session_state.origami) > 1 or st.session_state.origami[0]: # check that the origami has more than one line, or the line is not empty
                # If we want to add a new line, we stop the function here
                if st.session_state.line_index >= len(st.session_state.origami):
                    return 
                if len(st.session_state.origami) > 0:
                    del_line = st.button("Delete line", key = 'del_line', help="Delete the choosen line")
                    if del_line:
                        st.session_state.origami.pop(st.session_state.line_index) #remove choosen helix 
                        st.session_state.code.append(f'origami.pop({st.session_state.line_index})')
                        st.session_state.line_index -= 1
                        st.rerun() # rerun app
            # if there is the button to add a motif, add a empty button to align
            elif f_subcol3: 
                st.button('', type='tertiary')

    if st.session_state.line_index >= len(origami):
        return

    ### Motif index selection
    with col3:
        max_val = len(origami[line_index])
        if motif_index > max_val:
            motif_index = max_val
        motif_index = st.number_input("**Select motif index:**", min_value=0, max_value=max_val, value=motif_index, key = f'select_motif_index_{st.session_state.origami_click_detector_counter}')
        if motif_index != st.session_state.motif_index:
            st.session_state.motif_index = motif_index
            if clicked_indexes:
                st.session_state.origami_click_detector_counter += 1
            st.rerun()
        subcol1, subcol2 = st.columns(2)

        with subcol1:
            if f_subcol3:
                f_subcol3()
        with subcol2:
            ### Delete motif
            if len(origami[st.session_state.line_index]) > 0:
                delete_button = st.button(':red[Delete motif]', key = 'delete_motif')
                if delete_button:
                    if st.session_state.motif_index == len(origami[st.session_state.line_index]):
                        st.session_state.motif_index -= 1
                    origami.pop((st.session_state.line_index, st.session_state.motif_index))
                    st.session_state.code.append(f'origami.pop(({st.session_state.line_index}, {st.session_state.motif_index})) # Delete motif')
                    # st.session_state.redo = [] # clear the redo buffer
                    st.rerun() # rerun app

    ### duplicate the current line
    with col4:
        if st.button("Duplicate line", key='duplicate_line'):
            # warnings.filterwarnings("ignore") # ignore kl energy warning
            origami.duplicate_line(line_index, insert_idx=len(origami))
            st.session_state.code.append(f'origami.duplicate_line({line_index}, insert_idx={len(origami)}) # Duplicate line')
            st.session_state.line_index = len(origami) - 1
            st.rerun() # rerun app
    
    ### copy/pase motif
    with col5:
        subcol1, subcol2 = st.columns(2)
        with subcol1:
            copy_motif()
        # Paste motif
        with subcol2:
            paste_button = st.button('Paste motif', key = 'paste_motif')
            if paste_button:
                if not st.session_state.copied_motif:
                    st.warning('No motif copied')
                    return
                origami.insert((line_index, motif_index), st.session_state.copied_motif)
                st.session_state.code.append(st.session_state.copied_motif_text + f'\norigami.insert(({line_index}, {motif_index}), motif) # Paste motif')
                # st.session_state.redo = []
                st.rerun()

    # Uno action
    with col6:
        undo(key='motif_undo')


def add_motif(origami):
    ### Load the motif and prepare the motif buffer for custom motifs
    motif = st.session_state.motif
    if not motif:
        return
    elif type(motif) == pf.Motif: # If the motif is a generic motif, generate the motif text buffer
        ### generate the motif code
        st.session_state.motif_buffer = "strands = []\n"
        for s in st.session_state["custom_strands"]:
            st.session_state.motif_buffer += (f"strands.append(pf.Strand('{s}', "
                                              f"directionality='{s.directionality}', "
                                              f"start={tuple(s.start)}, "
                                              f"direction={tuple(s.direction)}, "
                                              f"coords=pf.Coords({s.coords.array.tolist()}, "
                                              f"dummy_ends=({s.coords.dummy_ends[0].tolist()}, "
                                              f"{s.coords.dummy_ends[1].tolist()}))))\n")
        st.session_state.motif_buffer += f"motif = pf.Motif(strands, structure = '{motif.structure}')"

    ### Show the preview of the motif
    def f_col1():
         # write Preview in small and in orange 
        st.write(':orange[Preview:]')
        # display the Motif to be added next
        # flip the motif if the top direction is 3' to 5' (unless it's a custom motif)
        if st.session_state.flip and motif.__class__.__name__ not in ('Motif', 'start_end_stem'):
            if motif.__class__.__name__ == 'Dovetail': # invert the top and bottom cross for dovetails
                motif.up_cross, motif.down_cross = motif.down_cross, motif.up_cross
            motif.flip(False, True)
        scrollable_text(motif_text_format(motif))

    ### add line
    def f_subcol2(linex_index):
        if len(st.session_state.origami) > 1 or st.session_state.origami[0]: # check that the origami has more than one line, or the line is not empty
            add_line = st.button("Add line", key = 'Insert_line', help=f"Add a new line n° {linex_index}.")
            if add_line:
                origami.insert(linex_index, []) #add empty line at the choosen position
                st.session_state.code.append(f'origami.insert({linex_index}, [])') #add empty line at the choosen position
                st.session_state.line_index += 1
                st.rerun()

    ### Add motif
    def f_subcol3():
            button_label = 'Insert motif'
            if st.session_state.motif_index == len(origami[st.session_state.line_index]):
                button_label = 'Add motif'
            add_button = st.button(f'**:green[{button_label}]**', key = 'insert_motif', help="Insert the motif at the choosen position")
            if add_button:
                origami.insert((st.session_state.line_index, st.session_state.motif_index), motif) # add motif into choosen position
                flip_text = ''
                if st.session_state.flip and motif.__class__.__name__ != 'Motif': 
                    if motif.__class__.__name__ == 'Dovetail': # invert the top and bottom cross for dovetails
                        flip_text += f'\nmotif.up_cross, motif.down_cross = motif.down_cross, motif.up_cross'
                    flip_text += f'\nmotif.flip(False, True) # flip the motif vertically'
                st.session_state.code.append(st.session_state.motif_buffer + f'{flip_text}\norigami.insert(({st.session_state.line_index}, {st.session_state.motif_index}), motif) # Add motif')
                st.session_state.motif_buffer = '' ### clear the buffer
                # st.session_state.redo = [] # clear the redo buffer
                st.session_state.motif_index += 1
                st.rerun()

    select_line(f_col1, f_subcol2, f_subcol3)


def copy_motif(key='', motif=None, motif_slice=None):
    copy_button = st.button('Copy motif', key = f'copy_motif{key}')
    if copy_button:
        if not motif:
            st.session_state.copied_motif = st.session_state.motif
            st.session_state.copied_motif_text = st.session_state.motif_buffer
        else:
            st.session_state.copied_motif = motif
            st.session_state.copied_motif_text = f"motif = origami[{motif_slice}].copy() # Copy motif at line {motif_slice[0]}, index {motif_slice[1]}"

###
# TAB FUNCTIONS
###    

@st.fragment
def generate_custom_motif_text(strand, x_size=50, y_size=10):
    ### update the content
    content = f"<div style='font-family: monospace; font-size: {st.session_state['origami_font_size'] + 8}px;'>"
    current_motif = st.session_state.motif
    m_pos_to_sym = current_motif.get_position_map()
    s_pos_to_sym = strand.get_position_map()
    for y in range(y_size):
        for x in range(x_size):
            if (x, y) in s_pos_to_sym:
                content += f'<a href="javascript:void(0);" id="{x},{y}" style="color: #D00000;">{s_pos_to_sym[(x, y)]}</a>'
            elif (x, y) in m_pos_to_sym: # strand not currently selected
                symbol = m_pos_to_sym[(x, y)]
                content += f'<a href="javascript:void(0);" id="{x},{y}" style="color: #00856A; opacity: 0.5;">{symbol}</a>'
            else:
                content += f'<a href="javascript:void(0);" id="{x},{y}" style="color: grey; opacity: 0.5;">•</a>'
        content += '<br />'
    content += "</div>"
    return content

def custom_text_input(current_custom_motif):
    # make a motif list adding a space around the motif
    motif_lines = str(current_custom_motif).split('\n')
    add_line = False
    add_char = False
    if current_custom_motif:
        if any(line[0] != ' ' for line in motif_lines):
            add_char = True

        motif_list = [[' '] * add_char + [char for char in line] + [' '] for line in motif_lines]

        if motif_lines[0].strip():
            add_line = True
            motif_list = [[' '] * (current_custom_motif.num_char + 2)] + motif_list
        if motif_lines[-1].strip():
            motif_list += [[' '] * (current_custom_motif.num_char + 2)]

    else:
        motif_list = [['']]
    # Add the  5 before the start of each motif.
    for s in current_custom_motif:
        if s.directionality == '53' and motif_list[s.prec_pos[1] + int(add_line)][s.prec_pos[0] + int(add_char)] == ' ': # add 5' to the start of the strand
            motif_list[s.prec_pos[1] + int(add_line)][s.prec_pos[0] + int(add_char)] = '5'
        elif s.directionality == '35' and motif_list[s.next_pos[1] + int(add_line)][s.next_pos[0] + int(add_char)] == ' ': # add 5' to the end of the strand
            motif_list[s.next_pos[1] + int(add_line)][s.next_pos[0] + int(add_char)] = '5'

    current_custom_motif_str = '\n'.join(''.join(line) for line in motif_list)

    strand_text = st.text_area("Motif text: draw a motif where each Strand starts with a 5", value=current_custom_motif_str, key="Motif_text", help='Draw a motif, where each strand has to start with "5". If you want to start a strand with 5, add an additional 5 at the beginning of the strand.')
    if strand_text != current_custom_motif_str:
        if '5' not in strand_text:
            st.warning('Don\'t forget to start the strand with "5"')
        else:
            new_motif = pf.Motif.from_text(strand_text).strip()
            current_custom_motif.replace_all_strands(new_motif._strands, copy=False)
            current_custom_motif.basepair = new_motif.basepair
            st.rerun()

def structure_converter(current_custom_motif):
    ### Add the structure converter
    col1, col2 = st.columns(2)
    with col1:
        structure = st.text_input("Dot-bracket structure:", 
                                 value=current_custom_motif.structure,
                                 key="Dot-bracket_structure", 
                                 help="Add a dot-bracket structure to convert it into a 3D motif.")
    with col2:
        default_seq = str(current_custom_motif.sequence)
        if len(structure) > len(default_seq):
            default_seq += 'N' * (len(structure) - len(default_seq))
        sequence = st.text_input("Sequence:", 
                                 value=default_seq, 
                                 key="Sequence", 
                                 help="Add the sequence of the motif.")
        if not sequence:
            sequence = None
    if (structure and structure != current_custom_motif.structure or
            sequence and sequence != current_custom_motif.sequence):
        try:
            new_motif = pf.Motif.from_structure(structure, sequence=sequence)
        except Exception as e:
            st.error(f"Error: {e}")
            return
        current_custom_motif.replace_all_strands(new_motif._strands, copy=False)
        current_custom_motif.basepair = new_motif.basepair
        st.rerun()

def upload_3d_interface(strand, strand_num, current_custom_motif):
    file_3d = st.file_uploader(f"3D coordinates (OxDNA format) of **strand {strand_num}**", type=['dat', 'pdb'], key=f'custom_{st.session_state.upload_key}', 
                            help='Upload an Oxview configuration ".dat" file with the 3D coordinates of one strand.')
    dummy_cols = st.columns(2)
    with dummy_cols[0]:
        dummy_start = st.toggle('Start dummy nucleotide', key=f'dummy_start_custom', help='The dummy base is a base in the coordinates that is not part of the sequence, but it is used to connect other strands to the beginning of the strand.')
    with dummy_cols[1]:
        dummy_end = st.toggle('End dummy nucleotide', key=f'dummy_end_custom', help='The dummy base is a base in the coordinates that is not part of the sequence, but it is used to connect other strands to the end of the strand.')
    if file_3d:
        pdb_format = file_3d.name.endswith('.pdb')
        strand.coords = pf.Coords.load_from_text(file_3d.getvalue().decode('utf-8'), 
                                                dummy_ends=(dummy_start, dummy_end), 
                                                pdb_format=pdb_format)
        update_strand(strand, strand_num, current_custom_motif, coords=True)

def update_strand(strand, strand_num, current_custom_motif, coords=None):
    st.session_state["custom_strands"][strand_num] = strand
    current_custom_motif.replace_all_strands(st.session_state["custom_strands"], copy=False)
    if coords:
        st.success("3D coordinates uploaded successfully!", icon=":material/view_in_ar:")
        sleep(2)
        update_file_uploader()
    st.rerun()

def strand_number_button(current_custom_motif, delete=True):
    strands = [s.copy() for s in st.session_state["custom_strands"]]
    strand_num = st.radio("Select strand:", list(range(len(strands))), index=len(strands)-1)
    strand = strands[strand_num]

    if delete and st.button("Delete strand", key="delete_strand") and current_custom_motif:
        st.session_state["custom_strands"].pop(strand_num)
        current_custom_motif.pop(strand_num)
        if not st.session_state["custom_strands"]:
            st.session_state["custom_strands"].append(pf.Strand(''))
        st.rerun()
    return strand, strand_num


@st.fragment
def custom(current_custom_motif):
    ### Silence the warnings
    warnings.filterwarnings("ignore")
    st.session_state.motif = current_custom_motif

    if not st.session_state["custom_edit"]:
        if st.button(":orange[Edit the motif]", key="edit_strand"):
            st.session_state["custom_edit"] = True
            st.rerun()
        return
    
    GeneralEditCommand().execute(current_custom_motif)
    
    st.session_state["custom_strands"] = [s.copy() for s in current_custom_motif]
    if not current_custom_motif:
        st.session_state["custom_strands"] = [pf.Strand('')]

    col1, col2 = st.columns([4, 1], vertical_alignment='bottom')
    with col1:
        method = st.segmented_control('Build the motif with:',
                                    ['Structure conveter',
                                    'Drawing Tool', 
                                    'Full text input',
                                    ],
                                    default='Structure conveter',
                                    help='Structure conveter: Convert a dot-bracket structure into a 3D motif automatically.\n\n'
                                    'Drawing Tool: Draw each strand on a grid by clicking on the positions where you want the strand to go;'
                                    'curves and crossings are calculated automatically.\n\n'
                                    'Full text input: Create a motif on a blank text area by typing the motif strand and basepairing.\n'
                                    )
        
    with col2:
        ### add instructions and common symbols
        with st.popover("Common symbols to copy"):
            common_pyroad_sym = ["─", "│", "╭", "╮", "╰", "╯", "^", "┼", '┊', "&"]
            cols = st.columns(len(common_pyroad_sym))
            for i, sym in enumerate(common_pyroad_sym):
                with cols[i]:
                    copy_to_clipboard(sym, sym)

    if method == 'Structure conveter':
        structure_converter(current_custom_motif)
    elif method == 'Full text input':
        custom_text_input(current_custom_motif)
    if method != 'Drawing Tool':

        with st.popover("Upload 3D coordinates"):
            subcol1, subcol2 = st.columns([1, 4])
            with subcol1:
                strand, strand_nr = strand_number_button(current_custom_motif, delete=False)
            with subcol2:
                upload_3d_interface(strand, strand_nr, current_custom_motif)

        if current_custom_motif:
            st.write('Current motif preview:')
            scrollable_text(motif_text_format(current_custom_motif))
            finish_editing(current_custom_motif)
        return

    cols = st.columns(4, vertical_alignment='bottom')
    with cols[0]:
        if st.button("Add strand", key="add_strand"):
            st.session_state["custom_strands"].append(pf.Strand(''))
            current_custom_motif.append(pf.Strand(''), copy=False, join=False)
    with cols[1]:
        x_dots = st.number_input("Canvas x size", min_value=1, value=64, key="x_size")
    with cols[2]:
        y_dots = st.number_input("Canvas y size", min_value=1, value=8, key="y_size")
    with cols[3]:
        if st.button("Clear", key="clear_strand"):
            st.session_state["custom_strands"] = [pf.Strand('')]
            current_custom_motif.replace_all_strands([pf.Strand('')], copy=False)

    # Update and display the motif text
    col1, col2 = st.columns([2, 7])
    with col1:
        strand, strand_num = strand_number_button(current_custom_motif)
    with col2:
        try:
            clicked = click_detector(generate_custom_motif_text(strand, x_size=x_dots, y_size=y_dots), key="custom_motif_click_detector") 
        except Exception as e:
            st.error(str(e))

    def move_strand(clicked, add_nucl=None):
        nonlocal strand
        pos = tuple([int(i) for i in clicked.split(',')])

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
                sym = '-'
                distance = (abs(direction[0]) + 1)
                if direction[0] > 0:
                    strand.direction = (1, 0)
                else:
                    strand.direction = (-1, 0)
            # the Start vertically
            else:
                sym = '|'
                distance = (abs(direction[1]) + 1)
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
                # calulcate the difference between the new direction and the end direction (use to determine the symbol to add)
                normal_direction = tuple([i // abs(i) if i != 0 else 0 for i in direction])
                direction_difference = (normal_direction[0] - strand.end_direction[0], normal_direction[1] - strand.end_direction[1])
                distance = max(abs(pos[0] - strand.end[0]), abs(pos[1] - strand.end[1]))
                # if the direction isn't changed, add the last symbol to match the strand
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
    if clicked and clicked != st.session_state["last_clicked_position"]:
        try:
            move_strand(clicked)
            st.session_state["last_clicked_position"] = clicked
        except pf.MotifStructureError as e:
            st.error(str(e))

    # ### Check the shortcut keys
    # if strand.positions:
    #     next_dir = list(strand.directions)[-1]
    # else:
    #     next_dir = (1, 0)
    # st.markdown("Click on a button to activate the shortcuts!")
    # cols = st.columns([1] * 5)
    # with cols[0]:
    #     button("Move right", "Shift+ArrowRight", shortcut_move, hint=True, args=[(1, 0)])
    # with cols[1]:
    #     button("Move down", "Shift+ArrowDown", shortcut_move, hint=True, args=[(0, 1)])
    # with cols[2]:
    #     button("Move left", "Shift+ArrowLeft", shortcut_move, hint=True, args=[(-1, 0)])
    # with cols[3]:
    #     button("Move up", "Shift+ArrowUp", shortcut_move, hint=True, args=[(0, -1)])
    # with cols[4]:
    #     button("N", f"Shift+N", shortcut_move, hint=True, args=[next_dir, "N"])
    # st.divider()

    ### Show the strand options
    col1, col2, col3, col4, col5 = st.columns([1, 1, 1, 1, 5])
    with col1:
        start_x = st.number_input("Start x:", min_value=0, value=strand.start[0], key=f'start_x_custom')
    with col2:
        start_y = st.number_input("Start y:", min_value=0, value=strand.start[1], key=f'start_y_custom')
    with col3:
        strand_direction_ind = [d for d in pf.Direction].index(strand.direction)
        new_dir = st.selectbox('Start direction:', direction_list, index=strand_direction_ind, key=f'dir_custom')
        new_dir_tuple = pf.Direction[new_dir]
    with col4:
        seq_dir = st.selectbox('Directionality:', ['35', '53'], index=['35', '53'].index(strand.directionality), key=f'seq_dir_custom')
    with col5:
        new_strand = st.text_input(f'New strand (strand directionality: {strand.directionality}) ', value=str(strand), key=f'strand_custom')
    
    ### update the strand
    strand.start = (start_x, start_y)
    strand.direction = new_dir_tuple
    strand.strand = new_strand
    strand.directionality = seq_dir
    with st.popover('Add strand 3D coordinates'):
        upload_3d_interface(strand, strand_num, current_custom_motif)

    # compare the strand with the version strand: if the strand is different, update it and rerun
    old_strand = st.session_state["custom_strands"][strand_num]
    if strand != old_strand:
        update_strand(strand, strand_num, current_custom_motif)

    ### check the base pair symbols of the motif
    current_structure = current_custom_motif.structure
    new_db = st.text_input('Add dot-bracket notation:', value=current_structure,  key=f'structure_custom', help='Add the dot-bracket notation of the motif for each strand, separated by a "&". If the paired bases are more than one position apart, the pairing symbol "┊" is not shown.')
    if new_db != current_structure:
        current_custom_motif.structure = new_db

    finish_editing(current_custom_motif)
    st.stop()

def finish_editing(current_custom_motif):
    if st.button(':green[Finish editing to add motif to the origami]', key="end_edit"):
        st.session_state["custom_edit"] = False
        current_custom_motif.replace_all_strands(st.session_state["custom_strands"], copy=False)
        st.rerun()

def update_code(code_text):
    # Check for forbidden keywords
    forbidden_keywords = ['exec', 'eval', 'compile', ' open', ' os', ' sys', ' __']
    for kw in forbidden_keywords:
        if kw in code_text:
            st.error(f"Forbidden keyword detected in the code: '{kw}'", icon=":material/sentiment_extremely_dissatisfied:")
            return False
    
    # Define the local environment in which the user's code will be executed
    local_context = {'origami': st.session_state.origami}
    # Attempt to execute the code safely
    # try:
    exec(code_text, {'__builtins__': __builtins__, 'pf': pf}, local_context)
    st.session_state.origami = local_context['origami']  # Retrieve the modified origami variable
    st.success("Nanostructure updated successfully!")

    # select the end of the origami
    st.session_state.line_index = len(st.session_state.origami) - 1
    st.session_state.motif_index = len(st.session_state.origami[-1])

    # except Exception as e:
    #     st.error(f"Error in executing the code: {e}")
    #     return False
    st.session_state.code = code_text.split('\n\n')
    st.rerun()

def code():
    if "last_code_id" not in st.session_state:
        st.session_state["last_code_id"] = ''
    code_text = '\n\n'.join(st.session_state.code)

    render_lines = st.slider("Number of lines to render:", min_value=1, max_value=200, value=15, key="render_lines")

    response_dict = code_editor(code_text, 
                                buttons=code_editor_buttons,
                                height=render_lines,
                                options={"showLineNumbers":True}, 
                                allow_reset=True, 
                                key="pyroad_code_editor")
    if response_dict['id'] != st.session_state["last_code_id"] and (response_dict['type'] == "submit" or response_dict['type'] == "selection") :
        st.session_state["last_code_id"] = response_dict['id']
        update_code(response_dict['text'])
    else:
        st.success("Structure updated!")

def undo(key=''):
    if not st.session_state.origami:
        return
    if len(st.session_state.code) <= 1:
        st.warning("Nothing to undo.")
        return
    undo_button = st.button("Undo", 
                            key=f'undo{key}', 
                            use_container_width=False,
                            help="Undo the last action.")
    if not undo_button:
        return
    # st.session_state.redo.append(st.session_state.code[-1])
    update_code('\n\n'.join(st.session_state.code[:-1]))

# def redo():
#     if not st.session_state.redo:
#         # st.warning("Nothing to redo.")
#         return
#     redo_button = st.button("Redo", help="Redo the last action.")
#     if not redo_button:
#         return
#     last_action = st.session_state.redo.pop()
#     update_code('\n\n'.join(st.session_state.code + [last_action]))
#     st.rerun()

def origami_select_display():
    option_data = {"Origami 2D view": "bi bi-square", 
                   "Origami 3D view": "bi bi-box", 
                   "Origami split view": "bi bi-window-split",
                   'Folding barrier': "bi bi-graph-up"}

    selected_display = option_menu(None,
                                list(option_data.keys()),
                                icons=list(option_data.values()),
                                orientation='horizontal',
                                key='DisplayMenu',
                                styles=main_menu_style
                                )
    return selected_display
    

@st.fragment
def origami_build_view(selected_display):
    ### Display the RNA origami structure with clickable elements and modify them in case
    warnings.filterwarnings("ignore") # ignore numpy warnings

    if selected_display == 'Origami 2D view':
        clicked = display_origami()
        clicked_options(clicked)

    elif selected_display == 'Origami 3D view':
        display3d()

    elif selected_display == 'Origami split view':
        col1, col2 = st.columns(2)
        with col1:
            clicked = display_origami()
        with col2:
            display3d()
        clicked_options(clicked)

    elif selected_display == 'Folding barrier':
        barriers, score = st.session_state.origami.folding_barriers()
        clicked = display_origami(barriers=barriers)
        clicked_options(clicked)
        
        help_message = "The folding barriers are a simplified calculation to check " \
        "if a structure is feasible to form co-transcriptionally. They are calculated considering that " \
        "pseudoknots lock the structure in a specific conformation; preventing the formation" \
        "of successive stems. A barrier happens when a stem is open before the formation of a " \
        "pseudoknots, and is closed after the formation of the pseudoknot (after a 150-bases delay). " \
        "The penalty is calculated as W (weak barrier): 1 point, S (strong barrier): 2 points, " \
        "and the total penalty is the sum of the points." 

        cols = st.columns(4, vertical_alignment='center')
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
                if st.button('Try to fix the folding pathway', key='fix_pathway',
                          help="Try to change the starting position of the sequence, "
                          "keeping the same structure, to reduce the folding pathway penalty."
                          "This function is designed for canonical RNA Origami structures, it"
                          "may not work for other complex structures."):
                    # calculate the simplfied dot bracket
                    with st.spinner("Optimizing the folding pathway..."):
                        code_lines = st.session_state.code
                        code_lines += ['\norigami = origami.improve_folding_pathway(kl_delay=150)\n']
                        update_code('\n'.join(code_lines))
                    
def display3d():
    origami = st.session_state.origami
    
    if st.session_state.gradient:
        seq = origami.sequence.replace('&', '')
        index_colors = list(range(len(seq)))

    else:
        m_slice = (st.session_state.line_index, st.session_state.motif_index)
        index_colors = ()
        if m_slice[0] < len(origami) and m_slice[1] < len(origami[m_slice[0]]):
            motif = origami[m_slice[0]][m_slice[1]]
            motif_shift = origami.index_shift_map[m_slice]
            index_colors = [0] * len(origami.sequence.replace('&', ''))
            for pos in motif.seq_positions:
                shifted = tuple([pos[0] + motif_shift[0], pos[1] + motif_shift[1]])
                index_colors[origami.seq_positions.index(shifted)] = 1

    if index_colors:
        for s in origami.strands:
            for protein in s.coords.proteins:
                index_colors += [0] * len(protein)

    conf, topo = origami.save_3d_model("origami", return_text=True)
    oxview_from_text(configuration=conf, # path to the configuration file
                    topology=topo,      # path to the topology file
                    width='99%',                # width of the viewer frame
                    colormap=st.session_state.oxview_colormap, # colormap for the viewer
                    height=st.session_state.oxview_frame_size, # height of the viewer frame
                    index_colors=index_colors, # color the bases in the viewer
                    key='display_nano')  

def build_origami_content(barriers=None):
    barriers_colors ={'▂':'#FFBA08', '▄':'#FFBA08', '█':'#D00000'}
    highlight_color = '#D00000'
    normal_color = '#333333' # 80% black

    origami = st.session_state.origami
    motif = origami.assembled

    if barriers:
        origami_lines = origami.barrier_repr(barriers=barriers, 
                                             return_list=True)
    else:
        if st.session_state.get('to_road'):
            origami_str = origami.to_road()

        else:
            origami_str = str(origami)
        origami_lines = origami_str.split('\n')

        
    # create color gradient
    if st.session_state.gradient:
        # create a dictionary from positions to index
        pos_to_index = {pos: ind for ind, pos in enumerate(motif.seq_positions)}
        tot_len = 0
        for s in origami.strands:
            tot_len += len(s.sequence)
            for protein in s.coords.proteins:
                tot_len += len(protein)
        cmap = plt.get_cmap(st.session_state.oxview_colormap)
        oxview_offset = tot_len // 5
        c_map = [mcolors.to_hex(cmap(i)) for i in np.linspace(0, 1, tot_len + oxview_offset)]
        c_map = c_map[oxview_offset:]
    
    # Prepare the string to add 5' and 3' symbols for the strands
    motif_list = [[' '] * (motif.num_char + 2)] + [[' '] + [char for char in line] + [' '] for line in origami_lines] + [[' '] * (motif.num_char + 2)]
    for s in motif:  # Add the 5' and 3' symbols to the motif as 1 and 2
        if not s.sequence:
            continue
        if s.sequence and s[0] not in '35' and motif_list[s.prec_pos[1] + 1][s.prec_pos[0] + 1] == ' ':
            motif_list[s.prec_pos[1] + 1][s.prec_pos[0] + 1] = '1' if s.directionality == '53' else '2'
        if s.sequence and s[-1] not in '35' and motif_list[s.next_pos[1] + 1][s.next_pos[0] + 1] == ' ':
            motif_list[s.next_pos[1] + 1][s.next_pos[0] + 1] = '2' if s.directionality == '53' else '1'
    origami_list = [''.join(line) for line in motif_list]

    content = "<div style='font-family: monospace; white-space: nowrap; overflow-x: auto;'>"
    content += f'<div style="display:inline-block; font-family: monospace; font-size: {st.session_state.origami_font_size}px;">Line:<br />'  # add a column with the line number
    line_nr = -2
    origami_list_len = len(origami_list)
    span_text = '<span style="font-family: monospace; '

    for y, line in enumerate(origami_list):
        line_color = normal_color
        motif_slice = None
        current_line_nr = [origami.pos_index_map[(x-1, y-1)] 
                                for x, _ in enumerate(line) 
                                    if (x-1, y-1) in origami.pos_index_map]
        
        next_line = origami_list[y+1] if  (y + 1) < origami_list_len else []
        
        # check the next line index
        next_line_nr = [origami.pos_index_map[(x-1, y)] 
                                for x, _ in enumerate(next_line) 
                                        if (x-1, y) in origami.pos_index_map] 

        if current_line_nr: 
            current_line_nr = current_line_nr[0][0]  # get the first line number
        if next_line_nr:
            next_line_nr = next_line_nr[0][0]

        if line_nr != current_line_nr and isinstance(current_line_nr, int):  # is a new origami line
            if current_line_nr == st.session_state.line_index:
                line_color = highlight_color
            content += span_text + f'color: {line_color}; line-height:1;">{current_line_nr})</span>' + (4 - len(str(current_line_nr))) * '&nbsp;' 
            line_nr = current_line_nr
        elif  isinstance(next_line_nr, int) and line_nr < st.session_state.line_index < next_line_nr:
            content += span_text + 'color: #D52919; line-height:1;">_____</span>'  # the origami line is empty
        else:
            content += span_text + 'line-height:1;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span>'  # is not a new origami line
        
        hit = False # to highlight the first symbol
        for x, char in enumerate(line):
            ori_pos = (x - 1, y - 1)
            color = normal_color
            if barriers_colors and char in barriers_colors:
                color = barriers_colors[char]

            if char == ' ':
                content += span_text + 'line-height:1;">&nbsp;</span>'
            elif char == '1':
                content += span_text + f'color: {highlight_color}; line-height:1;">5</span>'
            elif char == '2':
                content += span_text + f'color: {highlight_color}; line-height:1;">3</span>'
            elif char in pf.bp_symbols:
                content += span_text + f'color: {color}; line-height:1;">{char}</span>'  # do not highlight the base pair in red
            elif ori_pos in origami.pos_index_map:  # a motif symbol
                motif_slice = origami.pos_index_map[ori_pos]
                if st.session_state.gradient: 
                    index = pos_to_index.get(ori_pos)
                    if index is not None:
                        color = c_map[index]

                # This is the selected motif
                if motif_slice and motif_slice[0] == st.session_state.line_index and motif_slice[1] == st.session_state.motif_index:
                    if st.session_state.gradient: # Don't add color if the gradient is not active
                        color = normal_color
                    elif not hit: # the first symbol is yellow
                        color = '#FF8800'
                        hit = True
                    else:
                        color = highlight_color

                content += f'<a style="font-family: monospace; color: {color}; line-height:1;" href="javascript:void(0);" id="{motif_slice[0]},{motif_slice[1]},{x - 1},{y - 1}">{char}</a>'
            else:  # is a junction symbol 
                content += span_text + f'color: {color}; line-height:1;">{char}</span>'
        # Check if you wanna add the cursor
        if motif_slice and motif_slice[0] == st.session_state.line_index:
            if len(origami[motif_slice[0]]) == st.session_state.motif_index:
                content += span_text + f'color: {highlight_color}; line-height:1;">│</span>'
        content += '<br />'
    if line_nr < st.session_state.line_index:
        content += span_text + 'color: #D52919; line-height:1;">_____</span>'  # the origami line is empty
    content += '</div>'
    content += '</div>'
    return content


def display_origami(barriers=None):
    ### SHORTCUTS NOT READY YET!!!
    # if not split_view:
    #     def move_selection(x = 0, y = 0):
    #         st.session_state.line_index += y
    #         st.session_state.motif_index += x
    #         # st.rerun()

    #     with st.expander("Shortcuts to move over the Origami", expanded=False):
    #         st.write("Use the arrow keys to move the selection over the Origami lines and motifs. Not available in split view.")
    #         col1, col2, col3, col4 = st.columns(4)
    #         with col1:
    #             streamlit_shortcuts.button("Move Up", on_click=move_selection, shortcut="Shift+ArrowUp", hint=True, kwargs={"y" : -1})
    #         with col2:
    #             streamlit_shortcuts.button("Move Down", on_click=move_selection, shortcut="Shift+ArrowDown", hint=True, kwargs={"y" : 1})
    #         with col3:
    #             streamlit_shortcuts.button("Move Left", on_click=move_selection, shortcut="Shift+ArrowLeft", hint=True, kwargs={"x" : -1})
    #         with col4:
    #             streamlit_shortcuts.button("Move Right", on_click=move_selection, shortcut="Shift+ArrowRight", hint=True, kwargs={"x" : 1})

    try:
        content = build_origami_content(barriers=barriers)
        clicked = click_detector(content, key=f"origami_click_detector{st.session_state.origami_click_detector_counter}")
    except pf.MotifStructureError as e:
        st.error(f"Structure error: \n {e}", icon=":material/personal_injury:")
        undo('error', key='error')
        st.stop()
    except pf.AmbiguosStructure as e:
        st.error(f"Ambigouse structure: \n {e}", icon=":material/theater_comedy:")
        undo('warning', key='warning')
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
        y_slice, x_slice, x_pos, y_pos  = [int(i) for i in clicked.split(',')]
        motif_slice = (y_slice, x_slice)
        pos = (x_pos, y_pos)

        # Write the selected motif on screen
        try:
            motif = st.session_state.origami[motif_slice]
        except IndexError:
            return # this is a problem with updating the keys
        
        motif_class_name = motif.__class__.__name__
        if motif_class_name == 'Motif': ### is custom Motif
            motif_class_name = "Custom " + motif_class_name

        nucl_text = ''
        # Indicate the nucleotide index:
        if pos in st.session_state.origami.seq_positions:
            ind = st.session_state.origami.seq_positions.index(pos)
            nucl_text += f':orange[nucleotide index {ind}] in'

        ### highlight the selected motif
        if st.session_state.line_index != motif_slice[0] or st.session_state.motif_index != motif_slice[1]:
            st.session_state.line_index = motif_slice[0]
            st.session_state.motif_index = motif_slice[1]
            st.rerun()

        col1, col2, _ = st.columns([1, 1, 1], vertical_alignment='center')
        with col1:
            st.markdown(f'Selected {nucl_text} :orange[{motif_class_name}]: :orange[line {motif_slice[0]}], :orange[motif {motif_slice[1]}]')
        with col2:
            ### copy the motif
            copy_motif('selected', motif=motif, motif_slice=motif_slice)

def display_structure_sequence():
    st.write("#### Origami structure and sequence constraints:")
    origami = st.session_state.origami
    if origami.structure:
        structure_list = [char for char in origami.structure]
        sequence_list = [char for char in origami.sequence]
        motif_slice = (st.session_state.line_index, st.session_state.motif_index)
        indexes = []
        if motif_slice in origami.index_shift_map:
            motif = origami[motif_slice]
            shifts_x, shift_y = origami.index_shift_map[motif_slice]
            for pos in motif.seq_positions:
                indexes.append(origami.seq_positions.index((pos[0] + shifts_x, pos[1] + shift_y)))

        for i in indexes:
            structure_list[i] = '<span style="color: #D52919; line-height:1;">' + structure_list[i] + "</span>"
            sequence_list[i] = '<span style="color: #D52919; line-height:1;">' + sequence_list[i] + "</span>"
            # structure_list[i] = f":orange[{structure_list[i]}]"
            # sequence_list[i] = f":orange[{sequence_list[i]}]"

        ### Pseudoknots info
        pseudoknot_text = '; '.join([f'id: {k}, ' + str(val) 
                                        for k, val in origami.pseudoknots.items()]) + ';'
        remove_syms = str.maketrans('', '', "{}'\"")
        pseudoknot_text = pseudoknot_text.translate(remove_syms)

        st.markdown(f'**Structure length: :green[{str(len(origami.structure.replace("&", "")))}]**')
        scrollable_text(f"""
            >Origami</br>{''.join(sequence_list)}</br>{''.join(structure_list)}</br></br>Pseudoknots info:</br>{pseudoknot_text}""")

        st.markdown("<hr style='margin-top:-0em;margin-bottom:-1em;color:white;border-color:white' />", unsafe_allow_html=True)

        st.session_state.generate_structure = origami.structure
        st.session_state.generate_sequence = str(origami.sequence)
        st.session_state.generate_pseudoknots = pseudoknot_text
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.page_link("pages/2_Generate.py", label="**:orange[Generate the sequence]**", icon=":material/network_node:")
        with col2:
            copy_to_clipboard(origami.structure, 'Structure')
        with col3:
            copy_to_clipboard(origami.sequence, 'Sequence')
        with col4:
            copy_to_clipboard(pseudoknot_text, 'Pseudoknots')

def edit(x, y):
    origami = st.session_state.origami
    ### select the motif to change
    st.session_state.modified_motif_text = '' # initialize the modification
    try:
        motif_slice = (y, x)
        if y >= len(origami) or x >= len(origami[y]):
            select_line()
            return
        motif = origami[motif_slice]
    except KeyError:
        return # this is a problem with updating the keys
    motif_class_name = motif.__class__.__name__
    if motif_class_name == 'Motif': ### is custom Motif
        motif_class_name = "Custom " + motif_class_name
    elif isinstance(motif, pf.KissingLoop):
        motif_class_name = 'KissingLoops'

    ### Directly modify the motif and make a copy for the advanced feature
    if not st.session_state.flip: ### if the top strand direction is '35', the modifications don't work well
        command = f"{motif_class_name}Command"
        if command in globals():
            globals()[command]().execute(motif)
        else:
            GeneralEditCommand().execute(motif)
    else:
        st.warning("Attribute modifications are not available when the top strand direction is '35'.")

    if st.session_state.modified_motif_text:
        st.session_state.code.append(f'motif = origami[({motif_slice[0]}, {motif_slice[1]})] # motif slice: line {motif_slice[0]}, index {motif_slice[1]}' + st.session_state.modified_motif_text)
        st.session_state.modified_motif_text = ''
        st.rerun()

    select_line()
    advaced_edit(motif, motif_slice)

@st.fragment
def advaced_edit(motif, motif_slice):
    motif_copy = motif.copy()
    y, x = motif_slice

    ### try to change each strand
    with st. container():
        with st.popover('Advanced feature, modify the strands of the selected motif:',
                        use_container_width=True,):
            
            col1, col2 = st.columns(2)
            with col1:
                x_shift = st.number_input("Motif horizontal shift:", value=0)
            with col2:
                y_shift = st.number_input("Motif vertical shift:", value=0)

            shifted = False
            if x_shift != 0 or y_shift != 0:
                shifted = True
                motif_copy.shift((x_shift, y_shift))

            for i, s in enumerate(motif_copy):
                col1, col2, col3, col4, col5 = st.columns([1, 1, 1, 1, 5])
                with col1:
                    start_x = st.number_input("Start x:", min_value=0, value=s.start[0], key=f'start_x_{x}_{y}_{i}')
                with col2:
                        start_y = st.number_input("Start y:", min_value=0, value=s.start[1], key=f'start_y_{x}_{y}_{i}')
                with col3:
                    strand_direction_ind = [d for d in pf.Direction].index(s.direction)
                    new_dir = st.selectbox('Start direction:', direction_list, index=strand_direction_ind, key=f'dir_{x}_{y}_{i}')
                    new_dir_tuple = pf.Direction[new_dir]
                with col4:
                    seq_dir = st.selectbox('Directionality:', ['35', '53'], index=['35', '53'].index(s.directionality), key=f'seq_dir_{x}_{y}_{i}')
                with col5:
                    new_strand = st.text_input(f'New strand (strand directionality: {s.directionality}) ', value=str(s), key=f'strand_{x}_{y}_{i}')
                ### check if the strand is modified and update the motif
                stripped = new_strand.strip()
                if s.start != (start_x, start_y):
                    s.start = (start_x, start_y)
                    st.session_state.modified_motif_text += f'\nmotif[{i}].start = ({start_x}, {start_y})'
                if s.directionality != seq_dir:
                    s.directionality = seq_dir
                    st.session_state.modified_motif_text += f'\nmotif[{i}].directionality = "{seq_dir}"'
                if s.direction != new_dir_tuple:
                    s.direction = new_dir_tuple
                    st.session_state.modified_motif_text += f'\nmotif[{i}].direction = {new_dir_tuple}'
                if s.strand != stripped:
                    s.strand = stripped
                    st.session_state.modified_motif_text += f'\nmotif[{i}].strand = "{stripped}"'

            ### check the base pair symbols of the motif
            current_structure = motif_copy.structure
            new_db = st.text_input('Add dot-bracket notation:', value=current_structure,  key=f'structure_{x}_{y}', help='Add the dot-bracket notation of the motif for each strand, separated by a "&". If the paired bases are more than one position apart, the pairing symbol "┊" is not shown.')
            if new_db != current_structure:
                if not new_db:
                    motif_copy.autopairing = True
                    st.session_state.modified_motif_text += f'\nmotif.autopairing = True'
                else:
                    motif_copy.structure = new_db
                    st.session_state.modified_motif_text += f'\nmotif.structure = "{new_db}"'
            
            ### update the motif
            try:
                st.write(':orange[Preview:]')

                scrollable_text(motif_text_format(motif_copy))
                if not st.session_state.modified_motif_text and not shifted:
                    st.write(':green[Structure is updated]') 
                else:
                    st.warning('Be careful before updating, there is a high risk of breaking the structure.')
                    update = st.button('Update strands')
                    if update:
                        if shifted:
                            st.session_state.modified_motif_text += f'\nmotif.shift(({x_shift}, {y_shift}))'
                        st.session_state.code.append(f'motif = origami[({motif_slice[0]}, {motif_slice[1]})] # select the motif at index, line: {motif_slice}' + st.session_state.modified_motif_text)
                        st.session_state.origami[motif_slice] = motif_copy
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
            font-size: {st.session_state.origami_font_size + 1}px;
        }}
        </style>
        <div class="scroll-box">{text}</div>
        """, unsafe_allow_html=True)