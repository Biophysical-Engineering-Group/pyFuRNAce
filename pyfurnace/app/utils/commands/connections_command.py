import inspect
from pyfurnace.design import utils, start_end_stem  # Import the module with the variables

import streamlit as st
from streamlit_option_menu import option_menu
from streamlit_option_menu import option_menu
from .motif_command import MotifCommand
from .general_edit_command import GeneralEditCommand
from .tetraloop_command import TetraLoopCommand
from .. import second_menu_style
from ..motifs_icons import MOTIF_ICONS

# Filter and collect the motif utils
util_names = [ut_name for ut_name, obj in inspect.getmembers(utils.motif_lib) if inspect.isfunction(obj)]

# ignore the simple vertical connection
del util_names[util_names.index('vertical_link')]

# Add the tetraloop option
util_names = ['Tetraloop'] + util_names

class ConnectionsCommand(MotifCommand):
    
    def execute(self, motif=None):
        ### Modify motif
        if motif:
           GeneralEditCommand().execute(motif)

        ### Create new motif
        else:
            util_option = option_menu(
                                    None, 
                                    util_names,
                                    icons=[MOTIF_ICONS[name] for name in util_names],
                                    menu_icon="cast",  
                                    orientation="horizontal",
                                    styles=second_menu_style,
                                    key='UtilsOption',
                                    )

            if util_option == 'Tetraloop':
                TetraLoopCommand().execute()
                return
            if util_option == 'start_end_stem':
                start_end_stemCommand().execute()
                return
            
            name = util_option.replace(' ', '_')
            motif_util = getattr(utils.motif_lib, name)
            motif_text = f"motif = pf.{name}"
            flip_vert, flip_hor, rotate = self.interface()
            if flip_vert or flip_hor:
                motif_text += f".flip({flip_vert}, {flip_hor})"
            if rotate:
                motif_text += f".rotate({rotate})"
            # save the motif in the session state
            st.session_state.motif_buffer = motif = f"motif = pf.{name}(hflip={flip_hor}, vflip={flip_vert}, rotate={rotate})"
            st.session_state.motif = motif_util(hflip=flip_hor, vflip=flip_vert, rotate=rotate)
        

    def interface(self, key=''):
        return GeneralEditCommand.interface(key)
    

class start_end_stemCommand(MotifCommand):

    def execute(self, motif=None):
        kwargs = dict()
        if st.session_state.flip:
            kwargs = {'top_l_def': '5', 'top_r_def': '3', 'bot_l_def': '3', 'bot_r_def': '5', 'top_ind': 1, 'bot_ind': 0}
        t_l, t_r, b_l, b_r = self.interface(**kwargs)
        st.session_state.motif_buffer = f"motif = pf.start_end_stem(top_left='{t_l}', top_right='{t_r}', bot_left='{b_l}', bot_right='{b_r}')"
        st.session_state.motif = start_end_stem(top_left=t_l, top_right=t_r, bot_left=b_l, bot_right=b_r)


    def interface(self, key='', top_l_def='3', top_r_def='5', bot_l_def='5', bot_r_def='3', top_ind = 0, bot_ind = 1):
        _, _, col1, col2, _, _ =st.columns([1, 1, 1, 1, 1, 1])
        with col1:
            t_l = st.selectbox('Top left:', [top_l_def] + ['─', None], index=top_ind, key=f'start_end_stem_top_left{key}')
            b_l = st.selectbox('Bottom left:', [bot_l_def] +['─', None], index=bot_ind, key=f'start_end_stem_bot_left{key}')
        with col2:
            t_r = st.selectbox('Top right:', [top_r_def] + ['─', None], index=top_ind, key=f'start_end_stem_top_right{key}')
            b_r = st.selectbox('Bottom right:', [bot_r_def] + ['─', None], index=bot_ind, key=f'start_end_stem_bot_right{key}')
        if not t_l: t_l = ''
        if not t_r: t_r = ''
        if not b_l: b_l = ''
        if not b_r: b_r = ''
        return t_l, t_r, b_l, b_r