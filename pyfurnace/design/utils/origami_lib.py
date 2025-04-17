from IPython.display import display, HTML
try:
    from oxDNA_analysis_tools.UTILS.oxview import oxdna_conf
    from oxDNA_analysis_tools.UTILS.RyeReader import describe, get_confs
    oat_installed = True
except:
    oat_installed = False
import tempfile
from ..core import *
from ..motifs import *
from .motif_lib import *

angles_dt_dict = {26: -6,
                  58: - 5,
                  90 : -4,
                  122 : -3,
                  122: -3,
                  154: -2,
                  186: -1,
                  218: 0,
                  250: 1,
                  282: 2,
                  314: 3,
                  346: 4,
                  378: 5,
                  410: 6}

def convert_angles_to_dt(angles_list: list):
    angles_sanitize = [ang % 360 for ang in angles_list]
    # get the closest angle in the dict
    dt_list = [angles_dt_dict[min(angles_dt_dict, key=lambda x:abs(x-ang))] for ang in angles_sanitize]
    return dt_list
    

def simple_origami(dt_list: list, helix_kl: int = 1, main_stem: list = None, left_stem_kl: list = None, stem_pos: list = None, start: int=0, add_terminal_helix: bool = True, align: str = 'first',
                   end_helix_len = 8, use_angles=False, add_start_end = True) -> Origami:
    """ Build an origami structure with the given dovetail list.
    Use the helix_kl parameter to specify the number of KL repeats in the helix.

    Args:
        dt_list (list): 
            list of dovetail values
        helix_kl (int): 
            number of KL repeats in the helix
        main_stem (list): 
            List of the length of the consecutive stems for each KL (by default, is the minimum value for each dovetail)
        left_stem_kl (list):
            List of the length of the left stem for each KL (by default, the half of the main stem - abd(dt) - 8 (for KL))
        stem_pos (list):
            List of the helix position of the main stem at each kissing loop (by default at the first helix)
        start (int):
            The index of the main stem to start the origami (by default at the first main stem)
        add_terminal_helix (bool):
            Add a terminal helix at the start and end of the dovetail list (a helix with 0 dovetail)
        end_helix_len (int):
            Length of the end helix (default 8, the minimum length for a stable helix)
        use_angles (bool):
            Use helices angles insted of the dove-tail values
        ... to finish
        

    Returns:
        Origami: the origami structure built with the given parameters
    """

    # initialize the origami structure
    origami = Origami(align=align)

    if use_angles:
        dt_list = convert_angles_to_dt(dt_list)
    
    # add the start and end helix to the dovetail list
    if add_terminal_helix:
        dt_list = [0] + dt_list + [0]

    # if the main_stem list is not given, set it to the minimum value for each KL
    if main_stem is None:
        max_dt = max([abs(dt) for dt in dt_list], default=0)
        main_stem = [[11 * ((max_dt + 17) // 11 + 1)] * helix_kl] * len(dt_list)
    elif type(main_stem) == int:
        main_stem = [[main_stem for _ in range(helix_kl)] for _ in range(len(dt_list))]
    elif type(main_stem) == list and all(isinstance(x, int) for x in main_stem):
        main_stem = [main_stem for _ in range(len(dt_list))]
    elif type(main_stem) == list and all(isinstance(x, (tuple, list)) for x in main_stem):
        if not all(len(x) == helix_kl for x in main_stem):
            raise ValueError("The main_stem list should have the same length as the kissing loops repeats")
    else:
        raise ValueError("The main_stem can be an int, a list of int or a matrix of int")

    if left_stem_kl is None:
        left_stem_kl = [[None] * helix_kl for _ in range(len(dt_list))]
    elif type(left_stem_kl) == int:
        left_stem = [[left_stem_kl for _ in range(helix_kl)] for _ in range(len(dt_list))]
    elif type(left_stem_kl) == list and all(isinstance(x, int) for x in left_stem_kl):
        left_stem_kl = [[left_stem_kl[i]] * helix_kl for i in range(len(dt_list))]
    elif type(left_stem_kl) == list and all(isinstance(x, (tuple, list)) for x in left_stem_kl):
        if not all(len(x) == helix_kl for x in left_stem_kl):
            raise ValueError("The left_stem_kl list should have the same length as the kissing loops repeats")
    else:
        raise ValueError("The left_stem_kl can be an int, a list of int or a matrix of int")

    if stem_pos is None:
        stem_pos = [0 for _ in range(helix_kl)]
    elif type(stem_pos) == int:
        stem_pos = [stem_pos for _ in range(helix_kl)]

    # create an helix for each dovetail in the list
    for helix_in, dt in enumerate(dt_list):

        # create the start of the stem: a tetraloop and a stem of 5 bases
        helix = [TetraLoop(), Stem(end_helix_len), Dovetail(dt)]

        # add Kissing loops repeats to the helix
        for kl_index in range(helix_kl):
            stem_len = main_stem[helix_in][kl_index]
            left_stem = left_stem_kl[helix_in][kl_index]
            if left_stem is None:
                left_stem = (stem_len - 8 - abs(dt)) // 2
            right_stem = (stem_len - 8 - abs(dt)) - left_stem

            # if the helix position is in the stem_position list for the given KL index, add a stem
            if stem_pos[kl_index] == helix_in:
                if kl_index == start and add_start_end: # add the start motif after the first stem
                    half_l_stem = (stem_len - abs(dt)) // 2
                    half_r_stem = stem_len - abs(dt) - half_l_stem
                    helix += [Stem(half_l_stem)
                                .shift((1,0), extend=True), 
                              start_end_stem(), 
                              Stem(half_r_stem), Dovetail(dt)]
                else:
                    helix += [Stem(main_stem[helix_in][kl_index] - abs(dt))
                                .shift((6,0), extend=True),
                              Dovetail(dt)]
            # add a kissing normal loop repeat
            else:
                helix += [Stem(left_stem), KissingDimer(), Stem(right_stem), Dovetail(dt)]

        # add the end of the helix: a stem of 5 bases and a tetraloop
        helix += [Stem(end_helix_len), TetraLoop(open_left=True)]
        # add the helix to the origami
        origami.append(helix)

    # remove the top cross from the dovetails of the first helix
    for motif in origami[0]:
        if type(motif) == Dovetail:
            motif.up_cross = False

    # remove the bottom cross from the dovetails of the last helix
    for motif in origami[-1]:
        if type(motif) == Dovetail:
            motif.down_cross = False

    # return the origami structure
    return origami

def ipython_display_3D(origami, **kwargs):
    if not oat_installed:
        warnings.warn("The oxDNA_analysis_tools package is not installed, the 3D display is not available.")
        return
    # Create a temporary directory
    with tempfile.TemporaryDirectory() as tmpdirname:
        file_path = f"{tmpdirname}/origami"
        origami.save_3d_model(file_path)
        top_info, traj_info = describe(f'{file_path}.top', f'{file_path}.dat')
        conf = get_confs(top_info, traj_info, 0, 1)[0]
        oxdna_conf(top_info, conf, **kwargs)

def ipython_display_txt(origami_text, max_height='500'):
    # Convert your text to scrollable HTML
    scrollable_html = f"""
    <div style="max-height: {max_height}px; overflow-y: scroll; border: 1px solid #ccc; padding: 10px;">
    <pre>{origami_text}</pre>
    </div>
    """
    display(HTML(scrollable_html))
