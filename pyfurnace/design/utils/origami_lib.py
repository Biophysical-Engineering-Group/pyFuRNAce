from IPython.display import display, HTML
from typing import List, Optional, Union, Any
import tempfile
# try to import the oxDNA_analysis_tools package
try:
    from oxDNA_analysis_tools.UTILS.oxview import oxdna_conf
    from oxDNA_analysis_tools.UTILS.RyeReader import describe, get_confs
    oat_installed = True
except:
    oat_installed = False
# own imports
from ..core import *
from ..motifs import *
from .motif_lib import *

# A dictionary to convert angles to dovetail values
ANGLES_DT_DICT = {26: -6,
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

def convert_angles_to_dt(angles_list: List[float]) -> List[int]:
    """
    Convert a list of helix angles into corresponding dovetail values based 
    on a predefined mapping.

    Parameters
    ----------
    angles_list : list of float
        List of helix angles in degrees. Angles will be wrapped modulo 360.

    Returns
    -------
    list of int
        Corresponding dovetail values for each angle in the input list.
    """
    angles_sanitize = [ang % 360 for ang in angles_list]
    # get the closest angle in the dict
    dt_list = [ANGLES_DT_DICT[min(ANGLES_DT_DICT, 
                                key=lambda x:abs(x-ang))] 
                                    for ang in angles_sanitize]
    return dt_list
    

def simple_origami(
    dt_list: List[int],
    kl_columns: int = 1,
    main_stem: Optional[Union[int, List[int], List[List[int]]]] = None,
    left_stem_kl: Optional[Union[int, List[int], List[List[int]]]] = None,
    stem_pos: Optional[Union[int, List[int]]] = None,
    start: int = 0,
    add_terminal_helix: bool = True,
    end_helix_len: int = 8,
    use_angles: bool = False,
    add_start_end: bool = True,
    align: str = 'first'
    ) -> Origami:
    """
    Construct an RNA origami object based on a sequence of dovetail values and 
    kissing loop parameters.

    Parameters
    ----------
    dt_list : list of int
        List of dovetail values representing inter-helix connections.
    kl_columns : int, optional
        Number of kissing loop repeats in each helix (default is 1).
    main_stem : int or list of int or list of list of int, optional
        Length(s) of the main stem in each kissing loop.
        Can be a single int, a list (same for all loops), or a matrix for per-loop
          customization.
    left_stem_kl : int or list of int or list of list of int, optional
        Length(s) of the left stem for each kissing loop. Defaults to automatic
          computation.
    stem_pos : int or list of int, optional
        Position(s) of the main stem insertion among helices. Default is 0 for all.
    start : int, optional
        Index of the main stem where origami building starts (default is 0).
    add_terminal_helix : bool, default True
        Whether to prepend and append helices with no dovetails.
    end_helix_len : int, optional
        Length of the stems at the ends of the helices (default is 8).
    use_angles : bool, optional
        If True, interpret `dt_list` as helix angles and convert them to dovetail 
        values (default is False).
    add_start_end : bool, default True
        Whether to add a start-end motif in the initial helix.
    align : str, optional
        Alignment method for the origami object (default is 'first').

    Returns
    -------
    Origami
        The assembled Origami structure.
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
        main_stem = [[11 * ((max_dt + 17) // 11 + 1)] * kl_columns] * len(dt_list)
    elif type(main_stem) == int:
        main_stem = [[main_stem for _ in range(kl_columns)] for _ in range(len(dt_list))]
    elif type(main_stem) == list and all(isinstance(x, int) for x in main_stem):
        main_stem = [main_stem for _ in range(len(dt_list))]
    elif type(main_stem) == list and all(isinstance(x, (tuple, list)) for x in main_stem):
        if not all(len(x) == kl_columns for x in main_stem):
            raise ValueError("The main_stem list should have the same length as the kissing loops repeats")
    else:
        raise ValueError("The main_stem can be an int, a list of int or a matrix of int")

    if left_stem_kl is None:
        left_stem_kl = [[None] * kl_columns for _ in range(len(dt_list))]
    elif type(left_stem_kl) == int:
        left_stem = [[left_stem_kl for _ in range(kl_columns)] for _ in range(len(dt_list))]
    elif type(left_stem_kl) == list and all(isinstance(x, int) for x in left_stem_kl):
        left_stem_kl = [[left_stem_kl[i]] * kl_columns for i in range(len(dt_list))]
    elif type(left_stem_kl) == list and all(isinstance(x, (tuple, list)) for x in left_stem_kl):
        if not all(len(x) == kl_columns for x in left_stem_kl):
            raise ValueError("The left_stem_kl list should have the same length as the kissing loops repeats")
    else:
        raise ValueError("The left_stem_kl can be an int, a list of int or a matrix of int")

    if stem_pos is None:
        stem_pos = [0 for _ in range(kl_columns)]
    elif type(stem_pos) == int:
        stem_pos = [stem_pos for _ in range(kl_columns)]

    # create an helix for each dovetail in the list
    for helix_in, dt in enumerate(dt_list):

        # create the start of the stem: a tetraloop and a stem of 5 bases
        helix = [TetraLoop(), Stem(end_helix_len), Dovetail(dt)]

        # add Kissing loops repeats to the helix
        for kl_index in range(kl_columns):
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

def ipython_display_3D(origami: Origami, **kwargs: Any) -> None:
    """
    Display a 3D representation of an Origami structure within a J
    upyter notebook using oxDNA.

    Parameters
    ----------
    origami : Origami
        The Origami structure to visualize.
    **kwargs : dict, optional
        Additional keyword arguments passed to the `oxdna_conf` 
        visualization function.

    Returns
    -------
    None
    """
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

def ipython_display_txt(origami_text: str, max_height: str = '500') -> None:
    """
    Render plain text (e.g., a textual representation of an origami object) as a
    scrollable HTML block in Jupyter.

    Parameters
    ----------
    origami_text : str
        The content to display in scrollable format.
    max_height : str, optional
        Maximum height of the scrollable box in pixels (default is '500').

    Returns
    -------
    None
    """
    # Convert your text to scrollable HTML
    scrollable_html = f"""
    <div style="max-height: {max_height}px; overflow-y: scroll; border: 1px solid #ccc; padding: 10px;">
    <pre>{origami_text}</pre>
    </div>
    """
    display(HTML(scrollable_html))
