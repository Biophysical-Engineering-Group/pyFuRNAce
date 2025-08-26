from IPython.display import display, HTML
from matplotlib import colors as mcolors, pyplot as plt
import numpy as np
from typing import List, Optional, Union, Any
import tempfile

# try to import the oxDNA_analysis_tools package
try:
    from oxDNA_analysis_tools.UTILS.oxview import oxdna_conf
    from oxDNA_analysis_tools.UTILS.RyeReader import describe, get_confs

    oat_installed = True
except ModuleNotFoundError:
    oat_installed = False

# own imports
from ..core import *
from ..motifs import *
from .motif_lib import *

# A dictionary to convert angles to dovetail values
ANGLES_DT_DICT = {
    26: -6,
    58: -5,
    90: -4,
    122: -3,
    154: -2,
    186: -1,
    218: 0,
    250: 1,
    282: 2,
    314: 3,
    346: 4,
    378: 5,
    410: 6,
}


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
    dt_list = [
        ANGLES_DT_DICT[min(ANGLES_DT_DICT, key=lambda x: abs(x - ang))]
        for ang in angles_sanitize
    ]
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
    align: str = "first",
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

    # convert angles to dovetail values if needed
    if use_angles:
        dt_list = convert_angles_to_dt(dt_list)

    # add the start and end helix to the dovetail list
    if add_terminal_helix:
        dt_list = [0] + dt_list + [0]

    ### ADJUST THE MAIN STEM MATRIX

    if main_stem is None:  # set it to the minimum value for each KL
        max_dt = max([abs(dt) for dt in dt_list], default=0)
        main_stem = [[11 * ((max_dt + 17) // 11 + 1)] * kl_columns] * len(dt_list)

    elif isinstance(main_stem, int):
        main_stem = [
            [main_stem for _ in range(kl_columns)] for _ in range(len(dt_list))
        ]

    elif isinstance(main_stem, list) and all(isinstance(x, int) for x in main_stem):
        main_stem = [main_stem for _ in range(len(dt_list))]

    elif isinstance(main_stem, list) and all(
        isinstance(x, (tuple, list)) for x in main_stem
    ):

        if not all(len(x) == kl_columns for x in main_stem):
            raise ValueError(
                "The main_stem list should have the same length"
                " as the kissing loops repeats"
            )
    else:
        raise ValueError(
            "The main_stem can be an int, a list of int or a" " matrix of int"
        )

    ### ADJUST THE LEFT KL STEM MATRIX

    if left_stem_kl is None:
        left_stem_kl = [[None] * kl_columns for _ in range(len(dt_list))]

    elif isinstance(left_stem_kl, int):
        left_stem_kl = [
            [left_stem_kl for _ in range(kl_columns)] for _ in range(len(dt_list))
        ]

    elif isinstance(left_stem_kl, list) and all(
        isinstance(x, int) for x in left_stem_kl
    ):
        left_stem_kl = [[left_stem_kl[i]] * kl_columns for i in range(len(dt_list))]

    elif isinstance(left_stem_kl, list) and all(
        isinstance(x, (tuple, list)) for x in left_stem_kl
    ):

        if not all(len(x) == kl_columns for x in left_stem_kl):
            raise ValueError(
                "The left_stem_kl list should have the same length "
                "as the kissing loops repeats"
            )
    else:
        raise ValueError(
            "The left_stem_kl can be an int, a list of int or a " "matrix of int"
        )

    if stem_pos is None:
        stem_pos = [0 for _ in range(kl_columns)]
    elif isinstance(stem_pos, int):
        stem_pos = [stem_pos for _ in range(kl_columns)]

    ### BUILD THE ORIGAMI STRUCTURE, helix by helix
    for helix_in, dt in enumerate(dt_list):

        # create the start of the stem
        helix = [TetraLoop(), Stem(end_helix_len), Dovetail(dt)]

        # add Kissing loops repeats to the helix
        for kl_index in range(kl_columns):

            # calculate the stem lengths
            stem_len = main_stem[helix_in][kl_index]
            left_stem = left_stem_kl[helix_in][kl_index]
            if left_stem is None:
                left_stem = (stem_len - 8 - abs(dt)) // 2
            right_stem = (stem_len - 8 - abs(dt)) - left_stem

            # this is a position where to add a continuous stem
            if stem_pos[kl_index] == helix_in:

                # add the start position in this stem
                if kl_index == start and add_start_end:
                    half_l_stem = (stem_len - abs(dt)) // 2
                    half_r_stem = stem_len - abs(dt) - half_l_stem
                    helix += [
                        Stem(half_l_stem).shift((1, 0), extend=True),
                        start_end_stem(),
                        Stem(half_r_stem),
                        Dovetail(dt),
                    ]
                else:
                    stem_len = main_stem[helix_in][kl_index] - abs(dt)
                    helix += [Stem(stem_len).shift((6, 0), extend=True), Dovetail(dt)]

            # normal kissing loop repeat
            else:
                helix += [
                    Stem(left_stem),
                    KissingDimer(),
                    Stem(right_stem),
                    Dovetail(dt),
                ]

        # add the end of the helix
        helix += [Stem(end_helix_len), TetraLoop(open_left=True)]

        # add the helix to the origami
        origami.append(helix, copy=False)

    # remove the top cross from the dovetails of the first helix
    for motif in origami[0]:
        if isinstance(motif, Dovetail):
            motif.up_cross = False

    # remove the bottom cross from the dovetails of the last helix
    for motif in origami[-1]:
        if isinstance(motif, Dovetail):
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
        warnings.warn(
            "The oxDNA_analysis_tools package is not installed, "
            "the 3D display is not available."
        )
        return
    # Create a temporary directory
    with tempfile.TemporaryDirectory() as tmpdirname:
        file_path = f"{tmpdirname}/origami"
        origami.save_3d_model(file_path)
        top_info, traj_info = describe(f"{file_path}.top", f"{file_path}.dat")
        conf = get_confs(top_info, traj_info, 0, 1)[0]
        oxdna_conf(top_info, conf, **kwargs)


def ipython_display_txt(origami_text: str, max_height: str = "500") -> None:
    """
    Render plain text (e.g., a textual representation of an origami object) as a
    scrollable HTML block in Jupyter.

    Parameters
    ----------
    origami_text : str, Origami
        The content to display in scrollable format.
    max_height : str, optional
        Maximum height of the scrollable box in pixels (default is '500').

    Returns
    -------
    None
    """
    # Convert your text to scrollable HTML
    ori_txt = str(origami_text).replace("\n", "<br>")
    scrollable_html = (
        f'<div style="max-height: {max_height}px; white-space: nowrap;'
        f"overflow-x: auto; overflow-y: scroll;"
        f"font-family: monospace;"
        f'border: 1px solid #ccc; padding: 10px;">'
        f"{ori_txt.replace(' ', '&nbsp;')}"
        "</div>"
    )
    display(HTML(scrollable_html))


def ipython_clickable_txt(
    origami: Origami,
    max_height: Union[str, int] = "500",
    barriers: Optional[Any] = None,
    gradient: Union[bool, str] = False,
    font_size: int = 12,
) -> str:
    """
    Generate an interactive, scrollable HTML view of a RNA origami structure
    with clickable motifs that display their indexes in a JavaScript alert.

    Parameters
    ----------
    origami : Origami
        An Origami object representing the RNA structure.
    max_height : str or int, optional
        The maximum height of the scrollable view in pixels. Default is '500'.
    barriers : optional
        Optional barrier data to overlay on the origami representation.
    gradient : bool or str, optional
        Whether to color motifs using a gradient.
        If str, interpreted as colormap name from matplotlib.
        Default is False.
    font_size : int, optional
        Font size of the text in pixels. Default is 12.

    Returns
    -------
    str
        An HTML string rendered via IPython's display system.

    Notes
    -----
    - This function uses inline CSS and JavaScript for visual styling and interaction.
    - Clicking a motif will trigger a JavaScript alert showing its position.
    """
    barriers_colors = {"▂": "#FFBA08", "▄": "#FFBA08", "█": "#D00000"}
    high_color = "#D00000"
    normal_color = "inherit"

    motif = origami.assembled
    # create a dictionary from positions to index
    pos_to_index = {pos: ind for ind, pos in enumerate(motif.seq_positions)}

    if barriers:
        origami_lines = origami.barrier_repr(return_list=True)
    else:
        origami_str = str(origami)
        origami_lines = origami_str.split("\n")

    # create color gradient
    if gradient:
        tot_len = 0
        for s in origami.strands:
            tot_len += len(s.sequence)
            for protein in s.coords.proteins:
                tot_len += len(protein)
        cmap = plt.get_cmap(gradient)
        c_map = [mcolors.to_hex(cmap(i)) for i in np.linspace(0, 1, tot_len)]

    # Prepare the string to add 5' and 3' symbols for the strands
    motif_list = (
        [[" "] * (motif.num_char + 2)]
        + [[" "] + [char for char in line] + [" "] for line in origami_lines]
        + [[" "] * (motif.num_char + 2)]
    )

    for s in motif:  # Add the 5' and 3' symbols to the motif as 1 and 2
        if not s.sequence:
            continue
        if (
            s.sequence
            and s[0] not in "35"
            and motif_list[s.prec_pos[1] + 1][s.prec_pos[0] + 1] == " "
        ):
            if s.directionality == "53":
                motif_list[s.prec_pos[1] + 1][s.prec_pos[0] + 1] = "1"
            else:
                motif_list[s.prec_pos[1] + 1][s.prec_pos[0] + 1] = "2"
        if (
            s.sequence
            and s[-1] not in "35"
            and motif_list[s.next_pos[1] + 1][s.next_pos[0] + 1] == " "
        ):
            if s.directionality == "53":
                motif_list[s.next_pos[1] + 1][s.next_pos[0] + 1] = "2"
            else:
                motif_list[s.next_pos[1] + 1][s.next_pos[0] + 1] = "1"

    origami_list = ["".join(line) for line in motif_list]

    content = (
        f"<div style='font-family: monospace;"
        f"font-size: {font_size}px; "
        "white-space: nowrap; "
        "overflow-x: auto; "
        "overflow-y: scroll; "
        f"max-height: {max_height}px;'>"
    )
    span = '<span style="font-family: monospace; '

    for y, line in enumerate(origami_list):

        for x, char in enumerate(line):
            ori_pos = (x - 1, y - 1)
            color = normal_color
            if barriers_colors and char in barriers_colors:
                color = barriers_colors[char]

            if char == " ":
                content += span + 'line-height:1;">&nbsp;</span>'
            elif char == "1":
                content += span + f'color: {high_color}; line-height:1;">5</span>'
            elif char == "2":
                content += span + f'color: {high_color}; line-height:1;">3</span>'

            elif char in bp_symbols:  # do not highlight the base pair in red
                content += span + f'color: {color}; line-height:1;">{char}</span>'

            elif ori_pos in origami.pos_index_map:  # a motif symbol
                sl = origami.pos_index_map[ori_pos]
                index = pos_to_index.get(ori_pos)
                msg_text = f"Line {sl[0]}, Motif {sl[1]}"
                if index is not None:
                    msg_text = f"Base {index}, Line {sl[0]}, Motif {sl[1]}"
                    if gradient:
                        color = c_map[index]

                content += (
                    f'<a style="text-decoration: none;'
                    "font-family: monospace; "
                    f"color: {color}; "
                    'line-height:1;" '
                    f'href="#/" '
                    f"""onclick="alert('{msg_text}')" """
                    f'id="{sl[0]},{sl[1]},{x - 1},{y - 1}">'
                    f"{char}"
                    "</a>"
                )

            else:  # is a junction symbol
                content += span + f'color: {color}; line-height:1;">{char}</span>'
        content += "<br />"
    content += "</div>"

    return display(HTML(content))


def template_2_helix():
    """
    Generate the RNA origami template for a 2-helix structure.
    Reference:
    1. Krissanaprasit, A. et al. A functional RNA-origami as direct thrombin
        inhibitor with fast-acting and specific single-molecule reversal agents in
        vivo model. Molecular Therapy 32, 2286-2298 (2024).

    Returns
    -------
    Origami
        An Origami object representing the 2-helix RNA structure.
    """
    import pyfurnace as pf

    line_0 = [
        pf.TetraLoop(),
        pf.Stem(11),
        pf.Dovetail(0, up_cross=False),
        pf.Stem(11),
        pf.start_end_stem(),
        pf.Stem(11),
        pf.Dovetail(0, up_cross=False),
        pf.Stem(11),
        pf.TetraLoop(open_left=True),
    ]
    line_1 = [
        pf.TetraLoop(),
        pf.Stem(11),
        pf.Dovetail(0, down_cross=False),
        pf.Stem(7),
        pf.KissingDimer(),
        pf.Stem(7),
        pf.Dovetail(0, down_cross=False),
        pf.Stem(11),
        pf.TetraLoop(open_left=True),
    ]
    origami = pf.Origami([line_0, line_1], align="first")
    return origami


def template_rna_filament():
    """
    Generate the RNA origami template for an RNA filament structure.

    Reference:
    1. Tran, M. P. et al. Genetic encoding and expression of RNA origami
        cytoskeletons in synthetic cells. Nat. Nanotechnol. 20, 664-671 (2025).

    Returns
    -------
    Origami
        An Origami object representing the RNA filament structure.
    """
    import pyfurnace as pf

    origami = pf.simple_origami(
        dt_list=[-3],
        kl_columns=1,
        main_stem=[22],
        align="first",
    )

    # adjust middle helix terminal stem
    origami[(1, 7)].length = 5
    origami[(1, 1)].length = 5

    ### STEMS BEFORE THE EXTERNAL KL
    origami[(0, 7)].length = 6
    origami[(2, 7)].length = 6
    origami[(2, 1)].length = 8
    origami[(0, 1)].length = 7

    origami[(0, 0)] = pf.KissingLoop180(open_left=False, pk_index="2")
    origami[(-1, 0)] = pf.KissingLoop180(open_left=False, pk_index="1")
    origami[(0, -1)] = pf.KissingLoop180(open_left=True, pk_index="1'")
    origami[(-1, -1)] = pf.KissingLoop180(open_left=True, pk_index="2'")
    return origami


def template_rna_filament_ispinach():
    """
    Generate the RNA origami template for an RNA filament structure.
    The tiles contains an ispinach aptamer for fluorescence imaging.

    Reference:
    1. Tran, M. P. et al. Genetic encoding and expression of RNA origami
        cytoskeletons in synthetic cells. Nat. Nanotechnol. 20, 664-671 (2025).

    Returns
    -------
    Origami
        An Origami object representing the RNA filament structure.
    """
    import pyfurnace as pf

    origami = pf.simple_origami(
        dt_list=[-3],
        kl_columns=1,
        main_stem=[22],
        align="first",
    )

    # adjust middle helix terminal stem
    origami[(1, 7)].length = 5
    origami[(1, 1)].length = 5

    ### STEMS BEFORE THE EXTERNAL KL
    origami[(0, 7)].length = 6
    origami[(2, 7)].length = 6
    origami[(2, 1)].length = 8
    origami[(0, 1)].length = 7

    origami[(0, 0)] = pf.KissingLoop180(open_left=False, pk_index="2")
    origami[(-1, 0)] = pf.KissingLoop180(open_left=False, pk_index="1")
    origami[(0, -1)] = pf.KissingLoop180(open_left=True, pk_index="1'")
    origami[(-1, -1)] = pf.KissingLoop180(open_left=True, pk_index="2'")

    # add connector for ispinach line
    origami.insert((0, 4), pf.Dovetail(0, down_cross=False))

    # add ispinach line
    ispi_line = [
        pf.TetraLoop(),
        pf.Stem(3, strong_bases=False),
        pf.Ispinach().flip(),
        pf.Stem(3),
        pf.Motif.from_structure("...&", "UUU&").flip(),
        pf.stem_cap_link(hflip=True),
    ]
    origami.insert(0, ispi_line)
    return origami


def template_3_arms_droplet():
    """
    Generate the RNA origami template for a 3-arms droplet structure.

    References:
    1. Stewart, J. M. et al. Modular RNA motifs for orthogonal phase separated
        compartments. Nat Commun 15, (2024).
    2. Fabrini, G. et al. Co-transcriptional production of programmable RNA condensates
        and synthetic organelles. Nat. Nanotechnol. 19, 1665-1673 (2024).
    3. Monari, L., Braun, I., Poppleton, E. & Göpfrich, K. PyFuRNAce: An integrated
        design engine for RNA origami. (2025) doi:10.1101/2025.04.17.647389.

    Returns
    -------
    Origami
        An Origami object representing the 3-arms droplet structure.
    """
    import pyfurnace as pf

    palindr_kl = pf.KissingLoop(sequence="AUCGCGAAA")
    line_0 = [
        palindr_kl.copy(),
        pf.Stem(8),
        pf.start_end_stem(),
        pf.Stem(17),
        pf.Motif.from_structure(".&", "U&").flip(),
        pf.Dovetail(0, up_cross=False),
        pf.Motif.from_structure(".&.", "U&U"),
        pf.Stem(25),
        palindr_kl.copy().flip(),
    ]
    line_1 = [
        palindr_kl.copy(),
        pf.Stem(5),
        pf.Broccoli(),
        pf.Stem(8),
        pf.stem_cap_link().flip(),
    ]
    origami = pf.Origami([line_0, line_1], align="first")
    return origami


def template_4_arms_droplet():
    """
    Generate the RNA origami template for a 4-arms droplet structure.

    Reference:
    1. Fabrini, G. et al. Co-transcriptional production of programmable RNA condensates
        and synthetic organelles. Nat. Nanotechnol. 19, 1665-1673 (2024).

    Returns
    -------
    Origami
        An Origami object representing the 4-arms droplet structure.
    """
    import pyfurnace as pf

    uracil_connect = pf.Motif.from_structure(".&.", "U&U")
    palindr_kl = pf.KissingLoop(sequence="AUCGCGAAA")

    line_0 = [
        palindr_kl.copy(),
        pf.Stem(10),
        pf.MalachiteGreenShort(),
        pf.Stem(10),
        uracil_connect.copy(),
        pf.Dovetail(0, up_cross=False),
        pf.Stem(25),
        palindr_kl.copy().flip(),
    ]

    line_1 = [
        palindr_kl.copy(),
        pf.Stem(13),
        pf.start_end_stem(),
        pf.Stem(12),
        pf.Dovetail(0, down_cross=False),
        uracil_connect.copy(),
        pf.Stem(25),
        palindr_kl.copy().flip(),
    ]
    origami = pf.Origami([line_0, line_1], align="first")
    return origami


def template_pentagon_tile():
    """
    Generate the RNA origami template for a pentagon tile structure.
    This follows the 3H-4DT design from Geary et al. (2021).

    Reference:
    1. Geary, C., Grossi, G., McRae, E. K. S., Rothemund, P. W. K. & Andersen, E. S.
        RNA origami design tools enable cotranscriptional folding of kilobase-sized
        nanoscaffolds. Nat. Chem. 13, 549-558 (2021).

    Returns
    -------
    Origami
        An Origami object representing the pentagon tile structure.
    """
    import pyfurnace as pf

    origami = pf.simple_origami([-4], main_stem=33, end_helix_len=7)

    # adjust start position
    origami[0, 3].length = 12
    origami[0, 5].length = 21

    # adjust stems before the KL
    origami[0, 1].length = 5
    origami[-1, 1].length = 17
    origami[0, -2].length = 12
    origami[-1, -2].length = 21

    # put the KL in the right place
    origami[0, 0] = pf.KissingLoop120(pk_index=1)
    origami[0, -1] = pf.KissingLoop120(open_left=True, pk_index=-1)
    origami[-1, 0] = pf.KissingLoop120(pk_index=2)
    origami[-1, -1] = pf.KissingLoop120(open_left=True, pk_index=-2)

    return origami


def template_rectangle_10H_3X():
    """
    Generate the RNA origami template for a rectangle 10H-3X structure.

    References:
    1. Monari, L., Braun, I., Poppleton, E. & Göpfrich, K. PyFuRNAce: An integrated
        design engine for RNA origami. (2025) doi:10.1101/2025.04.17.647389.

    Returns
    -------
    Origami
        An Origami object representing the rectangle 10H-3X structure.
    """
    import pyfurnace as pf

    origami = pf.simple_origami(
        dt_list=[180] * 8,
        kl_columns=3,
        main_stem=33,
        add_terminal_helix=True,
        align="first",
        use_angles=True,
    )  # Create a simple origami

    origami = origami.improve_folding_pathway(kl_delay=150)

    origami.insert((0, 11), pf.Broccoli().flip(1, 1))  # Add motif

    origami.insert((0, 12), pf.Stem(5))  # Add motif

    return origami
