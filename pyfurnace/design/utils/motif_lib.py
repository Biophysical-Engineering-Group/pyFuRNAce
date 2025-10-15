from ..core import Motif, Strand


class Utils(Motif):
    """
    Utility class that extends the Motif class with optional flipping and rotation.

    Parameters
    ----------
    *args : tuple
        Positional arguments passed to the base `Motif` class.
    hflip : bool, optional
        If True, apply a horizontal flip to the motif. Default is False.
    vflip : bool, optional
        If True, apply a vertical flip to the motif. Default is False.
    rotate : int, optional
        Rotation in degrees (usually 90, 180, 270). Default is 0.
    **kwargs : dict
        Additional keyword arguments passed to the `Motif` base class.
    """

    def __init__(
        self, *args, hflip: bool = False, vflip: bool = False, rotate: int = 0, **kwargs
    ):
        kwargs.setdefault("lock_coords", False)
        super().__init__(*args, **kwargs)
        if hflip or vflip:
            self.flip(horizontally=hflip, vertically=vflip)
        if rotate:
            self.rotate(rotate)


def start_end_stem(
    up_left: str = "3",
    up_right: str = "5",
    down_left: str = "-",
    down_right: str = "-",
    **kwargs,
) -> Utils:
    """
    Creates a `Utils` motif representing the start or end of a stem with optional
    strand labels. For each position, the acceptable values are:
    '3', '5', '─', '-', '', or None.

    Parameters
    ----------
    up_left : str or None, optional. Default is '3'.
        Label for the top-left strand.
    up_right : str or None, optional
        Label for the top-right strand. Default is '5'.
    down_left : str or None, optional
        Label for the bottom-left strand. Default is '-'.
    down_right : str or None, optional
        Label for the bottom-right strand. Default is '-'.
    **kwargs : dict
        Additional keyword arguments passed to the `Utils` constructor.

    Returns
    -------
    Utils
        An instance of the `Utils` class with the appropriate strands.
    """
    accepted_values = ["3", "5", "─", "-", "", None]

    def _check_input(value):
        if value not in accepted_values:
            raise ValueError(
                f"Invalid value for input: {value}. "
                "The value must be '3', '5', '─', '-' or None."
            )

    for val in [up_left, up_right, down_left, down_right]:
        _check_input(val)

    # Normalize None values
    up_left = up_left or ""
    up_right = up_right or ""
    down_left = down_left or ""
    down_right = down_right or ""

    if down_left and down_right and down_left in "─-" and down_right in "─-":
        down_right += "─"
    if up_left and up_right and up_left in "─-" and up_right in "─-":
        up_left += "─"

    strands = kwargs.pop("strands", [])
    if not strands:
        if up_left:
            strands.append(Strand("-" + up_left))
        if up_right:
            strands.append(Strand(up_right + "-", start=(3, 0)))
        if down_left:
            strands.append(Strand(down_left + "-", start=(1, 2), direction=(-1, 0)))
        if down_right:
            strands.append(Strand("-" + down_right, start=(4, 2), direction=(-1, 0)))

    return Utils(strands=strands, **kwargs)


def single_strand(sequence: str = "", loop: bool = False, **kwargs) -> Utils:
    """
    Creates a single-stranded region motif, optionally forming a loop.

    Parameters
    ----------
    sequence : str, optional
        The nucleotide sequence for the single-stranded region.
        Default is an empty string.
    loop : bool, optional
        If True, the single-stranded region forms a loop.
        Default is False.
    **kwargs : dict
        Additional keyword arguments passed to the `Utils` constructor.

    Returns
    -------
    Utils
        A single-stranded region `Utils` object.
    """
    bottom = "─" * len(sequence)
    strands = kwargs.pop("strands", [])
    if loop:
        strands.append(Strand(sequence + "╮│╯" + bottom))
    else:
        strands.append(Strand(sequence))
        strands.append(Strand(bottom, start=(len(sequence) - 1, 2), direction=(-1, 0)))
    return Utils(strands=strands, **kwargs)


def vertical_link(*args, **kwargs) -> Utils:
    """
    Creates a vertical link motif represented by a single vertical strand.

    Returns
    -------
    Utils
        A vertical link `Utils` object.
    """
    kwargs["strands"] = [Strand("│", direction=(0, -1))]
    return Utils(*args, **kwargs)


def vertical_double_link(*args, **kwargs) -> Utils:
    """
    Creates a double vertical link motif using two parallel vertical strands.

    Returns
    -------
    Utils
        A vertical double link `Utils` object.
    """
    kwargs["strands"] = [
        Strand("│", direction=(0, -1)),
        Strand("│", direction=(0, 1), start=(1, 0)),
    ]
    return Utils(*args, **kwargs)


def stem_cap_link(*args, **kwargs) -> Utils:
    """
    Creates a stem cap motif with a curved connection and vertical segments.

    Returns
    -------
    Utils
        A stem cap `Utils` object.
    """
    kwargs["strands"] = (
        Strand("││╭─", start=(0, 2), direction=(0, -1)),
        Strand("╭", start=(1, 2), direction=(0, -1)),
    )
    return Utils(*args, **kwargs)


# def stem_cap(*args,**kwargs):
#     kwargs['strands'] = Strand('─╰│╭─', start=(1, 2), direction=(-1, 0))
#     return Utils(*args, **kwargs)
