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
    def __init__(self, 
                 *args, 
                 hflip: bool = False, 
                 vflip : bool=False, 
                 rotate : int =0, 
                 **kwargs):
        kwargs.setdefault('lock_coords', False)
        super().__init__(*args, **kwargs)
        if hflip or vflip:
            self.flip(horizontally=hflip, vertically=vflip)
        if rotate:
            self.rotate(rotate)


def start_end_stem(top_left: str = '3', 
                   top_right: str = '5', 
                   bot_left: str = '-', 
                   bot_right: str = '-', 
                   **kwargs) -> Utils:
    """
    Creates a `Utils` motif representing the start or end of a stem with optional
    strand labels. For each position, the acceptable values are:
    '3', '5', '─', '-', '', or None.

    Parameters
    ----------
    top_left : str or None, optional. Default is '3'.
        Label for the top-left strand.
    top_right : str or None, optional
        Label for the top-right strand. Default is '5'.
    bot_left : str or None, optional
        Label for the bottom-left strand. Default is '-'.
    bot_right : str or None, optional
        Label for the bottom-right strand. Default is '-'.
    **kwargs : dict
        Additional keyword arguments passed to the `Utils` constructor.

    Returns
    -------
    Utils
        An instance of the `Utils` class with the appropriate strands.
    """
    accepted_values = ['3', '5', '─', '-', '', None]

    def _check_input(value):
        if value not in accepted_values:
            raise ValueError(
                f"Invalid value for input: {value}. "
                "The value must be '3', '5', '─', '-' or None."
            )

    for val in [top_left, top_right, bot_left, bot_right]:
        _check_input(val)

    # Normalize None values
    top_left = top_left or ''
    top_right = top_right or ''
    bot_left = bot_left or ''
    bot_right = bot_right or ''

    if bot_left and bot_right and bot_left in '─-' and bot_right in '─-':
        bot_right += '─'
    if top_left and top_right and top_left in '─-' and top_right in '─-':
        top_left += '─'

    strands = kwargs.pop('strands', [])
    if not strands:
        if top_left:
            strands.append(Strand('-' + top_left))
        if top_right:
            strands.append(Strand(top_right + '-', start=(3, 0)))
        if bot_left:
            strands.append(Strand(bot_left + '-', start=(1, 2), direction=(-1, 0)))
        if bot_right:
            strands.append(Strand('-' + bot_right, start=(4, 2), direction=(-1, 0)))

    return Utils(strands=strands, **kwargs)


def vertical_link(*args, **kwargs) -> Utils:
    """
    Creates a vertical link motif represented by a single vertical strand.

    Returns
    -------
    Utils
        A vertical link `Utils` object.
    """
    kwargs['strands'] = [Strand('│', direction=(0, -1))]
    return Utils(*args, **kwargs)


def vertical_double_link(*args, **kwargs) -> Utils:
    """
    Creates a double vertical link motif using two parallel vertical strands.

    Returns
    -------
    Utils
        A vertical double link `Utils` object.
    """
    kwargs['strands'] = [
        Strand('│', direction=(0, -1)),
        Strand('│', direction=(0, -1), start=(1, 0))
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
    kwargs['strands'] = (
        Strand('││╭─', start=(0, 2), direction=(0, -1)),
        Strand('╭', start=(1, 2), direction=(0, -1))
    )
    return Utils(*args, **kwargs)

# def stem_cap(*args,**kwargs):
#     kwargs['strands'] = Strand('─╰│╭─', start=(1, 2), direction=(-1, 0))
#     return Utils(*args, **kwargs)