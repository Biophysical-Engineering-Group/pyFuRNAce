from ..core import Motif, Strand

class Utils(Motif):
    def __init__(self, *args, hflip=False, vflip=False, rotate=0, **kwargs):
        kwargs.setdefault('lock_coords', False)
        super().__init__(*args, **kwargs)
        if hflip or vflip:
            self.flip(horizontally=hflip, vertically=vflip)
        if rotate:
            self.rotate(rotate)

def start_end_stem(top_left='3', top_right='5', bot_left='-', bot_right='-', **kwargs):
    accepted_values = ['3', '5', '─', '-', '', None]
    
    def _check_input(value):
        if value not in accepted_values:
            raise ValueError(f"Invalid value for input: {value}. The value must be '3', '5', '─', '-' or None.")
        
    _check_input(top_left)
    _check_input(top_right)
    _check_input(bot_left)
    _check_input(bot_right)
    if top_left is None: top_left = ''
    if top_right is None: top_right = ''
    if bot_left is None: bot_left = ''
    if bot_right is None: bot_right = ''

    if bot_left and bot_right and bot_left in '─-' and bot_right in '─-':
        bot_right += '─'
    if top_left and top_right and top_left in '─-' and top_right in '─-':
        top_left += '─'

    if 'strands' in kwargs:
        strands = kwargs.pop('strands')
    else:
        strands = []
        if top_left:
            strands.append(Strand('-' + top_left))
        if top_right:
            strands.append(Strand(top_right + '-', start=(3, 0)))
        if bot_left:
            strands.append(Strand(bot_left + '-', start=(1, 2), direction=(-1, 0)))
        if bot_right:
            strands.append(Strand('-' + bot_right, start=(4, 2), direction=(-1, 0)))

    return Utils(strands=strands, **kwargs)

def vertical_link(*args,**kwargs):
    kwargs['strands'] = [Strand('│', direction=(0, -1))]
    return Utils(*args, **kwargs)

def vertical_double_link(*args,**kwargs):
    kwargs['strands'] = [Strand('│', direction=(0, -1)), Strand('│', direction=(0, -1), start= (1, 0))]
    return Utils(*args, **kwargs)

# def stem_cap(*args,**kwargs):
#     kwargs['strands'] = Strand('─╰│╭─', start=(1, 2), direction=(-1, 0))
#     return Utils(*args, **kwargs)

def stem_cap_link(*args,**kwargs):
    kwargs['strands'] = Strand('││╭─', start=(0, 2), direction=(0, -1)), Strand('╭', start=(1, 2), direction=(0, -1))
    return Utils(*args, **kwargs)