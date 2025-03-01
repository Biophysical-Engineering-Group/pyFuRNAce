
class Position(tuple):
    _dimension = 2
    """Handles 2D Position behavior."""

    @classmethod
    def set_dimension(cls, dimension):
        if dimension not in [2, 3]:
            raise ValueError("Dimension must be 2 or 3.")
        if dimension not in (2, 3):
            raise ValueError("Dimension must be 2 or 3.")
        if dimension == 2:
            Direction.UP = Position((0, -1))
            Direction.RIGHT = Position((1, 0))
            Direction.DOWN = Position((0, 1))
            Direction.LEFT = Position((-1, 0))
        else:
            Direction.UP = Position((0, -1, 0))
            Direction.RIGHT = Position((1, 0, 0))
            Direction.DOWN = Position((0, 1, 0))
            Direction.LEFT = Position((-1, 0, 0))
            Direction.IN = Position((0, 0, 1))
            Direction.OUT = Position((0, 0, -1))
        cls._dimension = dimension

    @classmethod
    def zero(cls):
        if Position._dimension == 3:
            return cls((0, 0, 0))
        return cls((0, 0))

    def __repr__(self):
        return 'Position' + super().__repr__()
    
    def __mul__(self, value):
        if isinstance(value, int):
            if self._dimension == 3:
                return Position((self[0] * value, self[1] * value, self[2] * value))
            return Position((self[0] * value, self[1] * value))
        if self._dimension == 3:
            return Position((self[0] * value[0], self[1] * value[1], self[2] * value[2]))
        return Position((self[0] * value[0], self[1] * value[1]))

    def __add__(self, other):
        if isinstance(other, int):
            if self._dimension == 3:
                return Position((self[0] + other, self[1] + other, self[2] + other))
            return Position((self[0] + other, self[1] + other))
        if self._dimension == 3:
            return Position((self[0] + other[0], self[1] + other[1], self[2] + other[2]))
        return Position((self[0] + other[0], self[1] + other[1]))
    
    def __sub__(self, other):
        if isinstance(other, int):
            if self._dimension == 3:
                return Position((self[0] - other, self[1] - other, self[2] - other))
            return Position((self[0] - other, self[1] - other))
        if self._dimension == 3:
            return Position((self[0] - other[0], self[1] - other[1], self[2] - other[2]))
        return Position((self[0] - other[0], self[1] - other[1]))
    
    def swap_xy(self):
        return self.replace(self[1], self[0])
    
    def change_sign_xy(self):
        return self.replace(-self[0], -self[1])
    
    def swap_change_sign_xy(self):
        return self.replace(-self[1], -self[0])
    
    def replace(self, x=None, y=None, z=None):
        if self._dimension == 3:
            return Position((x if x is not None else self[0],
                             y if y is not None else self[1],
                             z if z is not None else self[2]))
        return Position((x if x is not None else self[0],
                         y if y is not None else self[1]))
    

class DirectionMeta(type):
    def __iter__(cls):
        return (value for name, value in vars(cls).items() if not name.startswith("__"))

class Direction(metaclass=DirectionMeta):
    """Class to store and retrieve the Direction information."""
    UP = Position((0, -1))
    RIGHT = Position((1, 0))
    DOWN = Position((0, 1))
    LEFT = Position((-1, 0))