from typing import Literal


class Position(tuple):
    """
    A class to represent a 2D or 3D position.

    Attributes
    ----------
    _dimension : int
        The dimension of the position (2D or 3D).
    """

    _dimension = 2

    @classmethod
    def set_dimension(cls, dimension: Literal[2, 3]) -> None:
        """
        Sets the dimension of the position and direction classes.

        Parameters
        ----------
        dimension : int
            The dimension to set (must be 2 or 3).

        Raises
        ------
        ValueError
            If the dimension is not 2 or 3.
        """
        if dimension not in [2, 3]:
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
    def zero(cls) -> 'Position':
        """
        Returns a zero position vector for the current dimension.

        Returns
        -------
        Position
            A zero vector of appropriate dimension.
        """
        if Position._dimension == 3:
            return cls((0, 0, 0))
        return cls((0, 0))

    def __repr__(self):
        """ Repr of the position tuple. """
        return 'Position' + super().__repr__()
    
    def __mul__(self, value):
        """ Multiply the position by a scalar or another position element-wise. """
        if isinstance(value, int):
            if self._dimension == 3:
                return Position((self[0] * value, 
                                 self[1] * value, 
                                 self[2] * value))
            return Position((self[0] * value, 
                             self[1] * value))
        if self._dimension == 3:
            return Position((self[0] * value[0], 
                             self[1] * value[1], 
                             self[2] * value[2]))
        return Position((self[0] * value[0], 
                         self[1] * value[1]))

    def __add__(self, other):
        """ Add the position with a scalar or another position element-wise. """
        if isinstance(other, int):
            if self._dimension == 3:
                return Position((self[0] + other, 
                                 self[1] + other, 
                                 self[2] + other))
            return Position((self[0] + other, 
                             self[1] + other))
        if self._dimension == 3:
            return Position((self[0] + other[0], 
                             self[1] + other[1], 
                             self[2] + other[2]))
        return Position((self[0] + other[0], 
                         self[1] + other[1]))
    
    def __sub__(self, other):
        """ Subtract the position with a scalar or another position element-wise. """
        if isinstance(other, int):
            if self._dimension == 3:
                return Position((self[0] - other, 
                                 self[1] - other, 
                                 self[2] - other))
            return Position((self[0] - other, 
                             self[1] - other))
        if self._dimension == 3:
            return Position((self[0] - other[0], 
                             self[1] - other[1], 
                             self[2] - other[2]))
        return Position((self[0] - other[0], 
                         self[1] - other[1]))
    
    def swap_xy(self):
        """ Swap the x and y coordinates. """
        return self.replace(self[1], self[0])
    
    def change_sign_xy(self):
        """ Change the sign of the x and y coordinates. """
        return self.replace(-self[0], -self[1])
    
    def swap_change_sign_xy(self):
        """ Swap and change the sign of the x and y coordinates. """
        return self.replace(-self[1], -self[0])
    
    def replace(self, x=None, y=None, z=None):
        """
        Returns a new Position with modified coordinates.
        If a coordinate is None, it keeps the current value.

        Parameters
        ----------
        x : int, optional
            New x-coordinate.
        y : int, optional
            New y-coordinate. 
        z : int, optional
            New z-coordinate (only applies in 3D).

        Returns
        -------
        Position
            A new Position with updated coordinates.
        """
        if self._dimension == 3:
            return Position((x if x is not None else self[0],
                             y if y is not None else self[1],
                             z if z is not None else self[2]))
        return Position((x if x is not None else self[0],
                         y if y is not None else self[1]))
    

class DirectionMeta(type):
    """
    Metaclass for the Direction class to allow iteration over its attributes.
    """
    def __iter__(cls):
        """
        Iterates over the Direction class attributes.

        Returns
        -------
        generator
            A generator yielding direction values.
        """
        return (value for name, value in vars(cls).items() 
                if not name.startswith("__"))
    
    def __getitem__(cls, key):
        """
        Get an item from the Direction class.

        Returns
        -------
        Position
            A direction position.
        """
        return cls.__dict__[key]

class Direction(metaclass=DirectionMeta):
    """
    Class to store and retrieve the Direction information.

    Attributes
    ----------
    UP : Position
        The upward direction.
    RIGHT : Position
        The rightward direction.
    DOWN : Position
        The downward direction.
    LEFT : Position
        The leftward direction.
    """
    UP = Position((0, -1))
    RIGHT = Position((1, 0))
    DOWN = Position((0, 1))
    LEFT = Position((-1, 0))