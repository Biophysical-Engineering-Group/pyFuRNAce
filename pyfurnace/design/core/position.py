from typing import Literal
import numpy as np


class Position(tuple):
    """
    A class to represent a 2D or 3D position.
    The operations are hard-coded to be 2D or 3D to
    maximize performance.

    Attributes
    ----------
    _dimension : int
        The dimension of the position (2D or 3D).
    """

    _dimension = 2

    def __new__(cls, position):
        """
        Create a new Position instance.

        Parameters
        ----------
        position : tuple
            A tuple representing the coordinates of the position.

        Returns
        -------
        Position
            A new Position instance.
        """
        if cls._dimension == 3:
            position = list(position)
            diff = 3 - len(position)
            if diff > 0:
                position.extend([0] * diff)
        return super().__new__(cls, position)

    def __repr__(self):
        """Repr of the position tuple."""
        return "Position" + super().__repr__()

    def __mul__(self, value):
        """Multiply the position by a scalar or another position element-wise."""
        if isinstance(value, (int, np.int64)):
            if self._dimension == 3:
                return Position((self[0] * value, self[1] * value, self[2] * value))
            return Position((self[0] * value, self[1] * value))
        if self._dimension == 3:
            return Position(
                (self[0] * value[0], self[1] * value[1], self[2] * value[2])
            )
        return Position((self[0] * value[0], self[1] * value[1]))

    def __add__(self, other):
        """Add the position to another position element-wise."""
        if self._dimension == 3:
            return Position(
                (self[0] + other[0], self[1] + other[1], self[2] + other[2])
            )
        return Position((self[0] + other[0], self[1] + other[1]))

    def __sub__(self, other):
        """Subtract the position to another position element-wise."""
        if self._dimension == 3:
            return Position(
                (self[0] - other[0], self[1] - other[1], self[2] - other[2])
            )
        return Position((self[0] - other[0], self[1] - other[1]))

    def __neg__(self):
        """Return the opposite coordinate."""
        if self._dimension == 3:
            return Position((-self[0], -self[1], -self[2]))
        return Position((-self[0], -self[1]))

    ###
    ### PROPERTIES
    ###

    @property
    def x(self):
        """Returns the x-coordinate."""
        return self[0]

    @property
    def y(self):
        """Returns the y-coordinate."""
        return self[1]

    @property
    def z(self):
        """Returns the z-coordinate."""
        if self._dimension == 3:
            return self[2]
        raise AttributeError("z-coordinate is not available in 2D.")

    ###
    ### CLASS METHODS
    ###

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
        Direction._set_dimension(dimension)
        cls._dimension = dimension

    @classmethod
    def zero(cls) -> "Position":
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

    ###
    ### METHODS
    ###

    def swap_xy(self):
        """Swap the x and y coordinates."""
        return self.replace(self[1], self[0])

    def change_sign_xy(self):
        """Change the sign of the x and y coordinates."""
        return self.replace(-self[0], -self[1])

    def swap_change_sign_xy(self):
        """Swap and change the sign of the x and y coordinates."""
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
            return Position(
                (
                    x if x is not None else self[0],
                    y if y is not None else self[1],
                    z if z is not None else self[2],
                )
            )
        return Position(
            (x if x is not None else self[0], y if y is not None else self[1])
        )


class DirectionMeta(type):
    """
    Metaclass to make the Direction class an Enum. In particular, it will:
        - allow iteration over the directions attribute
        - prevent setting or deleting directions
        - get the names of the directions
        - set the dimension of the directions
    """

    def __getitem__(cls, key):
        """
        Get an item from the Direction class.

        Returns
        -------
        Position
            A direction position.
        """
        return cls.__dict__[key]

    def __setattr__(cls, name, value):
        """Removes the ability to set attributes to make Direction immutable."""
        raise AttributeError("Direction is immutable")

    def __delattr__(cls, name):
        """Removes the ability to delete attributes to make Direction immutable."""
        raise AttributeError("Direction is immutable")

    def __iter__(cls):
        """
        Iterates over the Direction class attributes.

        Returns
        -------
        generator
            A generator yielding direction values.
        """
        return (value for name, value in vars(cls).items() if not name.startswith("__"))

    def _set_dimension(cls, dimension: Literal[2, 3]) -> None:
        """
        Sets the dimension of the direction class.

        Parameters
        ----------
        dimension : int
            The dimension to set (must be 2 or 3).
        """
        if dimension == 2:
            type.__setattr__(cls, "RIGHT", Position((1, 0)))
            type.__setattr__(cls, "DOWN", Position((0, 1)))
            type.__setattr__(cls, "LEFT", Position((-1, 0)))
            type.__setattr__(cls, "UP", Position((0, -1)))
            if hasattr(cls, "IN"):
                type.__delattr__(cls, "IN")
            if hasattr(cls, "OUT"):
                type.__delattr__(cls, "OUT")
        else:
            type.__setattr__(cls, "RIGHT", Position((1, 0, 0)))
            type.__setattr__(cls, "DOWN", Position((0, 1, 0)))
            type.__setattr__(cls, "LEFT", Position((-1, 0, 0)))
            type.__setattr__(cls, "UP", Position((0, -1, 0)))
            type.__setattr__(cls, "IN", Position((0, 0, 1)))
            type.__setattr__(cls, "OUT", Position((0, 0, -1)))

    def names(cls):
        """
        Returns the names of the directions.

        Returns
        -------
        list
            A list of direction names.
        """
        return [name for name in cls.__dict__ if not name.startswith("__")]


class Direction(metaclass=DirectionMeta):
    """
    Enumerator to store and retrieve the Direction information.

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

    RIGHT = Position((1, 0))
    DOWN = Position((0, 1))
    LEFT = Position((-1, 0))
    UP = Position((0, -1))
