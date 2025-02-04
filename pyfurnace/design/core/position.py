from enum import Enum
from typing import Tuple, Union, List

class Position(tuple):
    """
    A tuple subclass representing a position in an n-dimensional space.
    
    Provides arithmetic operations such as addition, subtraction, and multiplication,
    with automatic adjustment of dimension sizes.
    """

    def __new__(cls, *args: Union[int, Tuple[int, ...], List[int]]) -> 'Position':
        """
        Create a new Position instance.

        Parameters
        ----------
        *args : Union[int, Tuple[int, ...], List[int]]
            Either individual integer arguments (e.g., Position(0, 0)) 
            or a single tuple/list containing the coordinates (e.g., Position((0, 0))).

        Returns
        -------
        Position
            A new instance of Position.

        Raises
        ------
        TypeError
            If the arguments provided are not integers, a tuple of integers, or a list of integers.
        """
        if len(args) == 1 and isinstance(args[0], (tuple, list)):
            pos_tuple = tuple(args[0])  # Handle single tuple or list argument
        else:
            # Handle individual integer arguments
            pos_tuple = tuple(args)
            
        if len(pos_tuple) <3:
            pos_tuple = pos_tuple + (0,)*(3-len(pos_tuple))

        if not all(isinstance(x, int) for x in pos_tuple):
            raise TypeError("Position arguments must be integers or a tuple/list of integers.")

        return super(Position, cls).__new__(cls, pos_tuple)
    
    def __getitem__(self, index: Union[int, slice]) -> Union[int, 'Position']:
        """
        Get the value or values at the specified index.

        Parameters
        ----------
        index : Union[int, slice]
            The index of the coordinate to retrieve. It can be an integer for a single value
            or a slice for a subset of the position.

        Returns
        -------
        Union[int, Position]
            The coordinate value at the specified index or a new Position if a slice is used.

        Raises
        ------
        IndexError
            If the index is out of range.
        TypeError
            If the index is not an integer or slice.
        """
        if isinstance(index, int):
            return super().__getitem__(index)
        elif isinstance(index, slice):
            return Position(super().__getitem__(index))
        else:
            raise TypeError("Index must be an integer or a slice.")


    def __str__(self) -> str:
        """
        Return a string representation of the Position.

        Returns
        -------
        str
            A string describing the Position with its dimensions.
        """
        dim_text = ', '.join([f"d{i}: {v}" for i, v in enumerate(self)])
        return f"Position({dim_text})"

    def __mul__(self, other: Union['Position', 'Direction', Tuple[int, ...], int]) -> 'Position':
        """
        Perform element-wise multiplication between this Position and another.

        Parameters
        ----------
        other : Union[Position, Direction, Tuple[int, ...], int]
            The other object to multiply with this Position.

        Returns
        -------
        Position
            The resulting Position after multiplication.
        """
        tuple1, tuple2 = self._adjust(self, other)
        return Position([a * b for a, b in zip(tuple1, tuple2)])

    def __add__(self, other: Union['Position', 'Direction', Tuple[int, ...]]) -> 'Position':
        """
        Perform element-wise addition between this Position and another.

        Parameters
        ----------
        other : Union[Position, Direction, Tuple[int, ...]]
            The other object to add to this Position.

        Returns
        -------
        Position
            The resulting Position after addition.
        """
        tuple1, tuple2 = self._adjust(self, other)
        return Position([a + b for a, b in zip(tuple1, tuple2)])

    def __iadd__(self, other: Union['Position', 'Direction', Tuple[int, ...]]) -> 'Position':
        """
        Perform in-place element-wise addition between this Position and another.

        Parameters
        ----------
        other : Union[Position, Direction, Tuple[int, ...]]
            The other object to add to this Position.

        Returns
        -------
        Position
            The modified Position after in-place addition.
        """
        return self.__add__(other)

    def __sub__(self, other: Union['Position', 'Direction', Tuple[int, ...]]) -> 'Position':
        """
        Perform element-wise subtraction between this Position and another.

        Parameters
        ----------
        other : Union[Position, Direction, Tuple[int, ...]]
            The other object to subtract from this Position.

        Returns
        -------
        Position
            The resulting Position after subtraction.
        """
        tuple1, tuple2 = self._adjust(self, other)
        return Position([a - b for a, b in zip(tuple1, tuple2)])

    @staticmethod
    def _adjust(iter1: Union[Tuple[int, ...], 'Position'], 
                iter2: Union[Tuple[int, ...], 'Position', 'Direction']) -> Tuple['Position', 'Position']:
        """
        Adjust two iterables to have the same length by padding with zeros.

        Parameters
        ----------
        iter1 : Union[Tuple[int, ...], Position]
            The first iterable.
        iter2 : Union[Tuple[int, ...], Position, Direction]
            The second iterable.

        Returns
        -------
        Tuple[Position, Position]
            Two Position objects of the same length.

        Raises
        ------
        TypeError
            If the provided iterables are not of the expected types.
        """
        for pos in (iter1, iter2):
            if not isinstance(pos, (tuple, list, Position, Direction)):
                raise TypeError(f"Cannot adjust Position with {type(pos)}")

        iter1, iter2 = tuple(iter1), tuple(iter2)
        max_length = max(len(iter1), len(iter2))
        iter1 = Position(iter1 + (0,) * (max_length - len(iter1)))
        iter2 = Position(iter2 + (0,) * (max_length - len(iter2)))
        return iter1, iter2


class Direction(Enum):
    """
    An enumeration representing cardinal directions in a 2D space.
    """
    UP = (0, -1)
    DOWN = (0, 1)
    LEFT = (-1, 0)
    RIGHT = (1, 0)

    def __str__(self) -> str:
        """
        Return a string representation of the Direction.

        Returns
        -------
        str
            The direction as a string showing the x and y components.
        """
        return f"Direction(x={self.value[0]}, y={self.value[1]})"

    def __repr__(self) -> str:
        """
        Return a short symbolic representation of the Direction.

        Returns
        -------
        str
            A symbolic string representing the direction.
        """
        if self == Direction.UP:
            return "^"
        if self == Direction.DOWN:
            return "v"
        if self == Direction.LEFT:
            return "<"
        if self == Direction.RIGHT:
            return ">"

    def __mul__(self, other: Union['Direction', 'Position', int]) -> 'Position':
        """
        Perform element-wise multiplication between this Direction and another object.

        Parameters
        ----------
        other : Union[Direction, Position, int]
            The object to multiply with this Direction.

        Returns
        -------
        Position
            The resulting Position after multiplication.
        """
        if isinstance(other, Direction):
            return Position((self.value[0] * other.value[0], self.value[1] * other.value[1]))
        if isinstance(other, Position):
            return Position((self.value[0] * other[0], self.value[1] * other[1]))
        if isinstance(other, int):
            return Position((self.value[0] * other, self.value[1] * other))
        return NotImplemented

    def __add__(self, other: Union['Direction', 'Position', int]) -> 'Position':
        """
        Perform element-wise addition between this Direction and another object.

        Parameters
        ----------
        other : Union[Direction, Position, int]
            The object to add to this Direction.

        Returns
        -------
        Position
            The resulting Position after addition.
        """
        if isinstance(other, Direction):
            return Position((self.value[0] + other.value[0], self.value[1] + other.value[1]))
        if isinstance(other, Position):
            return Position((self.value[0] + other[0], self.value[1] + other[1]))
        if isinstance(other, int):
            return Position((self.value[0] + other, self.value[1] + other))
        return NotImplemented

    def __sub__(self, other: Union['Direction', 'Position', int]) -> 'Position':
        """
        Perform element-wise subtraction between this Direction and another object.

        Parameters
        ----------
        other : Union[Direction, Position, int]
            The object to subtract from this Direction.

        Returns
        -------
        Position
            The resulting Position after subtraction.
        """
        if isinstance(other, Direction):
            return Position((self.value[0] - other.value[0], self.value[1] - other.value[1]))
        if isinstance(other, Position):
            return Position((self.value[0] - other[0], self.value[1] - other[1]))
        if isinstance(other, int):
            return Position((self.value[0] - other, self.value[1] - other))
        return NotImplemented

    def __hash__(self) -> int:
        """
        Return the hash value of the Direction.

        Returns
        -------
        int
            The hash value.
        """
        return hash(self.value)

    def __eq__(self, value: object) -> bool:
        """
        Check if this Direction is equal to another object.

        Parameters
        ----------
        value : object
            The object to compare with.

        Returns
        -------
        bool
            True if equal, otherwise False.
        """
        if isinstance(value, Direction):
            return self.value == value.value
        return self.value == value
