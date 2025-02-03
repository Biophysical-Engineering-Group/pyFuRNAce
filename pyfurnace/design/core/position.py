from enum import Enum

class Position(tuple):
    def __new__ (cls, pos_tuple):
        return super(Position, cls).__new__(cls, pos_tuple)

    def __str__(self):
        return f"Position(x={self[0]}, y={self[1]})"

    def __mul__(self, other):
        if isinstance(other, Direction):
            return Position((self[0] * other.value[0], self[1] * other.value[1]))
        if isinstance(other, Position):
            return Position((self[0] * other[0], self[1] * other[1]))
        if isinstance(other, int):
            return Position((self[0] * other, self[1] * other))
        return NotImplemented

    def __add__(self, other):
        if isinstance(other, Direction):
            return Position((self[0] + other.value[0], self[1] + other.value[1]))
        if isinstance(other, Position):
            return Position((self[0] + other[0], self[1] + other[1]))
        if isinstance(other, int):
            return Position((self[0] + other, self[1] + other))
        return NotImplemented

    def __sub__(self, other):
        if isinstance(other, Direction):
            return Position((self[0] - other.value[0], self[1] - other.value[1]))
        if isinstance(other, Position):
            return Position((self[0] - other[0], self[1] - other[1]))
        if isinstance(other, int):
            return Position((self[0] - other, self[1] - other))
        return NotImplemented

class Direction(Enum):
    UP = (0, -1)
    DOWN = (0, 1)
    LEFT = (-1, 0)
    RIGHT = (1, 0)

    def __str__(self):
        return f"Direction(x={self.value[0]}, y={self.value[1]})"
    
    def __repr__(self):
        if self == Direction.UP:
            return "^"
        if self == Direction.DOWN:
            return "v"
        if self == Direction.LEFT:
            return "<"
        if self == Direction.RIGHT:
            return ">"

    def __mul__(self, other):
        if isinstance(other, Direction):
            return Position((self.value[0] * other.value[0], self.value[1] * other.value[1]))
        if isinstance(other, Position):
            return Position((self.value[0] * other[0], self.value[1] * other[1]))
        if isinstance(other, int):
            return Position((self.value[0] * other, self.value[1] * other))
        return NotImplemented

    def __add__(self, other):
        if isinstance(other, Direction):
            return Position((self.value[0] + other.value[0], self.value[1] + other.value[1]))
        if isinstance(other, Position):
            return Position((self.value[0] + other[0], self.value[1] + other[1]))
        if isinstance(other, int):
            return Position((self.value[0] + other, self.value[1] + other))
        return NotImplemented

    def __sub__(self, other):
        if isinstance(other, Direction):
            return Position((self.value[0] - other.value[0], self.value[1] - other.value[1]))
        if isinstance(other, Position):
            return Position((self.value[0] - other[0], self.value[1] - other[1]))
        if isinstance(other, int):
            return Position((self.value[0] - other, self.value[1] - other))
        return NotImplemented
    
    def __hash__(self) -> int:
        return hash(self.value)
    
    def __eq__(self, value: object) -> bool:
        if isinstance(value, (Direction)):
            return self.value == value
        else:
            return self.value == value