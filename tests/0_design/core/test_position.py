from contextlib import contextmanager
import numpy as np
import pytest

from pyfurnace.design.core.position import Position, Direction


@contextmanager
def temp_dimension(dim: int):
    """Temporarily switch global dimension; always restore to 2D at the end of
    the context."""
    try:
        Position.set_dimension(dim)
        yield
    finally:
        Position.set_dimension(2)
        # Sanity: Direction must not keep IN/OUT in 2D
        assert not hasattr(Direction, "IN")
        assert not hasattr(Direction, "OUT")


def test_position_repr_and_properties_2d_and_z_error():
    with temp_dimension(2):
        p = Position((3, -4))
        assert repr(p) == "Position(3, -4)"
        assert p.x == 3
        assert p.y == -4
        with pytest.raises(AttributeError):
            _ = p.z  # z not available in 2D


def test_position_new_pads_to_3d_and_z_property_works():
    with temp_dimension(3):
        p = Position((5, 7))  # pads z with 0
        assert p == Position((5, 7, 0))
        assert p.z == 0
        # Also accept full 3-tuple unchanged
        q = Position((1, 2, 3))
        assert q == Position((1, 2, 3))


@pytest.mark.parametrize(
    "dim, a, b, expected",
    [
        (2, Position((2, 3)), 4, Position((8, 12))),  # scalar int
        (2, Position((2, 3)), np.int64(2), Position((4, 6))),  # numpy int64
        (2, Position((2, 3)), (10, 5), Position((20, 15))),  # elementwise 2D
        (3, Position((1, 2, 3)), 3, Position((3, 6, 9))),  # scalar 3D
        (3, Position((1, 2, 3)), (4, 5, 6), Position((4, 10, 18))),  # elementwise 3D
    ],
)
def test_mul_all_branches(dim, a, b, expected):
    with temp_dimension(dim):
        # Recreate 'a' in the current dimension so it's of the correct class config
        a_cur = Position(tuple(a))
        out = a_cur * b
        assert out == expected
        assert isinstance(out, Position)


@pytest.mark.parametrize(
    "dim, a, b, expected_add, expected_sub",
    [
        (2, Position((1, 2)), Position((3, 4)), Position((4, 6)), Position((-2, -2))),
        (
            3,
            Position((1, 2, 3)),
            Position((3, 4, 5)),
            Position((4, 6, 8)),
            Position((-2, -2, -2)),
        ),
    ],
)
def test_add_sub(dim, a, b, expected_add, expected_sub):
    with temp_dimension(dim):
        a_cur, b_cur = Position(tuple(a)), Position(tuple(b))
        assert a_cur + b_cur == expected_add
        assert a_cur - b_cur == expected_sub


@pytest.mark.parametrize(
    "dim, a, expected",
    [
        (2, Position((1, -2)), Position((-1, 2))),
        (3, Position((1, -2, 5)), Position((-1, 2, -5))),
    ],
)
def test_neg(dim, a, expected):
    with temp_dimension(dim):
        a_cur = Position(tuple(a))
        assert -a_cur == expected


@pytest.mark.parametrize("dim", [2, 3])
def test_replace_and_helpers(dim):
    with temp_dimension(dim):
        base = Position((10, -20, 7)) if dim == 3 else Position((10, -20))
        # replace x only
        r1 = base.replace(x=1)
        expect = (1, base[1], base[2]) if dim == 3 else (1, base[1])
        assert r1 == Position(expect)
        # replace y only
        r2 = base.replace(y=2)
        expect = (base[0], 2, base[2]) if dim == 3 else (base[0], 2)
        assert r2 == Position(expect)
        # replace z only (3D branch)
        if dim == 3:
            r3 = base.replace(z=99)
            assert r3 == Position((base[0], base[1], 99))
        # swap_xy works in both dimensions and keeps z intact
        s = base.swap_xy()
        if dim == 3:
            assert s == Position((base[1], base[0], base[2]))
        else:
            assert s == Position((base[1], base[0]))
        # change_sign_xy (explicit branch)
        cs = base.change_sign_xy()
        if dim == 3:
            assert cs == Position((-base[0], -base[1], base[2]))
        else:
            assert cs == Position((-base[0], -base[1]))
        # swap_change_sign_xy
        scs = base.swap_change_sign_xy()
        if dim == 3:
            assert scs == Position((-base[1], -base[0], base[2]))
        else:
            assert scs == Position((-base[1], -base[0]))


def test_zero_returns_dimension_aware_vector():
    with temp_dimension(2):
        assert Position.zero() == Position((0, 0))
    with temp_dimension(3):
        assert Position.zero() == Position((0, 0, 0))


def test_set_dimension_valid_and_invalid_and_direction_changes():
    # Invalid dimension raises
    with pytest.raises(ValueError):
        Position.set_dimension(4)  # not 2 or 3

    # Switch to 3D: IN/OUT appear, and direction basis changes to 3 components
    with temp_dimension(3):
        assert hasattr(Direction, "IN")
        assert hasattr(Direction, "OUT")
        assert Direction.RIGHT == Position((1, 0, 0))
        assert Direction.DOWN == Position((0, 1, 0))
        assert Direction.LEFT == Position((-1, 0, 0))
        assert Direction.UP == Position((0, -1, 0))
        assert Direction.IN == Position((0, 0, 1))
        assert Direction.OUT == Position((0, 0, -1))

    with temp_dimension(2):
        assert hasattr(Direction, "RIGHT")
        assert Direction.RIGHT == Position((1, 0))
        assert not hasattr(Direction, "IN")
        assert not hasattr(Direction, "OUT")


def test_direction_meta_getitem_iter_names_and_immutability():
    with temp_dimension(2):
        # __getitem__ on the metaclass (Direction["RIGHT"])
        assert Direction["RIGHT"] is Direction.RIGHT
        assert Direction["LEFT"] == Position((-1, 0))

        # __iter__ yields values of non-dunder attributes; should be Positions
        vals = list(iter(Direction))
        assert set(vals) == {
            Direction.RIGHT,
            Direction.DOWN,
            Direction.LEFT,
            Direction.UP,
        }
        for v in vals:
            assert isinstance(v, Position)

        # names() returns the attribute names (non-dunder)
        names = set(Direction.names())
        assert names == {"RIGHT", "DOWN", "LEFT", "UP"}

        # Immutability via metaclass: setting or deleting should raise
        with pytest.raises(AttributeError):
            Direction.NEW = Position((9, 9))
        with pytest.raises(AttributeError):
            del Direction.RIGHT
