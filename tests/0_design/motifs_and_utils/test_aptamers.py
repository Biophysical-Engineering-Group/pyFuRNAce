from pyfurnace.design import aptamers_list, aptamers, create_aptamer
from pyfurnace.design.core import Strand
from pyfurnace.design.motifs import Aptamer, Loop


# --------------------------------------------------------------------------- #
# create_aptamer: both branches (with/without inherit_from)
# --------------------------------------------------------------------------- #


def test_create_aptamer_branches():
    # plain Aptamer
    plain = create_aptamer([Strand("─")])
    assert isinstance(plain, Aptamer) and not isinstance(plain, Loop)
    assert len(plain.strands) == 1 and plain.strands[0].strand == "─"

    # mixed-in inheritance with Loop
    mixed = create_aptamer(sequence="A", inherit_from=Loop)
    assert isinstance(mixed, Aptamer) and isinstance(mixed, Loop)
    assert len(mixed.strands) == 1


# --------------------------------------------------------------------------- #
# Aptamer library
# --------------------------------------------------------------------------- #


def test_aptamer_library_contents():
    assert isinstance(aptamers_list, list) and len(aptamers_list) > 0
    for name in aptamers_list:
        assert name in aptamers.__dict__
        apt = aptamers.__dict__[name]
        assert hasattr(apt, "__call__")
        instance = apt()
        assert isinstance(instance, Aptamer)
        for s in instance:
            assert not s._coords.is_empty()
