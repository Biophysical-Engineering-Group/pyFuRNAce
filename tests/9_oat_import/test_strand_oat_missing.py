import pytest
import importlib
import sys
import types


def test_oat_missing(monkeypatch, tmp_path):
    """VERY IMPORTANT, THIS TEST MUST BE THE LAST ONE TO RUN IN THIS FILE
    because it messes with the import system to simulate oxDNA_analysis_tools
    not being installed. If other tests are after this one, they might fail.
    """
    # 1) Purge caches so the import block re-executes
    for name in list(sys.modules):
        if name.startswith("pyfurnace"):
            sys.modules.pop(name, None)
        if name.startswith("oxDNA_analysis_tools"):
            sys.modules.pop(name, None)
    importlib.invalidate_caches()

    # 2) Stub a non-package so submodule import fails
    fake_pkg = types.ModuleType("oxDNA_analysis_tools")
    monkeypatch.setitem(sys.modules, "oxDNA_analysis_tools", fake_pkg)

    # 3) Import the module(s) AFTER stubbing
    strand_mod = importlib.import_module("pyfurnace.design.core.strand")
    assert strand_mod.oat_installed is False

    design = importlib.import_module("pyfurnace.design")
    Strand = design.Strand
    RIGHT = design.RIGHT

    s = Strand("A\\U|G/C", directionality="53", start=(0, 0), direction=RIGHT)

    # 4) Now the call should warn because oat is "missing"
    with pytest.warns(UserWarning, match=r"oxDNA_analysis_tools is not installed"):
        s.save_3d_model(str(tmp_path), pdb=True)
