import pytest

from pyfurnace.design.core.basepair import BasePair


def test_getitem_reverse_and_keyerror():
    bp = BasePair({1: 2, 3: 4})
    # hit reverse lookup branch (lines 54–55)
    assert bp[2] == 1
    # hit KeyError (line 56)
    with pytest.raises(KeyError):
        _ = bp["missing"]


def test_setitem_overwrite_existing_key_updates_reverse():
    # exercise line 69 True -> line 70 delete reverse[old_value]
    bp = BasePair({1: 2})
    bp[1] = 3  # overwrite existing key
    assert bp[1] == 3
    # reverse updated
    assert bp[3] == 1
    # old reverse mapping removed -> KeyError on reverse access via old value
    with pytest.raises(KeyError):
        _ = bp[2]


def test_setitem_key_exists_as_value_branch():
    # exercise line 71 True -> line 72 delete store[reverse[key]]
    bp = BasePair({1: 2})
    bp[2] = 9  # '2' exists as a value in the mapping
    assert bp[2] == 9
    assert bp[9] == 2
    # original key->value removed
    with pytest.raises(KeyError):
        _ = bp[1]


def test_delitem_removes_both_directions():
    # cover lines 86–89 (__delitem__)
    bp = BasePair({1: 2, 3: 4})
    del bp[1]
    with pytest.raises(KeyError):
        _ = bp[1]
    with pytest.raises(KeyError):
        _ = bp[2]
    # remaining pair intact
    assert bp[3] == 4
    assert bp[4] == 3


def test_str_and_repr_match_underlying_store_str():
    # cover lines 93 and 97
    bp = BasePair({1: 2})
    assert str(bp) == str(bp._store)
    assert repr(bp) == str(bp._store)
    # and they match each other
    assert str(bp) == repr(bp)


def test_iter_and_len_and_contains_cover_iter_branch():
    # cover line 105 (__iter__)
    bp = BasePair({1: 2, 5: 6})
    assert list(iter(bp)) == list(bp._store.keys())
    assert len(bp) == 2
    # __contains__ already covered but assert both directions
    assert 1 in bp
    assert 2 in bp
    assert 999 not in bp


def test_eq_with_non_mapping_triggers_type_guard():
    # hit line 113 True -> line 114 return False
    bp = BasePair({1: 2})
    assert (bp == 123) is False


def test_values_exposes_underlying_values():
    # cover line 127
    bp = BasePair({1: 2, 3: 4})
    assert list(bp.values()) == [2, 4]


def test_update_with_kwargs_branch_and_callbacks_called():
    # cover lines 152–153 (kwargs loop)
    bp = BasePair({1: 2})
    # args dict and kwargs both used
    bp.update({3: 4}, five=5)
    assert bp[1] == 2
    assert bp[3] == 4
    assert bp["five"] == 5
