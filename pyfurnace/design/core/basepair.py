from collections.abc import MutableMapping
from .symbols import *
from .callback import Callback

class BasePair(MutableMapping, Callback):
    
    def __init__(self, *args, **kwargs):
        Callback.__init__(self, **kwargs)
        kwargs.pop('callback', None)
        kwargs.pop('callbacks', None)
        self._store = dict(*args, **kwargs)
        self._reverse = {v: k for k, v in self._store.items()}
        self._callbacks = []

    def __getitem__(self, key):
        if key in self._store:
            return self._store[key]
        if key in self._reverse:
            return self._reverse[key]
        raise KeyError(key)
    
    def __setitem__(self, key, value):
        if key in self._store:
            del self._reverse[self._store[key]]
        elif key in self._reverse:
            del self._store[self._reverse[key]]
        self._store[key] = value
        self._reverse[value] = key
        self._trigger_callbacks()

    def __delitem__(self, key):
        value = self._store[key]
        del self._store[key]
        del self._reverse[value]
        self._trigger_callbacks()

    def __str__(self) -> str:
        return str(self._store)

    def __repr__(self):
        return str(self._store)

    def __contains__(self, key):
        return key in self._store or key in self._reverse

    def __iter__(self):
        return iter(self._store)

    def __len__(self):
        return len(self._store)

    def keys(self):
        return self._store.keys()

    def values(self):
        return self._store.values()

    def items(self):
        return self._store.items()

    def update(self, *args, **kwargs):
        for new_dict in args:
            for k, v in new_dict.items():
                self[k] = v
        for k, v in kwargs.items():
            self[k] = v
        self._trigger_callbacks()
        return self
    
    def get(self, key, default=None):
        if key in self._store:
            return self._store[key]
        elif key in self._reverse:
            return self._reverse[key]
        return default

    def copy(self, **kwargs) -> dict:
        new_instance = BasePair(self._store, **kwargs)
        return new_instance
    
    def shift(self, shift):
        new_bp_dict = BasePair()
        for pos1, pos2 in self.items():
            new_bp_dict[(pos1[0] + shift[0], pos1[1] + shift[1])] = (pos2[0] + shift[0], pos2[1] + shift[1])
        return new_bp_dict
