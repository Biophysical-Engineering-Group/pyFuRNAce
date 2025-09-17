from collections.abc import MutableMapping
from typing import Any, Dict, Iterator

# pyFuRNAce IMPORTS
from .callback import Callback


class BasePair(MutableMapping, Callback):
    """
    A bidirectional dictionary that maintains a mapping between keys and values,
    ensuring that each value has a unique corresponding key and vice versa.

    Inherits from MutableMapping to provide dictionary-like behavior and from
    Callback to support event-driven notifications on updates.

    Parameters
    ----------
    *args : tuple
        Positional arguments passed to the internal dictionary.
    **kwargs : dict
        Keyword arguments passed to the internal dictionary. Special keys
        `callback` and `callbacks` are removed before processing.
    """

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        Callback.__init__(self, **kwargs)
        kwargs.pop("callback", None)
        kwargs.pop("callbacks", None)
        self._store: Dict[Any, Any] = dict(*args, **kwargs)
        self._reverse: Dict[Any, Any] = {v: k for k, v in self._store.items()}
        self._callbacks = []

    def __getitem__(self, key: Any) -> Any:
        """
        Retrieve the value associated with the given key.

        Parameters
        ----------
        key : hashable
            The key or value to retrieve its counterpart.

        Returns
        -------
        Any
            The corresponding value or key.

        Raises
        ------
        KeyError
            If the key is not found.
        """
        if key in self._store:
            return self._store[key]
        if key in self._reverse:
            return self._reverse[key]
        raise KeyError(key)

    def __setitem__(self, key: Any, value: Any) -> None:
        """
        Set a key-value pair in the dictionary, maintaining bidirectionality.

        Parameters
        ----------
        key : hashable
            The key to set.
        value : hashable
            The corresponding value.
        """
        if key in self._store:
            del self._reverse[self._store[key]]
        elif key in self._reverse:
            del self._store[self._reverse[key]]
        self._store[key] = value
        self._reverse[value] = key
        self._trigger_callbacks()

    def __delitem__(self, key: Any) -> None:
        """
        Delete a key-value pair from the dictionary.

        Parameters
        ----------
        key : hashable
            The key to remove.
        """
        value = self._store[key]
        del self._store[key]
        del self._reverse[value]
        self._trigger_callbacks()

    def __str__(self) -> str:
        """Return a string representation of the dictionary."""
        return str(self._store)

    def __repr__(self) -> str:
        """Return a string representation of the dictionary."""
        return str(self._store)

    def __contains__(self, key: Any) -> bool:
        """Check if the key or value exists in the dictionary."""
        return key in self._store or key in self._reverse

    def __iter__(self) -> Iterator[Any]:
        """Return an iterator over the keys of the dictionary."""
        return iter(self._store)

    def __len__(self) -> int:
        """Return the number of key-value pairs in the dictionary."""
        return len(self._store)

    def __eq__(self, other: Any) -> bool:
        """Check if two Basepair or dictionary have the same key-value pairs."""
        if not isinstance(other, (BasePair, dict)):
            return False
        return all(self[k] == other.get(k) for k in self.keys()) and all(
            other[k] == self.get(k) for k in other.keys()
        )

    ###
    ### PUBLIC METHODS
    ###

    def keys(self) -> Any:
        """Return the dictionary's keys."""
        return self._store.keys()

    def values(self) -> Any:
        """Return the dictionary's values."""
        return self._store.values()

    def items(self) -> Any:
        """Return the dictionary's key-value pairs."""
        return self._store.items()

    def update(self, *args: Dict[Any, Any], **kwargs: Any) -> "BasePair":
        """
        Update the dictionary with new key-value pairs.

        Parameters
        ----------
        *args : dict
            Dictionaries to merge into the current instance.
        **kwargs : dict
            Additional key-value pairs to update.

        Returns
        -------
        BasePair
            The updated instance.
        """
        for new_dict in args:
            for k, v in new_dict.items():
                self[k] = v
        for k, v in kwargs.items():
            self[k] = v
        self._trigger_callbacks()
        return self

    def get(self, key: Any, default: Any = None) -> Any:
        """
        Retrieve a value by key, returning a default if not found.

        Parameters
        ----------
        key : hashable
            The key to retrieve.
        default : Any, optional
            The value to return if key is not found (default is None).

        Returns
        -------
        Any
            The corresponding value or the default.
        """
        if key in self._store:
            return self._store[key]
        elif key in self._reverse:
            return self._reverse[key]
        return default

    def copy(self, **kwargs: Any) -> "BasePair":
        """
        Create a copy of the dictionary.

        Parameters
        ----------
        **kwargs : dict
            Additional parameters for the new instance.

        Returns
        -------
        BasePair
            A new instance of BasePair with the same data.
        """
        new_instance = BasePair(self._store, **kwargs)
        return new_instance

    def shift(self, shift: Any) -> "BasePair":
        """
        Apply a shift transformation to the key-value pairs.
        The shift is added to both the keys and values.

        Parameters
        ----------
        shift : hashable
            The amount to shift each key and value.

        Returns
        -------
        BasePair
            A new instance with shifted keys and values.
        """
        new_bp_dict = BasePair()
        for key, val in self.items():
            new_bp_dict[key + shift] = val + shift
        return new_bp_dict
