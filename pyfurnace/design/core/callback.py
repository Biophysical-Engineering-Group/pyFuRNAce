from typing import Callable, List, Any


class Callback:
    """
    A class to manage and trigger callback functions.

    This class allows the registration of functions (callbacks) that
    are executed when a specific event is triggered within the object.

    Attributes
    ----------
    _callbacks : List[Callable[..., None]]
        A list of registered callback functions.
    """

    def __init__(self, **kwargs: Any) -> None:
        """
        Initialize the Callback instance and register any initial callbacks if provided.

        Parameters
        ----------
        **kwargs : Any
            Arbitrary keyword arguments. Accepts:
            - 'callback': Callable[..., None]
                A single callback function to register.
            - 'callbacks': List[Callable[..., None]]
                A list of callback functions to register.
        """
        self._callbacks: List[Callable[..., None]] = []

        if "callback" in kwargs:
            self._callbacks.append(kwargs["callback"])  # Single callback registration

        if "callbacks" in kwargs:
            self._callbacks += kwargs["callbacks"]  # Multiple callback registration

    ###
    ### PROPERTIES
    ###

    @property
    def callbacks(self) -> List[Callable[..., None]]:
        """
        Get the list of registered callbacks.
        """
        return self._callbacks

    ###
    ### PROTECTED METHODS
    ###

    def _clear_callbacks(self) -> None:
        """
        Clear all registered callbacks from the object.

        Notes
        -----
        This method is useful for resetting the object's callback state,
        effectively removing all previously registered callbacks.
        """
        self._callbacks = []

    def _trigger_callbacks(self, **kwargs: Any) -> None:
        """
        Trigger all registered callback functions with the provided arguments.

        Parameters
        ----------
        **kwargs : Any
            Arbitrary keyword arguments passed to each callback function.

        Notes
        -----
        - If no callbacks are registered, the method simply returns.
        - Each callback will be executed with the provided keyword arguments.
        """
        if not hasattr(self, "_callbacks"):
            return

        for callback in self._callbacks:
            callback(**kwargs)

    ###
    ### PUBLIC METHODS
    ###

    def register_callback(self, callback: Callable[..., None]) -> None:
        """
        Register a new callback function to be executed when an event occurs.

        Parameters
        ----------
        callback : Callable[..., None]
            The callback function to register.

        Notes
        -----
        - A callback function should accept keyword arguments (`**kwargs`).
        - If the callback is already registered, it will not be added again.
        """
        if callback not in self._callbacks:
            self._callbacks.append(callback)
