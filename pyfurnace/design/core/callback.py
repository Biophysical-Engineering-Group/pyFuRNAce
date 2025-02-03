
class Callback:
    
    def __init__(self, **kwargs):
        self._callbacks = []
        if "callback" in kwargs:
            self._callbacks.append(kwargs["callback"])
        if 'callbacks' in kwargs:
            self._callbacks += kwargs['callbacks']

    @property
    def callbacks(self):
        return self._callbacks

    def register_callback(self, callback):
        """Register a callback function to be called when the object is modified."""
        if callback not in self._callbacks:
            self._callbacks.append(callback)

    def _trigger_callbacks(self, **kwargs):
        """Trigger all registered callbacks."""
        if not hasattr(self, '_callbacks'):
            return
        for callback in self._callbacks:
            callback(**kwargs)
        
    def _clear_callbacks(self):
        self._callbacks = []