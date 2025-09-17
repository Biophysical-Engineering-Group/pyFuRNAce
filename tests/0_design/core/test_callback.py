from pyfurnace.design.core.callback import Callback


def test_trigger_callbacks_early_return_when_attr_missing():
    """
    Cover line 78 -> 79: if self has no `_callbacks`, _trigger_callbacks returns early.
    """
    cb = Callback()
    # remove the attribute so hasattr(self, "_callbacks") is False
    del cb._callbacks
    # Should not raise and should just return None
    assert cb._trigger_callbacks(event="noop") is None


def test_init_with_single_and_multiple_callbacks_and_register_deduplicates():
    calls = []

    def a(**kw):
        calls.append(("a", kw))

    def b(**kw):
        calls.append(("b", kw))

    # Provide both 'callback' and 'callbacks' on init
    cb = Callback(callback=a, callbacks=[b])
    # Registering the same function again must not duplicate it
    cb.register_callback(a)
    cb.register_callback(b)

    assert cb.callbacks.count(a) == 1
    assert cb.callbacks.count(b) == 1

    # Fire them; order should be [a, b]
    cb._trigger_callbacks(x=1)
    assert calls == [("a", {"x": 1}), ("b", {"x": 1})]


def test_clear_callbacks_then_trigger_is_noop():
    fired = []

    def f(**kw):
        fired.append(kw)

    cb = Callback(callback=f)
    cb._clear_callbacks()
    cb._trigger_callbacks(x=42)  # nothing should happen
    assert fired == []
    assert cb.callbacks == []
