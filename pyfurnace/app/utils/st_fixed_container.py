from typing import Literal
import streamlit as st

# Solution adapted from:
# https://discuss.streamlit.io/t/is-it-possible-to-create-a-sticky-header/33145/3

MARGINS = {
    "top": "2.875rem",
    "bottom": "1rem",
}

STICKY_CONTAINER_HTML = """
<style>
div[data-testid="stVerticalBlock"] div:has(div.fixed-header-{i}) {{
    position: sticky;
    {position}: {margin};
    background-color: white;
    z-index: 1000000;
}}
</style>
<div class='fixed-header-{i}'/>
""".strip()

# Not to apply the same style to multiple containers
count = 0


def sticky_container(
    *,
    # height: str = 'content',
    border: bool = False,
    mode: Literal["top", "bottom"] = "top",
    margin: str | None = None,
):
    if margin is None:
        margin = MARGINS[mode]

    global count
    html_code = STICKY_CONTAINER_HTML.format(position=mode, margin=margin, i=count)
    count += 1

    # I never use the height option, so commented out for now
    # it could be added later if needed, but changes in streamlit versions
    container = st.container(border=border)  # height=height,
    container.markdown(html_code, unsafe_allow_html=True)
    return container


if __name__ == "__main__":
    st.title("Sticky Container")
    st.write("This is a demonstration of a sticky container.")
    for i in range(30):
        st.write(f"Line {i}")
    with sticky_container(mode="top", border=True):
        st.write("This is a sticky container at the top.")
        st.write("This is a sticky container at the top.")
    for i in range(30):
        st.write(f"Line {i}")
    with sticky_container(mode="bottom", border=True):
        st.write("This is a sticky container at the bottom.")
        st.write("This is a sticky container at the bottom.")
    for i in range(30):
        st.write(f"Line {i}")
