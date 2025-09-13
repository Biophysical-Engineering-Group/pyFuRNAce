from pathlib import Path
import pytest
from streamlit.testing.v1.app_test import AppTest

# repo/
# ├─ pyfurnace/
# │  └─ app/
# │     ├─ streamlit_app.py
# │     ├─ utils.py           <-- if pages do "from utils import ..."
# │     └─ pages/...
# └─ tests/
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
APP_PATH = REPO_ROOT / "pyfurnace" / "app"

# Pages to smoke-test
PAGES = [
    "pages/0_Home.py",
    "pages/1_Design.py",
    "pages/2_Generate.py",
    "pages/3_Convert.py",
    "pages/4_Prepare.py",
]


@pytest.fixture(autouse=True)
def run_from_app_path(monkeypatch):
    """
    Ensure we run *from* the app path and that imports resolve the same
    way they do when launching Streamlit manually.
    """
    # 1) CWD = app dir (so relative paths in the app work)
    monkeypatch.chdir(APP_PATH)

    # 2) sys.path: allow "import utils" (module in app dir)
    monkeypatch.syspath_prepend(str(APP_PATH))

    # 3) sys.path: allow absolute package imports like "from pyfurnace...."
    monkeypatch.syspath_prepend(str(REPO_ROOT))


def test_home_page_renders_without_exception():
    at = AppTest.from_file("streamlit_app.py").run()

    assert not at.exception, f"Streamlit raised: {at.exception}"
    assert at.markdown[0].value == "# Hello and Welcome to pyFuRNAce!"


@pytest.mark.parametrize("page_path", PAGES)
def test_each_page_loads_without_exception(page_path):
    # Launch main app
    at = AppTest.from_file("streamlit_app.py").run()

    # Switch page and re-run
    at.switch_page(page_path)
    at.run()

    assert not at.exception, f"Exception on {page_path}: {at.exception}"
