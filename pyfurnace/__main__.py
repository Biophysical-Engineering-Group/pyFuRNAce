import os
import sys
import pathlib
import subprocess


def main():
    config_options = {
        "theme.base": "light",
        "theme.font": "monospace",
        "theme.primaryColor": "#00856A",
        "server.enableStaticServing": True,
        "logger.level": "error",
        "global.showWarningOnDirectExecution": False,
        "client.toolbarMode": "minimal",
    }
    config_str = " ".join(
        [f'--{key} "{value}"' for key, value in config_options.items()]
    )

    app_path = str(pathlib.Path(__file__).parent.resolve() / "app" / "streamlit_app.py")

    ### COPY THE PYTHON PATH AND ADD RNAfold
    python_path = sys.executable  # Get path to the current Python interpreter
    python_dir = os.path.dirname(python_path)

    # Prepend it to PATH, to make sure the correct Python is used
    env = os.environ.copy()
    env["PATH"] = python_dir + os.pathsep + os.environ["PATH"]

    subprocess.run(f'streamlit run "{app_path}" {config_str}', shell=True, env=env)


if __name__ == "__main__":
    main()
