import pathlib
import subprocess

def main():
    config_options = {'theme.base': 'light',
                      'theme.font': "monospace",
                      'theme.primaryColor': "#00856A",
                      'server.enableStaticServing': True,
                      'logger.level': 'error',
                      'global.showWarningOnDirectExecution': False,
                      'client.toolbarMode': "minimal",
                    }
    config_str = ' '.join([f'--{key} "{value}"' for key, value in config_options.items()])
    app_path = str(pathlib.Path(__file__).parent.resolve() / 'app' / 'Home.py')
    subprocess.run(f'streamlit run {app_path} {config_str}', shell=True)

if __name__ == "__main__":
    main()