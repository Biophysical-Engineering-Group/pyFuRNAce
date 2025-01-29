import os
import pkgutil
import inspect
from importlib import import_module

# Get the current package's directory
current_package = __name__
current_path = os.path.dirname(__file__)

# Iterate through all modules in the current package
for module_info in pkgutil.iter_modules([current_path]):
    module_name = module_info.name  # Get the module name
    if module_name == "__init__":
        continue

    # Import the module
    module = import_module(f".{module_name}", package=current_package)

    # Iterate through all classes in the module
    for name, obj in inspect.getmembers(module, inspect.isclass):
        # Check if the class ends with "Command" and is defined in the current module
        if name.endswith("Command") and obj.__module__ == module.__name__:
            # Add the class to the package's namespace
            globals()[name] = obj
