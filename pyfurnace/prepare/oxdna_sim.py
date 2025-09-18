from typing import Union
import os
import zipfile
import tempfile
import warnings
from ..design.core import Origami


def oxdna_simulations(
    origami: Origami,
    oxdna_directory: str,
    directory: Union[None, str] = None,
    zip_directory: bool = False,
    temperature: Union[int, float] = 37,
    mc_relax_steps: Union[int, float] = 5e3,
    md_relax_steps: Union[int, float] = 1e7,
    md_equil_steps: Union[int, float] = 1e8,
    md_prod_steps: Union[int, float] = 1e9,
) -> Union[str, None]:
    """
    Prepare the files for the OxDNA simulations. The function saves the
    files for oxdna simulations in the specified directory (created if it does not
    exist). The function also creates a zip file with the files for the simulations
    if the flag zip_directory is set to True. If the directory is set to None, the
    program automatically creates a zip file and returns the path to the zip file.

    Parameters
    ----------
    origami : Origami
        The origami object.
    oxdna_directory : str
        The path to the OxDNA executable. It is important to load sequence dependent
        parameters.
    directory : str
        The output directory for the MD simulations files. If None, the zip_directory
        flag is set to True.
    zip_directory : bool
        If True, the function creates a zip file with the files for the simulations
        and returns the path to the zip file. If False, the function returns None.
    temperature : Union[int, float], optional
        The temperature in Celsius, by default 37.
    mc_relax_steps : Union[int, float], optional
        The number of Monte Carlo relaxation steps, by default 5e3.
    md_relax_steps : Union[int, float], optional
        The number of MD relaxation steps, by default 1e7.
    md_equil_steps : Union[int, float], optional
        The number of MD equilibration steps, by default 1e8.
    md_prod_steps : Union[int, float], optional
        The number of MD production steps, by default 1e9.

    Returns
    -------
    Union[str, None]
        If zip_directory is True or directory is None, the function returns the path
        to the zip file. If zip_directory is False and directory is not None, the
        function returns None.
    """

    inputs_dir = __file__.replace("oxdna_sim.py", "oxdna_inputs")

    ### CHECK THE TEMPORARY DIRECTORY
    if directory is None:
        zip_directory = True
        tempdir = tempfile.TemporaryDirectory()
        directory = tempdir.name
    else:
        if not os.path.exists(directory):
            warnings.warn(f"Directory {directory} does not exist. Creating it.")
            os.makedirs(directory)
        tempdir = None
        directory = directory

    # save configuration and forces
    origami.save_3d_model(
        directory + os.sep + "start",
        config=True,
        topology=False,
        forces=True,
        pk_forces=True,
    )

    # rename the forces file
    forces_file = os.path.join(directory, "start_forces.txt")
    new_forces_file = os.path.join(directory, "forces.txt")
    os.rename(forces_file, new_forces_file)

    # rename the pk_forces file
    pk_forces_file = os.path.join(directory, "start_pk_forces.txt")
    new_pk_forces_file = os.path.join(directory, "pk_forces.txt")
    os.rename(pk_forces_file, new_pk_forces_file)

    # save topology
    origami.save_3d_model(
        directory + os.sep + "topology",
        config=False,
        topology=True,
    )

    for file in os.listdir(inputs_dir):
        path = os.path.join(inputs_dir, file)
        with open(path, "r") as f:
            text = f.read()
        text = text.replace("TEMPERATURE", str(temperature))
        text = text.replace("OXDNA_DIRECTORY", oxdna_directory)

        match file:
            case "MC_relax.txt":
                steps = mc_relax_steps
            case "MD_relax.txt":
                steps = md_relax_steps
            case "MD_equil.txt":
                steps = md_equil_steps
            case "MD_prod.txt":
                steps = md_prod_steps
        text = text.replace("STEPS", str(int(steps)))

        with open(os.path.join(directory, file), "w", encoding="utf-8") as f:
            f.write(text)

    files_to_include = [
        "start.dat",
        "topology.top",
        "forces.txt",
        "pk_forces.txt",
        "MC_relax.txt",
        "MD_relax.txt",
        "MD_equil.txt",
        "MD_prod.txt",
    ]

    if zip_directory:
        # create a temporary zip file
        temp_zip = tempfile.NamedTemporaryFile(suffix=".zip", delete=False)
        temp_zip.close()  # Close so we can write to it

        # Create the zip and add selected files
        with zipfile.ZipFile(temp_zip.name, "w", zipfile.ZIP_DEFLATED) as zipf:
            for file in files_to_include:
                # Store relative to directory for cleaner archive structure
                zipf.write(os.path.join(directory, file), arcname=file)

    if tempdir is not None:
        # Close the temporary directory
        tempdir.cleanup()

    # Return the path to the zip file
    if zip_directory:
        return temp_zip.name
