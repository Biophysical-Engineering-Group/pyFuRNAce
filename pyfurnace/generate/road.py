import os
import sys
import subprocess
import tempfile
import zipfile
import shutil
import warnings
from contextlib import contextmanager
from typing import Callable, Optional, Tuple, Union
import multiprocessing

### Local imports
from ..design.core.symbols import dot_bracket_to_pair_map, nucl_to_none
from .utils import find_stems_in_multiloop
from .pk_utils import parse_pseudoknots
from .viennarna import fold_p


@contextmanager
def _road_cwd(path: str):
    """Context manager to change the current working directory and revert it back."""
    oldpwd = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(oldpwd)


def generate_road(
    structure: str,
    sequence: str,
    pseudoknots: Union[str, dict] = "",
    name: str = "origami",
    initial_sequence: Optional[str] = None,
    callback: Optional[Callable[[str, str, str, str, float], None]] = None,
    verbose: bool = False,
    timeout: int = 7200,
    directory: Optional[str] = None,
    zip_directory: bool = False,
    origami_code: Optional[str] = None,
) -> Union[str, Tuple[str, str]]:
    """
    Generate an RNA origami design using the ROAD Revolvr algorithm.

    This function orchestrates the structure processing, pseudoknot embedding,
    and Perl-based design process for RNA origami, with optional support for
    real-time updates and zipped output.

    Parameters
    ----------
    structure : str
        RNA secondary structure in dot-bracket notation.
    sequence : str
        Initial sequence to use for the design (must match structure length).
    pseudoknots : str or dict, optional
        Pseudoknot definitions in ROAD-compatible format (string or parsed dict).
    name : str, optional
        Base name for the design file (default is 'origami').
    initial_sequence : str, optional
        Optional initial sequence to pre-fill as the starting point of the
        optimization.
    callback : callable, optional
        A function called during ROAD execution with progress updates:
        `callback(structure, sequence, line, stage_name, stage_progress)`
    verbose : bool, optional
        If True, prints progress updates to the console.
        If `callback` is provided, this is ignored. Default is False.
    timeout : int, optional
        Timeout in seconds for the ROAD optimization
        (default is 7200 seconds = 2 hours).
    directory : str, optional
        Directory where files should be written. If None, a temporary directory
        is used.
    zip_directory : bool, optional
        Whether to generate a `.zip` file with all intermediate and result files.
    origami_code : str, optional
        Source code (e.g., Python) representing the origami design, to be saved
        and zipped.

    Returns
    -------
    str or tuple
        If `zip_directory` is False: returns the designed RNA sequence as a string.
        If `zip_directory` is True: returns
        the tuple `(designed_sequence, path_to_zip_file)`.

    Raises
    ------
    ValueError
        If the structure and sequence lengths do not match, or if multistrand
        '&' is used.

    Warnings
    --------
    - If the working directory doesn't exist, it is created.
    - If the final design file is missing, a warning is issued.
    """

    road_dir = __file__.replace("road.py", "road_bin")

    files_to_include = []

    # Sanity check
    if "&" in sequence or "&" in structure:
        raise ValueError(
            "The ROAD Revolvr algorithm does not support multistranded structures."
        )

    if len(sequence) != len(structure):
        raise ValueError(
            "The length of the sequence and structure must match."
            " Got {} and {}.".format(len(sequence), len(structure))
        )

    if verbose and callback is None:

        def callback(_, __, spool_line, stage_name, *args):
            print(f"{stage_name}: {spool_line.strip()}", flush=True)

    ### Fix the short stems with '{' ROAD symbols
    struct_list = list(structure)
    pair_map = dot_bracket_to_pair_map(structure)
    for dt in find_stems_in_multiloop(structure):
        # force pairing in dovetails that are shorter than 3
        if dt[-1] - dt[0] + 1 <= 2:
            for i in range(dt[0], dt[1] + 1):
                struct_list[i] = "{"
                struct_list[pair_map[i]] = "}"

    ### ADD THE PSEUDOKNOTS ROAD NOTATION
    avg_pk_E = 0
    avg_pk_dE = 0
    if not pseudoknots:
        avg_pk_E = -9.0
        avg_pk_dE = 1.81

    if isinstance(pseudoknots, dict):
        pk_dict = pseudoknots
    else:
        pk_dict = parse_pseudoknots(pseudoknots)

    road_pk_notation = {
        "A": "1",
        "B": "2",
        "C": "3",
        "D": "4",
        "E": "5",
        "F": "6",
        "G": "7",
        "H": "8",
        "I": "9",
    }
    external_pk_count = 0

    # sort the pseudoknots by their index
    pk_dict = {
        k: v
        for k, v in sorted(
            pk_dict.items(),
            key=lambda item: min(
                inds[0] for inds in item[1]["ind_fwd"] + item[1]["ind_rev"]
            ),
        )
    }
    # Calculate the average energy and dE of the pseudoknots
    for pk_info in pk_dict.values():
        used = False
        avg_pk_E += pk_info["E"]
        avg_pk_dE += abs(pk_info["dE"])
        pk_sym = list(road_pk_notation)[external_pk_count]

        # Find the pk with the lowest index, between fwd and rev
        min_fwd = min([inds[0] for inds in pk_info["ind_fwd"]], default=float("inf"))
        min_rev = min([inds[0] for inds in pk_info["ind_rev"]], default=float("inf"))
        if min_fwd < min_rev:
            first_pk = pk_info["ind_fwd"]
            second_pk = pk_info["ind_rev"]
        else:
            first_pk = pk_info["ind_rev"]
            second_pk = pk_info["ind_fwd"]

        for start, end in first_pk:
            if struct_list[start] not in ".()":
                continue
            for i in range(start, end + 1):
                struct_list[i] = pk_sym
            used = True
        for start, end in second_pk:
            if struct_list[start] not in ".()":
                continue
            for j in range(start, end + 1):
                struct_list[j] = road_pk_notation[pk_sym]
            used = True
        if used:
            external_pk_count += 1

    if pk_dict:
        avg_pk_E /= len(pk_dict)
        avg_pk_dE /= len(pk_dict)
    structure = "".join(struct_list)

    ### COPY THE PYTHON PATH AND ADD RNAfold
    python_path = sys.executable  # Get path to the current Python interpreter
    python_dir = os.path.dirname(python_path)

    # Prepend it to PATH, to make sure the correct Python is used
    env = os.environ.copy()
    env["PATH"] = python_dir + os.pathsep + os.environ["PATH"]

    ### CHECK THE TEMPORARY DIRECTORY
    if directory is None:
        tempdir = tempfile.TemporaryDirectory()
        directory = tempdir.name
    else:
        if not os.path.exists(directory):
            warnings.warn(f"Directory {directory} does not exist. Creating it.")
            os.makedirs(directory)
        tempdir = None
        # take the absolute path of the directory
        directory = os.path.abspath(directory)

    ### WORK IN THE DIRECTORY

    # save the origami code
    if origami_code is not None:
        origami_code_path = os.path.join(directory, f"{name}.py")
        with open(origami_code_path, "w", encoding="utf-8") as f:
            f.write(origami_code)
        # include it in the zip
        files_to_include.append(origami_code_path)

    # create the target input file
    target_path = os.path.join(directory, "target.txt")
    with open(target_path, "w", encoding="utf-8") as f:
        f.write(f"{name}\n{structure}\n{sequence}\n")

        if initial_sequence is not None and len(initial_sequence) == len(structure):
            f.write(f"{initial_sequence}\n")

    # include it in the zip
    files_to_include.append(target_path)

    # read the revolvr file
    revolvr_local_path = os.path.join(road_dir, "revolvr.pl")
    with open(revolvr_local_path, "r") as f:
        revolvr_text = f.read()

    # replace the KL energy parameters
    revolvr_text = revolvr_text.replace(
        "my $MinKL = -7.2;", f"my $MinKL = {avg_pk_E + avg_pk_dE};"
    )
    revolvr_text = revolvr_text.replace(
        "my $MaxKL = -10.8;", f"my $MaxKL = {avg_pk_E - avg_pk_dE};"
    )
    revolvr_text = revolvr_text.replace(
        "my $timeout_length = 7200;", f"my $timeout_length = {int(timeout)};"
    )

    # create the revolvr file with specific KL parameters
    out_revolvr = os.path.join(directory, "revolvr.pl")
    with open(out_revolvr, "w", encoding="utf-8") as f:
        f.write(revolvr_text)
    # include it in the zip
    files_to_include.append(out_revolvr)

    vienna_out_path = os.path.join(directory, "viennarna_funcs.py")
    shutil.copyfile(os.path.join(road_dir, "viennarna_funcs.py"), vienna_out_path)
    # include it in the zip
    files_to_include.append(vienna_out_path)

    # move to the directory
    with _road_cwd(directory):
        command = f'perl revolvr.pl "{directory}"'
        process = subprocess.Popen(
            command,
            shell=True,
            cwd=directory,
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,  # makes output strings
        )

        # Read output in real time
        last_seq = ""
        last_struct = ""
        prev_line = ""
        n_stage = 0
        stages = [
            "Designing",
            "GC-Reduction & GC/UA-Rich Reduction",
            "Mapping Kissing Loops",
            "Optimization",
        ]

        for line in process.stdout:
            line = line.strip()

            # update the stage
            if stages[n_stage] in line and n_stage < len(stages) - 1:
                n_stage += 1

            # update the last sequence and structure
            if prev_line and not prev_line.translate(nucl_to_none):
                last_seq = prev_line
                if line and any(s in line for s in ".()"):
                    last_struct = line

            if line and last_seq and last_struct and callback:
                callback(
                    last_struct,
                    last_seq,
                    line,
                    stages[n_stage - 1],
                    n_stage / len(stages),
                )

            prev_line = line

        # wait for process to finish
        process.wait()

    # spool file
    files_to_include.append(os.path.join(directory, f"{name}_spool.txt"))

    # read the results
    design_path = os.path.join(directory, f"{name}_design.txt")
    if not os.path.exists(design_path):
        last_seq = ""
        warnings.warn(
            f"Design file {design_path} not found. Optimization failed. "
            "Please check the ROAD algorithm output."
        )
    else:
        with open(os.path.join(directory, f"{name}_design.txt"), "r") as f:
            lines = f.readlines()
            last_seq = lines[2].strip()

    # include it in the zip
    files_to_include.append(design_path)

    if zip_directory:
        # create a temporary zip file
        temp_zip = tempfile.NamedTemporaryFile(suffix=".zip", delete=False)
        temp_zip.close()  # Close so we can write to it

        # Create the zip and add selected files
        with zipfile.ZipFile(temp_zip.name, "w", zipfile.ZIP_DEFLATED) as zipf:
            for file in files_to_include:
                try:
                    # Store relative to directory for cleaner archive structure
                    zipf.write(file, arcname=os.path.basename(file))
                except Exception as e:
                    warnings.warn(f"Failed to add {file} to zip: {e}")

    if tempdir is not None:
        # Close the temporary directory
        tempdir.cleanup()

    # Return the path to the zip file
    if zip_directory:
        return last_seq, temp_zip.name

    return last_seq


def _worker_single_road(
    structure: str,
    sequence: str,
    pseudoknots,
    name: str,
    callback: callable = None,
    verbose: bool = False,
    timeout: int = 7200,
    zip_directory: bool = False,
    origami_code: Optional[str] = None,
    initial_sequence: Optional[str] = None,
    result_queue: multiprocessing.Queue = None,
    stop_event: "multiprocessing.Event" = None,
):
    """
    Internal worker function that runs a single ROAD design trial.

    This function wraps `generate_road` for parallel execution. It checks
    whether an early stop is requested via `stop_event`, runs the design
    pipeline, and returns the result via a shared `result_queue`.

    Parameters
    ----------
    structure : str
        RNA secondary structure in dot-bracket notation.
    sequence : str
        Input sequence, same length as the structure.
    pseudoknots : str or dict
        Pseudoknot annotations in ROAD-compatible format.
    name : str
        Base name for generated files.
    callback : callable, optional
        Optional function for progress updates.
    verbose : bool, optional
        If True, prints updates to the console.
    timeout : int, optional
        Max time in seconds to allow ROAD to run.
    zip_directory : bool, optional
        Whether to store intermediate and output files in a ZIP archive.
    origami_code : str, optional
        Python code to be stored with the output files.
    initial_sequence : str, optional
        Optional prefilled sequence as the starting point.
    result_queue : multiprocessing.Queue, optional
        Queue for storing results from each trial.
    stop_event : multiprocessing.Event, optional
        Used to signal early stopping when a successful result is found.
    """
    try:
        if stop_event.is_set():
            return

        result = generate_road(
            structure=structure,
            sequence=sequence,
            pseudoknots=pseudoknots,
            name=name,
            initial_sequence=initial_sequence,
            callback=callback,
            verbose=verbose,
            timeout=timeout,
            zip_directory=zip_directory,
            origami_code=origami_code,
        )

        if result:
            result_queue.put(result)
            stop_event.set()

    except Exception as e:
        print(f"Error in ROAD design process: {e}", file=sys.stderr)


def parallel_road(
    structure: str,
    sequence: str,
    pseudoknots="",
    name: str = "origami",
    callback=None,
    verbose=False,
    timeout: int = 7200,
    zip_directory: bool = True,
    origami_code: str = None,
    initial_sequence: str = None,
    n_trials: int = 8,
    save_to: Optional[str] = None,
    wait_for_all: bool = False,
) -> Union[str, Tuple[str, str]]:
    """
    Run multiple ROAD optimization trials in parallel and return results.

    You can choose between returning early on the first successful result
    or waiting for all trials to complete. Results can optionally be
    ZIP-archived and saved. If the wait_for_all flag is set to True,
    the sequences are returned sorted by their MFE frequency in the ensemble
    (highest to lowest).

    Parameters
    ----------
    structure : str
        Dot-bracket notation of the target RNA secondary structure.
    sequence : str
        Reference sequence with the same length as the structure.
    pseudoknots : str or dict, optional
        ROAD-style or pyFuRNAce-style pseudoknot definition.
    name : str, optional
        Base name for output files (default: "origami").
    callback : callable, optional
        Optional function that receives real-time progress updates.
    verbose : bool, optional
        If True, prints progress to console.
    timeout : int, optional
        Timeout in seconds per trial (default: 7200).
    zip_directory : bool, optional
        Whether to save each trial's output as a ZIP archive.
    origami_code : str, optional
        Python source code to include with the design output.
    initial_sequence : str, optional
        Optional sequence to use as a starting point in optimization.
    n_trials : int, optional
        Number of parallel optimization trials to run.
    save_to : str, optional
        Directory to store ZIP outputs if `zip_directory=True`.
    wait_for_all : bool, optional
        If True, waits for all trials to finish and returns results sorted
        by MFE frequency. If False, returns as soon as the first trial succeeds.

    Returns
    -------
    str or Tuple[str, str] or List[str]
        - If `wait_for_all=False` returns the first designed sequence.
        - If `wait_for_all=True` returns a list of designed sequences.

    Raises
    ------
    RuntimeError
        If no trials successfully generate a design.
    """
    manager = multiprocessing.Manager()
    result_queue = manager.Queue()
    stop_event = manager.Event()
    processes = []

    for i in range(n_trials):
        p = multiprocessing.Process(
            target=_worker_single_road,
            args=(
                structure,
                sequence,
                pseudoknots,
                name,
                callback,
                verbose,
                timeout,
                zip_directory,
                origami_code,
                initial_sequence,
                result_queue,
                stop_event,
            ),
        )
        processes.append(p)
        p.start()

    results = []

    try:
        if wait_for_all:
            # Wait for all trials to finish and collect all successful results
            for _ in range(n_trials):
                try:
                    result = result_queue.get(timeout=timeout)
                    results.append(result)
                except Exception:  # likely empty queue error
                    pass  # skip failed or timed-out results
        else:
            # Wait for the first successful result
            try:
                result = result_queue.get(timeout=timeout)
                results.append(result)
            except Exception:  # likely empty queue error
                pass  # no result within timeout
            stop_event.set()
    finally:
        for p in processes:
            p.terminate()
        for p in processes:
            p.join()

    if not results:
        raise RuntimeError("No successful ROAD design found in any trial.")

    # order the sequences from highest frequency of MFE in the ensemble to lowest
    results = sorted(results, key=lambda x: -fold_p(x[0])[2])

    if zip_directory and isinstance(result, tuple) and save_to:
        # If saving zips and save_to is set, move them all
        if not os.path.exists(save_to):
            os.makedirs(save_to)
            warnings.warn(f"Directory {save_to} did not exist. Created it.")
        for i, res in enumerate(results):
            if isinstance(res, tuple):
                final_zip_path = os.path.join(save_to, f"{name}_trial{i + 1}.zip")
                shutil.move(res[1], final_zip_path)

        if wait_for_all:
            # return all sequences in a list
            return [res[0] for res in results if isinstance(res, tuple)]
        else:
            # return the first sequence
            return results[0][0]

    # Return appropriately
    if wait_for_all:
        return results
    else:
        return results[0]
