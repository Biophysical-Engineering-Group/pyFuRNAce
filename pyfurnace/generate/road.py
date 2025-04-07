import os
import sys
import subprocess
import tempfile
import zipfile
import shutil
import warnings

from ..design.core.symbols import *
from .utils import find_stems_in_multiloop
from .pk_utils import parse_pseudoknots

def generate_road(structure, 
                  sequence, 
                  pseudoknots, 
                  name='origami', 
                  callback=None, 
                  directory=None, 
                  zip_directory=False,
                  origami_code : str = None):
    
    road_dir = __file__.replace('road.py', 'road_bin')
    
    files_to_include = []
        
    # Sanity check
    if '&' in sequence or '&' in structure:
        raise ValueError("The ROAD algorithm does not support multistranded structures.")
    
    ### Fix the short stems with '{' ROAD symbols
    struct_list = list(structure)
    pair_map = dot_bracket_to_pair_map(structure)
    for dt in find_stems_in_multiloop(structure):
        # force pairing in dovetails that are shorter than 3
        if dt[-1] - dt[0] + 1 <= 2: 
            for i in range(dt[0], dt[1] + 1):
                struct_list[i] = '{'
                struct_list[pair_map[i]] = '}'

    ### ADD THE PSEUDOKNOTS ROAD NOTATION
    if type(pseudoknots) == dict:
        pk_dict = pseudoknots
    else:
        pk_dict = parse_pseudoknots(pseudoknots)

    road_pk_notation = {'A' : '1', 'B' : '2', 'C' : '3', 'D' : '4', 'E' : '5', 
                        'F' : '6', 'G' : '7', 'H' : '8', 'I' : '9'}
    external_pk_count = 0
    avg_pk_E = 0
    avg_pk_dE = 0
    for pk_info in pk_dict.values():
        used = False
        avg_pk_E += pk_info['E']
        avg_pk_dE += abs(pk_info['dE'])
        pk_sym = list(road_pk_notation)[external_pk_count]
        for (start, end) in pk_info['ind_fwd']:
            if struct_list[start] not in '.()':
                continue
            for i in range(start, end + 1):
                struct_list[i] = pk_sym
            used = True
        for (start, end) in pk_info['ind_rev']:
            if struct_list[start] not in '.()':
                continue
            for j in range(start, end + 1):
                struct_list[j] = road_pk_notation[pk_sym]
            used = True
        if used:
            external_pk_count += 1

    if pk_dict:
        avg_pk_E /= len(pk_dict)
        avg_pk_dE /= len(pk_dict)
    structure = ''.join(struct_list)

    ### COPY THE PYTHON PATH AND ADD RNAfold
    python_path = sys.executable # Get path to the current Python interpreter
    python_dir = os.path.dirname(python_path)

    # Prepend it to PATH
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
        directory = directory

    ### WORK IN THE DIRECTORY

    # save the origami code
    if origami_code is not None:
        origami_code_path = os.path.join(directory, f'{name}.py')
        with open(origami_code_path, 'w') as f:
            f.write(origami_code)
        # include it in the zip
        files_to_include.append(origami_code_path)
        
    # create the target input file
    target_path = os.path.join(directory, 'target.txt')
    with open(target_path, 'w') as f:
        f.write(f"{name}\n{structure}\n{sequence}\n")
    # include it in the zip
    files_to_include.append(target_path)

    # read the revolvr file
    revolvr_local_path = os.path.join(road_dir, 'revolvr.pl')
    with open(revolvr_local_path, 'r') as f:
        revolvr_text = f.read()
        
    # replace the KL energy parameters
    revolvr_text = revolvr_text.replace('my $MinKL = -7.2;',
                                        f'my $MinKL = {avg_pk_E + avg_pk_dE};')
    revolvr_text = revolvr_text.replace('my $MaxKL = -10.8;',
                                        f'my $MaxKL = {avg_pk_E - avg_pk_dE};')

    # create the revolvr file with specific KL parameters
    out_revolvr = os.path.join(directory, 'revolvr.pl')
    with open(out_revolvr, 'w') as f:
        f.write(revolvr_text)
    # include it in the zip
    files_to_include.append(out_revolvr)

    vienna_out_path = os.path.join(directory, 'viennarna_funcs.py')
    shutil.copyfile(os.path.join(road_dir, 'viennarna_funcs.py'),
                    vienna_out_path
                    )
    # include it in the zip
    files_to_include.append(vienna_out_path)
    
    command = f'perl "{out_revolvr}" "{directory}"'
    process = subprocess.Popen(command,
                               shell=True,
                               cwd=directory,
                               env=env,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT,
                               text=True  #makes output strings
                               )

    # Read output in real time
    last_seq = '' 
    last_struct = ''
    prev_line = ''
    n_stage = 0
    stages = ['Designing', 
              'GC-Reduction & GC/UA-Rich Reduction', 
              'Mapping Kissing Loops', 
              'Optimization']

    for line in process.stdout:
        line = line.strip()

        # update the stage
        if stages[n_stage] in line and n_stage < len(stages) - 1:
            n_stage += 1

        # update the last sequence and structure
        if prev_line and not prev_line.translate(nucl_to_none):
            last_seq = prev_line
            if line and any(s in line for s in '.()'):
                last_struct = line
        
        if line and last_seq and last_struct and callback:
            callback(last_struct, 
                     last_seq, 
                     line, 
                     stages[n_stage - 1], 
                     n_stage / len(stages))

        prev_line = line

    # wait for process to finish
    process.wait()

    # spool file
    files_to_include.append(os.path.join(directory, f'{name}_spool.txt'))
    
    # read the results
    design_path = os.path.join(directory, f'{name}_design.txt')
    with open(os.path.join(directory, f'{name}_design.txt'), 'r') as f:
        lines = f.readlines()
        last_seq = lines[2].strip()

    # include it in the zip
    files_to_include.append(design_path)

    if zip_directory:
        # create a temporary zip file
        temp_zip = tempfile.NamedTemporaryFile(suffix=".zip", delete=False)
        temp_zip.close()  # Close so we can write to it

        # Create the zip and add selected files
        with zipfile.ZipFile(temp_zip.name, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for file in files_to_include:
                # Store relative to directory for cleaner archive structure
                zipf.write(file, arcname=os.path.basename(file))

    if tempdir is not None:
        # Close the temporary directory
        tempdir.cleanup()
    
    # Return the path to the zip file
    if zip_directory:
        return last_seq, temp_zip.name
    
    return last_seq
