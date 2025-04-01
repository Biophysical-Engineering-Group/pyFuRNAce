import os
import sys
import subprocess
import tempfile
import shutil
import warnings

from ..design.core.symbols import *
from .utils import find_stems_in_multiloop
from .pk_utils import parse_pseudoknots

def generate_road(structure, sequence, pseudoknots, name='origami', callback=None, directory=None, zip_directory=False):
    road_dir = __file__.replace('road.py', 'road_bin')
        
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
    print(f"Python path: {python_path}")
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
        
    # create the target input file
    with open(os.path.join(directory, 'target.txt'), 'w') as f:
        f.write(f"{name}\n{structure}\n{sequence}\n")

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

    shutil.copyfile(os.path.join(road_dir, 'viennarna_funcs.py'),
                    os.path.join(directory, 'viennarna_funcs.py')
                    )
    
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
        print(line)

        # update the stage
        if stages[n_stage] in line and n_stage < len(stages) - 1:
            n_stage += 1

        # update the last sequence and structure
        if prev_line and not prev_line.translate(nucl_to_none):
            last_seq = prev_line
            if line and any(s in line for s in '.()'):
                last_struct = line
        
        if line and last_seq and last_struct and callback:
            callback(last_struct, last_seq, line, stages[n_stage - 1], n_stage / len(stages))

        prev_line = line

    # Wait for process to finish
    process.wait()

    try:
        with open(os.path.join(directory, f'{name}_design.txt'), 'r') as f:
            lines = f.readlines()
            last_seq = lines[2].strip()
    except:
        return last_seq, structure, None


    # create a temporary zip file with all the info
    temp_zip = tempfile.NamedTemporaryFile(suffix=".zip", delete=False)
    temp_zip.close()  # Close so shutil can write to it

    # Create a zip archive from the directory
    shutil.make_archive(base_name=temp_zip.name.replace(".zip", ""), 
                        format='zip', 
                        root_dir=directory)

    if tempdir is not None:
        # Close the temporary directory
        tempdir.cleanup()
    
    # Return the path to the zip file
    if zip_directory:
        return last_seq, temp_zip.name
    return last_seq


# out = subprocess.run("cd road_bin; perl RNAbuild.pl pattern.txt", shell=True, capture_output=True)
# vrna_path = '/home/adminuser/.conda/bin/'
# subprocess.run(f"cd road_bin; perl trace_pattern.pl pattern.txt > target.txt", shell=True)
# subprocess.run(f"export PATH=$PATH:{vrna_path}; cd road_bin; perl batch_revolvr.pl 1", shell=True)
