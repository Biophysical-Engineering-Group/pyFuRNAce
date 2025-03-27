import os
import subprocess
import tempfile

from ..design.core.symbols import *
from .utils import find_stems_in_multiloop
from .pk_utils import parse_pseudoknots

def generate_road(structure, sequence, pseudoknots):
    try:
        # check if RNAfold is installed
        rnafold_path = subprocess.check_output(['which', 'RNAfold']).decode("utf-8").strip()
    except subprocess.CalledProcessError:
        rnafold_path = '/home/adminuser/.conda/bin/RNAfold'
        if not os.path.exists(rnafold_path):
            raise ValueError("RNAfold not found. Please install it or provide the path to the executable.")
        
    if '&' in sequence or '&' in structure:
        raise ValueError("The ROAD algorithm does not support multistranded structures.")
    
    ### PREPARE THE STUCTURE FILE

    ### Fix the short stems
    struct_list = list(structure)
    pair_map = dot_bracket_to_pair_map(structure)
    for dt in find_stems_in_multiloop(structure):
        # force pairing in dovetails that are shorter than 4
        if dt[-1] - dt[0] + 1 <= 3: 
            for i in range(dt[0], dt[1] + 1):
                struct_list[i] = '{'
                struct_list[pair_map[i]] = '}'

    ### Add the pseudoknots
    if type(pseudoknots) == dict:
        pk_dict = pseudoknots
    else:
        pk_dict = parse_pseudoknots(pseudoknots)

    road_pk_notation = {'A' : '1', 'B' : '2', 'C' : '3', 'D' : '4', 'E' : '5', 'F' : '6', 'G' : '7', 'H' : '8', 'I' : '9'}
    external_pk_count = 0
    for pk_info in pk_dict.values():
        used = False
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

    struct_list = ''.join(struct_list)
    print(''.join(struct_list))
    print(sequence)
        
    # create temporary directory to store files
    # with tempfile.TemporaryDirectory() as tmpdirname:
    #     with tempfile.NamedTemporaryFile(dir=self.temp_folder, suffix='.stl', delete=False) as temp_file:
    #         if isinstance(text, bytes):
    #             temp_file.write(text)
    #         elif isinstance(text, str):
    #             temp_file.write(text.encode("utf-8"))  # Write the text content to the file
    #         else:
    #             raise ValueError(f"Invalid text type for the stl file")
    #         temp_file.flush()  # Ensure all data is written to disk
    #         file_path = temp_file.name.split(os.sep)[-1]  # Store the relative path
    #         self.current_temp_files.append(temp_file.name)  # Keep track of the file for cleanup
    

# out = subprocess.run("cd road_bin; perl RNAbuild.pl pattern.txt", shell=True, capture_output=True)
# vrna_path = '/home/adminuser/.conda/bin/'
# subprocess.run(f"cd road_bin; perl trace_pattern.pl pattern.txt > target.txt", shell=True)
# subprocess.run(f"export PATH=$PATH:{vrna_path}; cd road_bin; perl batch_revolvr.pl 1", shell=True)
