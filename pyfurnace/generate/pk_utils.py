import ast
import warnings
from ..design import db_pairs, dot_bracket_to_pair_map, dot_bracket_to_stacks

def parse_pseudoknots(input_string, default_energy : float = -9.0, default_energy_tolerance : float = 1.0):
    
    # Split the input string into pseudoknot entries
    all_pk_dicts = input_string.split(';')
    final_pk_dicts = dict()
    
    for pk_idct in all_pk_dicts:
        stripped = pk_idct.strip()
        if not stripped:
            continue
        attribute_list = pk_idct.split(',')
        attr_fixed = []
        # adjust the atrribute: you may have lists there
        for attr in attribute_list:
            if ":" in attr:
                attr_fixed.append([x.strip() for x in attr.split(":")])
            else:
                attr_fixed[-1][1] += "," + attr
        pk_id = None
        ind_fwd = None
        ind_rev = None
        pk_energy = None
        pk_en_tolerance = None
        for k, v in attr_fixed:
            try:
                if k == 'id':
                    pk_id = v
                elif k == 'ind_fwd':
                    ind_fwd = ast.literal_eval(v)
                elif k == 'ind_rev':
                    ind_rev = ast.literal_eval(v)
                elif k == 'E':
                    pk_energy = float(v)
                elif k == 'dE':
                    pk_en_tolerance = float(v)
            except Exception as e:
                print(f"Error in parsing pseudoknots with key {k} and value {v}. Error: {e}")
        if not pk_id:
            raise ValueError("Pseudoknot id is missing")
        if pk_energy is None:
            pk_energy = default_energy
        if pk_en_tolerance is None:
            pk_en_tolerance = default_energy_tolerance
        if not ind_fwd or not ind_rev:
            warnings.warn(f"Skipping pseudoknot with id {pk_id} due to missing indices", stacklevel=3)
            continue
        try:
            all_ind = ind_fwd + ind_rev
            pk_seq_len = all_ind[0][1] - all_ind[0][0]
            for ind in all_ind:
                if ind[1] - ind[0] != pk_seq_len:
                    raise ValueError(f"Invalid indices for pseudoknot {pk_id}. All indices should have the same sequence length")
        except Exception as e:
            warnings.warn(f"Skipping pseudoknot with id {pk_id} due to invalid indices. Exception: {e}", stacklevel=3)
            continue
            
        final_pk_dicts[pk_id] = {'ind_fwd': ind_fwd, 'ind_rev': ind_rev, 'E': pk_energy, 'dE': pk_en_tolerance}

    return final_pk_dicts

def add_untracked_pseudoknots(pk_dict, structure, energy = -9.0, energy_tolerance = 1.0, pair_map = None):
    """ Add pseudoknots presents in the structure but not in the pk_dict to the pk_dict """
    if pair_map is None:
        pair_map = dot_bracket_to_pair_map(structure)

    reduced_db, stacks = dot_bracket_to_stacks(structure, only_opening = True)
    extra_pk_ind = 0
    all_pk_indexes = {x for pk_info in pk_dict.values() for x in pk_info['ind_fwd'] + pk_info['ind_rev']}
    for sym, indexes in zip(reduced_db, stacks):
        if sym not in '.(' and sym in db_pairs:

            # # check if the sequence is not constrained
            # if all(seq_constraints[i] in 'AUCG' for i in range(indexes[0], indexes[1] + 1)):
            #     continue

            if indexes not in all_pk_indexes:
                paired_indexes = [pair_map[x] for x in indexes[::-1]]
                pk_id = f"extra_{extra_pk_ind}"
                pk_dict[pk_id] = {  'id': pk_id,
                                    'ind_fwd': [indexes], 
                                    'ind_rev': [paired_indexes], 
                                    'E': energy, 
                                    'dE': energy_tolerance}
                extra_pk_ind += 1

    return pk_dict
