import configparser
import itertools
import os
from copy import deepcopy
from multiprocessing import Pool
from typing import List, Optional, Dict

import numpy as np
# noinspection PyUnresolvedReferences
from pyscf import gto, scf

from classes import Molecule
from key_fns import add_path_sub_library, parse_input_file, check_options_in_config, extract_program_keyword_dict, rotate_a_to_b, rotate_about_k, \
    print_to_file
from program_dependencies import AVAILABLE_FRAGS, SUPPORTED_PROGRAMS, generate_program_strings, write_qc_input_file


def generate_possible_combinations(config_info: configparser.ConfigParser) -> List:
    frag_lists = []
    nb_substitutions = sum(1 for s in config_info.sections() if 'SUBSTITUTION' in s)
    for sub_pos_idx in range(1, nb_substitutions + 1):
        core_position_sub = config_info[f'SUBSTITUTION {sub_pos_idx}']['CORE_SUB_POS']
        list_of_fragments = config_info[f'SUBSTITUTION {sub_pos_idx}']['FRAGMENT_LIST']
        path_to_fragment_library = config_info['GENERAL']['PATH_TO_FRAGMENT_LIBRARY']
        print('------------------------------------------------------------------------')
        print(f"Substitution position {sub_pos_idx} (Position {core_position_sub}): {list_of_fragments}")
        subs_at_curr_pos_idx = create_frag_list(path_to_fragment_library, parsed_options=config_info[f'SUBSTITUTION {sub_pos_idx}'])
        frag_lists.append(subs_at_curr_pos_idx)
    print('------------------------------------------------------------------------')
    print(frag_lists)

    list_of_combinations = list(itertools.product(*frag_lists))
    return list_of_combinations


def create_frag_list(path_to_fragment_library: str, parsed_options: Optional[configparser.SectionProxy] = None):
    fraglist = []
    for fragment_name in parsed_options['FRAGMENT_LIST'].replace(' ', '').split(','):
        try:
            fraglist.append(select_fragment(
                parsed_options['SUBTYPE'],
                fragment_name.upper(),
                int(parsed_options['CORE_AT_TO_REM']),
                int(parsed_options['CORE_SUB_POS']),
                path_to_fragment_library,
            ))
        except ValueError:
            exit('Check input file - substitution positions must be integers.')
    return fraglist


def select_fragment(fragtype: str, subst_name: str, h_rem: int, core_pos: int, frag_library_loc: str) -> Optional[Dict]:
    frag = {
        'FragType': fragtype,
        'H_rem': h_rem,
        'CorePos': core_pos,
    }
    poss_sub_list = [fname.upper().replace('.XYZ', '') for fname in os.listdir(f'{frag_library_loc}/{fragtype}')]
    if subst_name.upper() in poss_sub_list:
        frag['Name'] = subst_name
    else:
        raise Exception(f'Invalid input for type {fragtype} - must be one of {poss_sub_list}')
    return frag


def add_substitutions(inputfile: str, scriptloc: str):
    # Add PATH_TO_FRAGMENT_LIBRARY if needed
    add_path_sub_library(inputfile, scriptloc)
    # Read configuration info
    config_info = parse_input_file(inputfile, scriptloc)
    check_options_in_config(config_info, SUPPORTED_PROGRAMS, AVAILABLE_FRAGS, scriptloc)

    list_of_frag_comb = generate_possible_combinations(config_info)
    nb_core_geoms = sum(1 for s in config_info.sections() if 'CORE_INFO' in s)

    print(f'{nb_core_geoms} Core Geometries: [{", ".join(config_info[f"CORE_INFO {s}"]["CORE_NAME"] for s in range(1, nb_core_geoms + 1))}]')
    print('------------------------------------------------------------------------')

    if not os.path.isdir('xyz_files'):
        os.system(f'mkdir xyz_files')
    if not os.path.isdir('input_files'):
        os.system(f'mkdir input_files')

    # Run generation script for each individual core geometry
    for core_number in range(1, nb_core_geoms + 1):
        core_name = config_info[f'CORE_INFO {core_number}']['CORE_NAME']
        out_logfile = f'{core_name}_log.out'
        path_to_core = config_info['GENERAL']['PATH_TO_CORE']
        core_type = config_info[f'CORE_INFO {core_number}']['CORE_TYPE']
        core_state = config_info[f'CORE_INFO {core_number}']['STATE']

        program_keyword_dict = extract_program_keyword_dict(config_info, core_number)
        core = Molecule(filename=f"{path_to_core}/{core_name}")

        print(f"Core molecule: {core_name} - {len(core.coords)} atoms.  ")
        print_to_file(
            f"Core molecule: {core_name} - {len(core.coords)} atoms.  "
            f"Type = {core_type}, State = {core_state}", out_logfile
        )

        generate_molecules_for_one_core(core_number, config_info, list_of_frag_comb, program_keyword_dict, out_logfile)


def generate_molecules_for_one_core(core_number: int, config_info: configparser.ConfigParser, fragmentlist: List[Dict],
                                    program_keyword_dict: Dict, out_logfile: str):
    core_name = str(config_info[f'CORE_INFO {core_number}']['CORE_NAME'])
    path_to_frag_library = config_info['GENERAL']['PATH_TO_FRAGMENT_LIBRARY']
    path_to_core = config_info['GENERAL']['PATH_TO_CORE']
    program = config_info['GENERAL']['PROGRAM'].upper()
    print_to_file(f'Keywords used for CORE {core_name} : {program_keyword_dict}', out_logfile)

    program_strings = generate_program_strings(program_keyword_dict['CORETYPE'].upper(), program_keyword_dict['STATE'], program,
                                               program_keyword_dict['MEM'], program_keyword_dict['N_CPU'], program_keyword_dict['SOLVENT'])
    if not os.path.isdir(f'input_files/{core_name}'):
        os.system(f'mkdir input_files/{core_name}')

    print_to_file('Input files generated: ', out_logfile)
    parallel_process_attach_substitutions(fragmentlist, path_to_core, core_name, path_to_frag_library, out_logfile, program_keyword_dict, program,
                                          program_strings)


def parallel_process_attach_substitutions(fragmentlist: List[Dict], path_to_core: str, core_name: str, path_to_frag_library: str,
                                          out_logfile: str, program_keyword_dict: Dict, program: str, program_strings: Dict):
    with Pool() as pool:
        # Using starmap because we need to pass multiple arguments
        pool.starmap(attach_and_write_substitutions,
                     [(sublist_frags, path_to_core, core_name, path_to_frag_library, out_logfile, program_keyword_dict, program,
                       program_strings)
                      for sublist_frags in fragmentlist])


def attach_and_write_substitutions(sublist_frags: List[Dict], path_to_core: str, core_name: str, path_to_frag_library: str, out_logfile: str,
                                   program_keyword_dict: Dict, program: str, program_strings: Dict):
    # Sort such that the highest numbered H atom is removed first from core to avoid atom numbering conflicts
    sorted_fraglist = sorted(sublist_frags, key=lambda x: x['H_rem'], reverse=True)
    # reload core cleanly
    core = Molecule(filename=f'{path_to_core}/{core_name}')
    # Successively add each substitution for this particular combination
    for idx, fragment in enumerate(sorted_fraglist):
        frag_path = path_to_frag_library + f'/{fragment["FragType"]}/'
        substitution_name = fragment["Name"].upper()
        core = connect_geoms(core, substitution_name, frag_path, fragment["CorePos"], fragment["H_rem"])

    all_frag_string = '_'.join([fragment['Name'] for fragment in sorted_fraglist])
    generated_geo_name = f"{core_name}_{all_frag_string}_{program_keyword_dict['STATE']}"
    core.write_mol(outfile=f"xyz_files/{generated_geo_name}.xyz")
    program_keyword_dict['NAME'] = generated_geo_name
    write_qc_input_file(program, core, program_strings, core_name, all_frag_string, program_keyword_dict)
    print_to_file(f'{generated_geo_name}.com', out_logfile)


def connect_geoms(core: Molecule, fragment_name: str, path_to_fragment_library: str, core_sub_pos: int, core_at_to_remove: int):
    # print_to_file(f'Attaching {fragment_name} to position {core_sub_pos}', out_logfile)
    path_to_fragment = path_to_fragment_library + f'{fragment_name}'
    frag = Molecule(filename=path_to_fragment)

    ## Connect substituent to core
    core_atom = core_sub_pos
    core_h = core_at_to_remove
    fragment_atom = int(frag.label.split(';')[1].split()[-1])  # Fragment sub. position
    fragment_h = int(frag.label.split(';')[0].split()[-1])  # Fragment H to remove
    if len(frag.coords) > 2:
        adjacent_atom_fragment = int(frag.label.split(';')[2].split()[-1])  # Fragment adjacent atom
    else:
        adjacent_atom_fragment = False

    # Check if rotation of substitution will be necessary

    if len(frag.coords) == 2:  # Single atom substitution
        rotation_check = False
    elif len(frag.coords) == 3 and abs(
            frag.angle(fragment_h, fragment_atom, adjacent_atom_fragment)) > 179.5:  # Linear substitution (e.g. C=N group)
        rotation_check = False
    else:
        rotation_check = True

    ## Shift geometries such that core_atom (Core substitution position) and fragment_h (H to remove from Substitution) are superposed at the origin

    translate_core = core.coords[core_atom].xyz_coords
    translate_frag = frag.coords[fragment_h].xyz_coords

    for atom in core.coords:
        core.coords[atom].xyz_coords = core.coords[atom].xyz_coords - translate_core
    for atom in frag.coords:
        frag.coords[atom].xyz_coords = frag.coords[atom].xyz_coords - translate_frag

    ## Get unit vectors for C-H and H-F
    uvector_CH = core.coords[core_h].xyz_coords / np.linalg.norm(core.coords[core_h].xyz_coords)
    uvector_HF = frag.coords[fragment_atom].xyz_coords / np.linalg.norm(frag.coords[fragment_atom].xyz_coords)

    ## Get rotation matrix
    vec_F1_rot_norm, rot_matrix = rotate_a_to_b(uvector_HF, uvector_CH)

    ## Rotate fragment to align core_atom-fragment_h and fragment_atom-fragment_h
    for atom in frag.coords:
        if atom == fragment_atom:
            frag.coords[atom].xyz_coords = vec_F1_rot_norm * np.linalg.norm(frag.coords[fragment_atom].xyz_coords)
        else:
            frag.coords[atom].xyz_coords = np.matmul(rot_matrix, frag.coords[atom].xyz_coords)

    ## Remove excess hydrogens
    core_no_h = deepcopy(core)
    del core_no_h.coords[core_h]
    del frag.coords[fragment_h]
    core_no_h.coords[core_atom].set_label('CORE_SUB_POS')
    frag.coords[fragment_atom].set_label('FRAG_SUB_POS')
    if adjacent_atom_fragment:
        frag.coords[adjacent_atom_fragment].set_label('adjacent_atom_fragment')

    # Single atom substitution - no rotation test
    if not rotation_check:  # Single atom substitution
        combined_coord = [core_no_h.coords[idx] for idx in core_no_h.coords]
        combined_coord += [frag.coords[idx] for idx in frag.coords]
        joined_mol = Molecule(coord_list=combined_coord)
        # print_to_file('Skipping rotation check - linear or single atom substitution.', out_logfile)
        return joined_mol

    else:
        # print_to_file('Running relative substitution rotation check.', out_logfile)
        ## Generate fragment rotations
        frag_orients = [frag]
        for a in range(1, 12):
            rot_fragment = deepcopy(frag)
            theta = a * (np.pi / 6)
            fragment_rot_matrix = rotate_about_k(uvector_CH, theta)
            for atom in rot_fragment.coords:
                rot_fragment.coords[atom].xyz_coords = np.matmul(fragment_rot_matrix,
                                                                 rot_fragment.coords[atom].xyz_coords)
            frag_orients += [rot_fragment]

        ## Join fragments and cores and find lowest energy (HF, STO-3G) orientation

        combined_coords = []
        for fr in frag_orients:
            cc = [core_no_h.coords[idx] for idx in core_no_h.coords]
            cc += [fr.coords[idx] for idx in fr.coords]
            combined_coords += [cc]

        joined_molecules = [Molecule(coord_list=a) for a in combined_coords]

        lowest_e = 0
        selected_orientation = False
        for frg_or in joined_molecules:
            xyz_str = frg_or.gen_xyz_string()
            mol_scf_form = gto.M(atom=xyz_str, basis='STO-3G')
            mol_scf_form.verbose = 0
            HF_en = scf.RHF(mol_scf_form).kernel()
            if HF_en < lowest_e:
                lowest_e = HF_en
                selected_orientation = deepcopy(frg_or)

        if not selected_orientation: raise Exception(
            'Script failed to find lowest energy orientation of fragment molecule')

        return selected_orientation
