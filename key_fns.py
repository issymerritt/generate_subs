import configparser
import itertools
import os
from copy import deepcopy
from multiprocessing import Pool, cpu_count
from typing import List, Optional, Dict, Tuple, Union

import numpy as np
from numpy import ndarray
# noinspection PyUnresolvedReferences
from pyscf import gto, scf
from tqdm import tqdm

from classes import Molecule
from program_dependencies import SUPPORTED_PROGRAMS, AVAILABLE_FRAGS, DEFAULT_CPU, DEFAULT_MEMORY, generate_program_strings, write_qc_input_file


######################### GENERIC KEY FUNCTIONS #########################

def set_script_location():
    return os.path.realpath(os.path.dirname(__file__))


def run_input_checks(inputfile: str, scriptloc: str):
    config_info = parse_input_file(inputfile)
    check_options_in_config(config_info, SUPPORTED_PROGRAMS, AVAILABLE_FRAGS, scriptloc)
    print('Input file is valid.')


def list_substitution_library(script_location: str):
    substitution_library_path = os.path.join(script_location, 'substituent_library')
    if not os.path.exists(substitution_library_path):
        raise Exception(f"Substitution library directory not found: {substitution_library_path}")
    print("Available substitutions:")
    for subfolder in os.listdir(substitution_library_path):
        subfolder_path = os.path.join(substitution_library_path, subfolder)
        if os.path.isdir(subfolder_path):
            # Get files in the subfolder
            files = [a.split('.')[0] for a in os.listdir(subfolder_path)]
            # Print subfolder and its files
            print(f"[{subfolder}]: {', '.join(files)}")


def parse_input_file(file_path: str) -> Optional[configparser.ConfigParser]:
    if not os.path.isfile(file_path):
        raise Exception(f'Input file {file_path} not found. ')

    config = configparser.ConfigParser()
    config.read(file_path)

    return config


def add_path_sub_library(file_path: str, script_location: str):
    """ Add the path to the substitution library (where scripts are saved) to the input file, if no path is given in the input file."""
    with open(file_path) as input_file:
        input_contents = input_file.readlines()

    frag_library_line = f'PATH_TO_FRAGMENT_LIBRARY = {script_location}/substituent_library/\n'
    frag_path_exists = any('PATH_TO_FRAGMENT_LIBRARY' in line for line in input_contents)

    if not frag_path_exists:
        input_contents.insert(2, frag_library_line)
        with open(file_path, 'w') as input_file:
            input_file.writelines(input_contents)


def check_options_in_config(config_info: configparser.ConfigParser, supported_programs: List, available_frag_type: List, scriptloc: str):
    # Check if all required options exist
    req_option_list_gen = ['PATH_TO_CORE', 'PROGRAM']
    nb_core_geoms = sum(1 for s in config_info.sections() if 'CORE_INFO' in s)
    nb_sub_positions = sum(1 for s in config_info.sections() if 'SUBSTITUTION' in s)
    req_option_list_core = ['CORE_NAME', 'STATE', 'CORE_TYPE']
    req_option_list_qc = ['BASIS SET', 'FUNCTIONAL', 'CHARGE', 'SPIN']
    if any(option is False for option in [config_info.get('GENERAL', a) for a in req_option_list_gen]):
        raise Exception('Required inputs in GENERAL section: ', req_option_list_gen)
    elif any(option is False for option in [config_info.get('PROG_PARAMS', a) for a in req_option_list_qc]):
        raise Exception('Required inputs in PROG_PARAMS section: ', req_option_list_qc)
    else:
        for core_nb in range(1, nb_core_geoms + 1):
            if any(option is False for option in [config_info.get(f'CORE_INFO {core_nb}', a) for a in req_option_list_core]):
                raise Exception(f'Required inputs in CORE_INFO {core_nb} section: ', req_option_list_core)

    # Logical Checks
    path_to_core_library = config_info.get('GENERAL', 'PATH_TO_CORE', fallback=False)
    available_cores = [os.path.splitext(file)[0] for file in os.listdir(path_to_core_library)]
    for core_nb in range(1, nb_core_geoms + 1):
        core_name = config_info[f'CORE_INFO {core_nb}']['CORE_NAME']
        core_type = config_info[f'CORE_INFO {core_nb}']['CORE_TYPE'].upper()
        core_state = config_info[f'CORE_INFO {core_nb}']['STATE']
        if core_type not in ['MIN', 'TS']:
            raise Exception(f'Core type (CORE_INFO {core_nb}) must be either MIN or TS')
        try:
            int(core_state)
        except ValueError:
            raise Exception(f'Core state (CORE_INFO {core_nb}) must be an integer')
        if core_name not in available_cores:
            raise Exception(f'Core {core_name} not found at location {path_to_core_library}.\n Available cores: {available_cores} (cae sensitive)')

    path_to_frag_library = config_info.get('GENERAL', 'PATH_TO_FRAGMENT_LIBRARY', fallback=False)
    if path_to_frag_library and path_to_frag_library != f'{scriptloc}/substituent_library/':
        print('WARNING : PATH_TO_FRAGMENT_LIBRARY defined in input file as a different directory to standard substitution library. \n'
              '          Remove this line unless using personal substitution library.')
    elif not path_to_frag_library:
        path_to_frag_library = f'{scriptloc}/substituent_library/'

    for sub_nb in range(1, nb_sub_positions + 1):
        substitution_type = config_info[f'SUBSTITUTION {sub_nb}']['SUBTYPE'].upper()
        if substitution_type not in available_frag_type:
            raise Exception(f'Unrecognised substitution type selected - available substitution types are : ', available_frag_type)
        try:
            int(config_info[f'SUBSTITUTION {sub_nb}']['CORE_AT_TO_REM'])
            int(config_info[f'SUBSTITUTION {sub_nb}']['CORE_SUB_POS'])
        except ValueError:
            raise Exception(f'Substitution position variables (SUBSTITUTION {sub_nb}) must be integers')
        available_substitutions = [a.split('.')[0] for a in os.listdir(f"{path_to_frag_library}/{substitution_type}")]
        for sub_name in config_info[f'SUBSTITUTION {sub_nb}']['FRAGMENT_LIST'].replace(',', ' ').split():
            if sub_name.upper() not in available_substitutions:
                raise Exception(f"Substitution {sub_name} not found at location {path_to_frag_library}/{substitution_type}/{sub_name}.")

    if config_info.get('PROG_PARAMS', 'SOLV_EPS', fallback=False) and config_info.get('PROG_PARAMS', 'SOLVENT', fallback=False):
        raise Exception('Either solvent parameters or a solvent name must be given in input file, not both.')
    elif not config_info.get('PROG_PARAMS', 'SOLV_EPS', fallback=False) and not config_info.get('PROG_PARAMS', 'SOLVENT', fallback=False):
        raise Exception('One of either solvent parameters or a solvent name must be given in input file.')
    elif not config_info.get('PROG_PARAMS', 'SOLV_EPS', fallback=False) or not config_info.get('PROG_PARAMS', 'SOLV_EPSINF', fallback=False):
        if any(val is True for val in [
            config_info.get('PROG_PARAMS', 'SOLV_EPS', fallback=False),
            config_info.get('PROG_PARAMS', 'SOLV_EPSINF', fallback=False),
        ]):
            raise Exception('Both SOLV_EPS and SOLV_EPSINF must be given when defining own solvent parameters.')
    elif config_info['GENERAL']['PROGRAM'].upper() not in supported_programs:
        raise Exception('PROGRAM must be one of : ', supported_programs)
    if config_info.get('PROG_PARAMS', 'TDDFT', fallback='ON') not in ['OFF', 'ON']:
        raise Exception('TDDFT keyword must be followed by either OFF or ON (default = ON for GS calculations)')

    try:
        int(config_info['PROG_PARAMS']['CHARGE'])
        int(config_info['PROG_PARAMS']['SPIN'])
    except:
        raise Exception('Charge and spin must be integers.')

    if not config_info.get('PROG_PARAMS', 'SOLVENT', fallback=False):
        try:
            float(config_info['PROG_PARAMS']['SOLV_EPS'])
            float(config_info['PROG_PARAMS']['SOLV_EPSINF'])
        except ValueError:
            raise ValueError('Solvent parameters (SOLV_EPS and SOLV_EPSINF) must be floats')

    if config_info.get('PROG_PARAMS', 'N_CPU', fallback=False):
        try:
            int(config_info['PROG_PARAMS']['N_CPU'])
        except ValueError:
            raise ValueError('N_CPU must be an integer value.')

    if config_info.get('PROG_PARAMS', 'MEM', fallback=False):
        try:
            int(config_info['PROG_PARAMS']['MEM'])
        except ValueError:
            raise ValueError('MEM must be an integer value.')


def extract_program_keyword_dict(config_info: configparser.ConfigParser, core_number: int) -> Dict:
    program_keyword_dict = dict(FNAL=config_info['PROG_PARAMS']['FUNCTIONAL'], BAS_SET=config_info['PROG_PARAMS']['BASIS SET'],
                                CHARGE=int(config_info['PROG_PARAMS']['CHARGE']), SPIN=int(config_info['PROG_PARAMS']['SPIN']),
                                STATE=int(config_info[f'CORE_INFO {core_number}']['STATE']),
                                CORETYPE=config_info[f'CORE_INFO {core_number}']['CORE_TYPE'])

    program_keyword_dict['TDDFT'] = config_info.get('PROG_PARAMS', 'TDDFT', fallback='ON')
    solvent = config_info.get('PROG_PARAMS', 'SOLVENT', fallback=False)
    if solvent:
        program_keyword_dict['SOLVENT'] = solvent
    else:
        program_keyword_dict['SOLV_EPS'] = float(config_info['PROG_PARAMS']['SOLV_EPS'])
        program_keyword_dict['SOLV_EPSINF'] = float(config_info['PROG_PARAMS']['SOLV_EPSINF'])
        program_keyword_dict['SOLVENT'] = 'Generic, Read'

    # Read optional parameters (mem/cpu)
    if config_info.get('PROG_PARAMS', 'N_CPU', fallback=False):
        program_keyword_dict['N_CPU'] = int(config_info['PROG_PARAMS']['N_CPU'])
    else:
        program_keyword_dict['N_CPU'] = DEFAULT_CPU

    if config_info.get('PROG_PARAMS', 'MEM', fallback=False):
        program_keyword_dict['MEM'] = int(config_info['PROG_PARAMS']['MEM'])
    else:
        program_keyword_dict['MEM'] = DEFAULT_MEMORY

    return program_keyword_dict


def rotate_a_to_b(A: np.array, B: np.array) -> Tuple[np.array, np.array]:
    """ Rotates vector A onto a vector B """

    # Basis for rotation
    u = A
    v = (B - np.dot(A, B) * A) / np.linalg.norm(B - np.dot(A, B) * A)
    w = np.cross(B, A) / np.linalg.norm(np.cross(B, A))

    basis = np.column_stack((u, v, w))
    F = np.linalg.inv(basis)

    # Rotation matrix
    G = np.array([
        [np.dot(A, B), -np.linalg.norm(np.cross(A, B)), 0],
        [np.linalg.norm(np.cross(A, B)), np.dot(A, B), 0],
        [0, 0, 1]
    ])

    # Transformation matrix
    trans_matrix = np.matmul(basis, np.matmul(G, F))

    # Rotate A onto B
    rot_A = np.matmul(trans_matrix, A)

    return rot_A, trans_matrix


def rotate_about_k(k: ndarray, theta: Union[int, float]):
    """Returns rotation matrix (Rodriguez Formula) to rotate a vector by angle theta (in radians) about unit vector k"""
    if type(k) != np.ndarray or len(k) != 3: raise Exception('k must be of type ndarray, size 3')
    if type(theta) not in [int, float]: raise Exception('theta must either integer or floating point number')

    if np.linalg.norm(k) != 1:
        # print('Renormalising axis of rotation vector')
        k = k / np.linalg.norm(k)

    K_mat = np.array([
        [0, -k[2], k[1]],
        [k[2], 0, -k[0]],
        [-k[1], k[0], 0]
    ])

    r_mat = np.identity(3) + np.sin(theta) * K_mat + (1 - np.cos(theta)) * np.matmul(K_mat, K_mat)

    return r_mat


def print_to_file(writestring: str, filename: str):
    with open(filename, 'a') as outfile:
        outfile.writelines(writestring)
        outfile.writelines('\n')


################## KEY FUNCTIONS CONNECTING GEOMETRIES ##################


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
    for fragment_name in parsed_options['FRAGMENT_LIST'].replace(',', ' ').split():
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
        raise Exception(f'{subst_name.upper()} is invalid input for type {fragtype} - must be one of {poss_sub_list}')
    return frag


def add_substitutions(inputfile: str, scriptloc: str):
    # Add PATH_TO_FRAGMENT_LIBRARY if needed
    add_path_sub_library(inputfile, scriptloc)
    # Read configuration info
    config_info = parse_input_file(inputfile)
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
    total_cpus = cpu_count()
    cpus_per_process = min(8, total_cpus)  # Use up to 8 CPUs per process, but not more than available
    num_processes = max(1, total_cpus // cpus_per_process)  # Ensure at least 1 process runs
    with Pool(processes=num_processes, initializer=set_threads, initargs=(cpus_per_process,)) as pool:
        # Using starmap because we need to pass multiple arguments
        pool.starmap(attach_and_write_substitutions,
                     [(sublist_frags, path_to_core, core_name, path_to_frag_library, out_logfile, program_keyword_dict, program,
                       program_strings)
                      for sublist_frags in fragmentlist])


def set_threads(cpus_per_process):
    """Set environment variables to limit CPU usage per process."""
    os.environ["OMP_NUM_THREADS"] = str(cpus_per_process)  # OpenMP-based libraries
    os.environ["MKL_NUM_THREADS"] = str(cpus_per_process)  # Intel MKL
    os.environ["NUMEXPR_NUM_THREADS"] = str(cpus_per_process)  # NumExpr
    os.environ["OPENBLAS_NUM_THREADS"] = str(cpus_per_process)  # OpenBLAS
    os.environ["VECLIB_MAXIMUM_THREADS"] = str(cpus_per_process)  # macOS Accelerate


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
        for frg_or in tqdm(joined_molecules):
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
