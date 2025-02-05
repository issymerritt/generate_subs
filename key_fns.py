import configparser
import os
from typing import List, Optional, Dict, Tuple, Union

import numpy as np
from numpy import ndarray

from program_dependencies import SUPPORTED_PROGRAMS, AVAILABLE_FRAGS, DEFAULT_CPU, DEFAULT_MEMORY


def run_input_checks(inputfile: str, scriptloc: str):
    config_info = parse_input_file(inputfile, scriptloc)
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


def parse_input_file(file_path: str, script_location: str) -> Optional[configparser.ConfigParser]:
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


def check_options_in_config(config_info: configparser.ConfigParser, supported_programs: List, available_frags: List, scriptloc: str):
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
            raise Exception(f'Core {core_name} not found at location {path_to_core_library}.\n Available cores: {available_cores} (care sensitive)')

    for sub_nb in range(1, nb_sub_positions + 1):
        if config_info[f'SUBSTITUTION {sub_nb}']['SUBTYPE'].upper() not in available_frags:
            raise Exception(f'Unrecognised substitution type selected - available substitution types are : ', available_frags)
        try:
            int(config_info[f'SUBSTITUTION {sub_nb}']['CORE_AT_TO_REM'])
            int(config_info[f'SUBSTITUTION {sub_nb}']['CORE_SUB_POS'])
        except ValueError:
            raise Exception(f'Substitution position variables (SUBSTITUTION {sub_nb}) must be integers')

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

    path_to_frag_library = config_info.get('GENERAL', 'PATH_TO_FRAGMENT_LIBRARY', fallback=False)
    if path_to_frag_library and path_to_frag_library != f'{scriptloc}/substituent_library/':
        print('WARNING : PATH_TO_FRAGMENT_LIBRARY defined in input file as a different directory to standard substitution library. \n'
              '          Remove this line unless using personal substitution library.')


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
