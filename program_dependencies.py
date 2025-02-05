from enum import Enum
from typing import Literal, Optional, Dict

from classes import Molecule


class QCprog(Enum):
    GAUSSIAN16 = 'GAUSSIAN16'
    ORCA = 'ORCA'


DEFAULT_CPU: int = 32
DEFAULT_MEMORY: int = 63900
SUPPORTED_PROGRAMS = list(QCprog)
AVAILABLE_FRAGS = ['EDG', 'EWG']  # TODO - fixed list ?


# TODO xyz extracvtor
# def extract_xyz_files(logname: str, script_location: str):
#     detect_program
#     if program == QCprog.GAUSSIAN16.value:
#
#     return


def generate_program_strings(coretype: Literal['MIN', 'TS'], coreState: int, program: str, memory: int = 63900,
                             cpus: int = 32, pcm_solvent: str = 'Generic, Read') -> Optional[Dict]:
    if program == QCprog.GAUSSIAN16.value:
        return get_gauss_strings(coretype, coreState, memory, cpus, pcm_solvent)
    elif program == QCprog.ORCA:
        return get_orca_strings(coretype, coreState, memory, cpus, pcm_solvent)  # TODO add other options ?


def write_qc_input_file(program: str, core: Molecule, program_strings: Dict, NAME: str, all_frag_string: str,
                        program_keyword_dict: Dict):
    if program == QCprog.GAUSSIAN16.value:
        write_gaussian_input(f"input_files/{NAME}/{NAME}_{all_frag_string}_{program_keyword_dict['STATE']}.com",
                             program_strings, program_keyword_dict, core)
    elif program == QCprog.ORCA.value:  # TODO make orca bit
        write_orca_input(f"input_files/{NAME}_{all_frag_string}_{program_keyword_dict['STATE']}.input", program_strings,
                         program_keyword_dict, core)


# Program specific functions
####### Gaussian 16 #######

def get_gauss_strings(Coretype: Literal['MIN', 'TS'], CoreState: int, memory: int = DEFAULT_MEMORY, cpus: int = DEFAULT_CPU,
                      pcm_solvent: str = 'Generic, Read') -> Dict:
    out_strings = {}
    if pcm_solvent != 'Generic, Read':
        out_strings['Solvent'] = None
    else:
        out_strings['Solvent'] = '\neps={SOLV_EPS}\nepsinf={SOLV_EPSINF}\n\n'
    if CoreState == 0:
        if Coretype == 'MIN':
            out_strings[
                'OptString'] = f'%mem={memory}MB\n%nprocshared={cpus}' + '\n%chk={NAME}.chk\n#P {FNAL}/{BAS_SET} ' + f'Opt Freq SCRF(PCM,Solvent={pcm_solvent})' + '\n \n Title \n \n {CHARGE} {SPIN}\n'
        elif Coretype == 'TS':
            out_strings[
                'OptString'] = f'%mem={memory}MB\n%nprocshared={cpus}' + '\n%chk={NAME}.chk\n#P {FNAL}/{BAS_SET} ' + f'Opt=(TS,NoEigen,calcfc) Freq SCRF(PCM,Solvent={pcm_solvent}) ' + '\n \n Title \n \n {CHARGE} {SPIN}\n'
        out_strings[
            'TDDFT_String'] = '\n--Link1--\n' + f'%mem={memory}MB\n%nprocshared={cpus}' + '\n%chk={NAME}.chk\n#P {FNAL}/{BAS_SET} Geom=AllCheck Guess=Read ' + f' SCRF(PCM,Solvent={pcm_solvent}) TD(NStates=3, Root=1) \n '
        return out_strings
    else:
        if Coretype == 'MIN':
            out_strings[
                'OptString'] = f'%mem={memory}MB\n%nprocshared={cpus}' + '\n%chk={NAME}.chk\n#P {FNAL}/{BAS_SET} TD(Nstates=3, Root={STATE}) ' + f'Opt Freq SCRF(PCM,Solvent={pcm_solvent})' + '\n \n Title \n \n {CHARGE} {SPIN}\n'
        elif Coretype == 'TS':
            out_strings[
                'OptString'] = f'%mem={memory}MB\n%nprocshared={cpus}' + '\n%chk={NAME}.chk\n#P {FNAL}/{BAS_SET} TD(Nstates=3, Root={STATE}) ' + f'Opt=(TS, noeigen, calcfc) Freq SCRF(PCM,Solvent={pcm_solvent})' + '\n \n Title \n \n {CHARGE} {SPIN}\n'
        return out_strings


def write_gaussian_input(filename: str, strings_gauss: dict, program_keyword_dict, geometry: Molecule):
    with open(filename, 'w') as gaussinp:
        gaussinp.write(strings_gauss['OptString'].format(**program_keyword_dict))
        gaussinp.write(geometry.gen_xyz_string().replace(';', '\n'))
        if strings_gauss['Solvent']:
            gaussinp.write(strings_gauss['Solvent'].format(**program_keyword_dict))
        if program_keyword_dict['STATE'] == 0 and program_keyword_dict.get('TDDFT', 'ON') != 'OFF':
            gaussinp.write(strings_gauss['TDDFT_String'].format(**program_keyword_dict))
            if strings_gauss['Solvent']:
                gaussinp.write(strings_gauss['Solvent'].format(**program_keyword_dict))
    return


####### ORCA #######

def get_orca_strings(Coretype: Literal['MIN', 'TS'], CoreState: int, memory: int = 63900, cpus: int = 32,
                     pcm_solvent: str = 'Water') -> Dict[str, str]:
    # TODO
    return {}


def write_orca_input(filename: str, strings_orca: dict, program_keyword_dict, geometry: Molecule):
    with open(filename, 'w') as orcainp:
        # TODO
        pass
    return
