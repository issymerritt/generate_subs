import math
from enum import Enum
from typing import Literal, Optional, Dict

from classes import Molecule


class QCprog(Enum):
    GAUSSIAN16 = 'GAUSSIAN16'
    ORCA = 'ORCA'


def round_down_to_1000(n):
    return math.floor(n / 1000) * 1000


DEFAULT_CPU: int = 32
DEFAULT_MEMORY: int = 63900
SUPPORTED_PROGRAMS = list(QCprog)
AVAILABLE_FRAGS = ['EDG', 'EWG']  # TODO - fixed list ?


def generate_program_strings(coretype: Literal['MIN', 'TS'], coreState: int, program: str, memory: int = 63900,
                             cpus: int = 32, pcm_solvent: str = 'Generic, Read') -> Optional[Dict]:
    if program == QCprog.GAUSSIAN16.value:
        return get_gauss_strings(coretype, coreState, memory, cpus, pcm_solvent)
    elif program == QCprog.ORCA.value:
        return get_orca_strings(coretype, coreState, memory, cpus, pcm_solvent)  # TODO add other options ?


def write_qc_input_file(program: str, core: Molecule, program_strings: Dict, NAME: str, all_frag_string: str,
                        program_keyword_dict: Dict):
    if program == QCprog.GAUSSIAN16.value:
        write_gaussian_input(f"input_files/{NAME}/{NAME}_{all_frag_string}_{program_keyword_dict['STATE']}.com",
                             program_strings, program_keyword_dict, core)
    elif program == QCprog.ORCA.value:  # TODO make orca bit
        write_orca_input(f"input_files/{NAME}/{NAME}_{all_frag_string}_{program_keyword_dict['STATE']}.inp", program_strings,
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
            if pcm_solvent.upper() == 'NONE':
                out_strings[
                    'OptString'] = f'%mem={memory}MB\n%nprocshared={cpus}' + '\n%chk={NAME}.chk\n#P {FNAL}/{BAS_SET} ' + f'Opt Freq' + '\n \n Title \n \n {CHARGE} {SPIN}\n'
            else:
                out_strings[
                    'OptString'] = f'%mem={memory}MB\n%nprocshared={cpus}' + '\n%chk={NAME}.chk\n#P {FNAL}/{BAS_SET} ' + f'Opt Freq SCRF(PCM,Solvent={pcm_solvent})' + '\n \n Title \n \n {CHARGE} {SPIN}\n'
        elif Coretype == 'TS':
            if pcm_solvent.upper() == 'NONE':
                out_strings[
                    'OptString'] = f'%mem={memory}MB\n%nprocshared={cpus}' + '\n%chk={NAME}.chk\n#P {FNAL}/{BAS_SET} ' + f'Opt=(TS,NoEigen,calcfc) Freq' + '\n \n Title \n \n {CHARGE} {SPIN}\n'
            else:
                out_strings[
                    'OptString'] = f'%mem={memory}MB\n%nprocshared={cpus}' + '\n%chk={NAME}.chk\n#P {FNAL}/{BAS_SET} ' + f'Opt=(TS,NoEigen,calcfc) Freq SCRF(PCM,Solvent={pcm_solvent}) ' + '\n \n Title \n \n {CHARGE} {SPIN}\n'
        if pcm_solvent.upper() == 'NONE':
            out_strings[
                'TDDFT_String'] = '\n--Link1--\n' + f'%mem={memory}MB\n%nprocshared={cpus}' + '\n%chk={NAME}.chk\n#P {FNAL}/{BAS_SET} Geom=AllCheck Guess=Read ' + f' TD(NStates=3, Root=1) \n '
        else:
            out_strings[
                'TDDFT_String'] = '\n--Link1--\n' + f'%mem={memory}MB\n%nprocshared={cpus}' + '\n%chk={NAME}.chk\n#P {FNAL}/{BAS_SET} Geom=AllCheck Guess=Read ' + f' SCRF(PCM,Solvent={pcm_solvent}) TD(NStates=3, Root=1) \n '
        return out_strings
    else:
        if Coretype == 'MIN':
            if pcm_solvent.upper() == 'NONE':
                out_strings[
                    'OptString'] = f'%mem={memory}MB\n%nprocshared={cpus}' + '\n%chk={NAME}.chk\n#P {FNAL}/{BAS_SET} TD(Nstates=3, Root={STATE}) ' + f'Opt Freq' + '\n \n Title \n \n {CHARGE} {SPIN}\n'
            else:
                out_strings[
                    'OptString'] = f'%mem={memory}MB\n%nprocshared={cpus}' + '\n%chk={NAME}.chk\n#P {FNAL}/{BAS_SET} TD(Nstates=3, Root={STATE}) ' + f'Opt Freq SCRF(PCM,Solvent={pcm_solvent})' + '\n \n Title \n \n {CHARGE} {SPIN}\n'
        elif Coretype == 'TS':
            if pcm_solvent.upper() == 'NONE':
                out_strings[
                    'OptString'] = f'%mem={memory}MB\n%nprocshared={cpus}' + '\n%chk={NAME}.chk\n#P {FNAL}/{BAS_SET} TD(Nstates=3, Root={STATE}) ' + f'Opt=(TS, noeigen, calcfc) Freq' + '\n \n Title \n \n {CHARGE} {SPIN}\n'
            else:
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
    out_strings = {}
    if CoreState == 0:  # GS Geometry optimisation (DFT) followed by TD-DFT single point calculation (if requested)
        if Coretype == 'MIN':  # Minimum
            if pcm_solvent.upper() == 'NONE':
                out_strings[
                    'OptString'] = '! {FNAL} {BAS_SET} OPT FREQ' + f'\n\n%maxcore {round_down_to_1000(memory / cpus)} \n%PAL NPROCS {cpus} END \n' + '\n* xyz {CHARGE} {SPIN}\n'
            elif pcm_solvent == 'Generic, Read':
                out_strings[
                    'OptString'] = '! {FNAL} {BAS_SET} OPT FREQ CPCM \n\n%CPCM \n EPSILON    {SOLV_EPS}\n REFRAC    {SOLV_EPSINF}\nEND\n' + f'\n%maxcore {round_down_to_1000(memory / cpus)} \n%PAL NPROCS {cpus} END \n' + '\n* xyz {CHARGE} {SPIN}\n'
            else:
                out_strings[
                    'OptString'] = '! {FNAL} {BAS_SET} OPT FREQ' + f' CPCM({pcm_solvent}) \n\n%maxcore {round_down_to_1000(memory / cpus)} \n%PAL NPROCS {cpus} END \n' + '\n* xyz {CHARGE} {SPIN}\n'

        elif Coretype == 'TS':  # Transition state
            if pcm_solvent.upper() == 'NONE':
                out_strings[
                    'OptString'] = '! {FNAL} {BAS_SET} OptTS FREQ' + f'\n\n%maxcore {round_down_to_1000(memory / cpus)} \n%PAL NPROCS {cpus} END \n' + '\n* xyz {CHARGE} {SPIN}\n'
            elif pcm_solvent == 'Generic, Read':
                out_strings[
                    'OptString'] = '! {FNAL} {BAS_SET} OptTS FREQ CPCM \n\n%CPCM \n EPSILON    {SOLV_EPS}\n REFRAC    {SOLV_EPSINF}\nEND\n' + f'\n%maxcore {round_down_to_1000(memory / cpus)} \n%PAL NPROCS {cpus} END \n' + '\n* xyz {CHARGE} {SPIN}\n'
            else:
                out_strings[
                    'OptString'] = '! {FNAL} {BAS_SET} OptTS FREQ' + f' CPCM({pcm_solvent}) \n\n%maxcore {round_down_to_1000(memory / cpus)} \n%PAL NPROCS {cpus} END\n' + '\n* xyz {CHARGE} {SPIN}\n'

        if pcm_solvent.upper() == 'NONE':
            out_strings[
                'TDDFT_String'] = '! {FNAL} {BAS_SET} \n \n %TDDFT \n NROOTS 5 \n END' + f'\n\n%maxcore {round_down_to_1000(memory / cpus)} \n%PAL NPROCS {cpus} END\n '
        elif pcm_solvent == 'Generic, Read':
            out_strings[
                'TDDFT_String'] = '! {FNAL} {BAS_SET} CPCM \n\n%CPCM \n EPSILON    {SOLV_EPS}\n REFRAC    {SOLV_EPSINF}\nEND' + ' \n\n %TDDFT NROOTS 5 END' + f'\n\n%maxcore {round_down_to_1000(memory / cpus)} \n%PAL NPROCS {cpus} END\n '
        else:
            out_strings[
                'TDDFT_String'] = '! {FNAL} {BAS_SET} ' + f' CPCM({pcm_solvent}) ' + ' \n \n %TDDFT \nNROOTS 5 \nEND' + f'\n\n%maxcore {round_down_to_1000(memory / cpus)} \n%PAL NPROCS {cpus} END\n '
        return out_strings

    else:  # ES Geometry optimisation (TD-DFT)
        if Coretype == 'MIN':
            if pcm_solvent.upper() == 'NONE':
                out_strings[
                    'OptString'] = '! {FNAL} {BAS_SET} OPT FREQ\n\n%TDDFT\nNROOTS 5 \n IRoot {STATE}\n END \n' + f'\n%maxcore {round_down_to_1000(memory / cpus)} \n%PAL NPROCS {cpus} END\n ' + '\n* xyz {CHARGE} {SPIN}\n'
            elif pcm_solvent == 'Generic, Read':
                out_strings[
                    'OptString'] = '! {FNAL} {BAS_SET} OPT FREQ CPCM \n\n%TDDFT\nNROOTS 5 \n IRoot {STATE}\n END  \n\n%CPCM \n EPSILON    {SOLV_EPS}\n REFRAC    {SOLV_EPSINF}\nEND' + f'\n\n%maxcore {round_down_to_1000(memory / cpus)} \n%PAL NPROCS {cpus} END \n' + '\n* xyz {CHARGE} {SPIN}\n'
            else:
                out_strings[
                    'OptString'] = '! {FNAL} {BAS_SET} OPT FREQ' + f' CPCM({pcm_solvent})' + '\n%TDDFT\nNROOTS 5 \n IRoot {STATE}\n END \n' + f'\n%maxcore {round_down_to_1000(memory / cpus)} \n%PAL NPROCS {cpus} END \n' + '\n* xyz {CHARGE} {SPIN}\n'
        elif Coretype == 'TS':
            if pcm_solvent.upper() == 'NONE':
                out_strings[
                    'OptString'] = '! {FNAL} {BAS_SET} OptTS FREQ\n\n%TDDFT\nNROOTS 5 \n IRoot {STATE}\n END \n' + f'\n\n%maxcore {round_down_to_1000(memory / cpus)} \n%PAL NPROCS {cpus} END \n' + '\n* xyz {CHARGE} {SPIN}\n'
            elif pcm_solvent == 'Generic, Read':
                out_strings[
                    'OptString'] = '! {FNAL} {BAS_SET} OptTS FREQ CPCM \n\n%TDDFT\nNROOTS 5 \n IRoot {STATE}\n END  \n\n%CPCM \n EPSILON    {SOLV_EPS}\n REFRAC    {SOLV_EPSINF}\nEND' + f'\n\n%maxcore {round_down_to_1000(memory / cpus)} \n%PAL NPROCS {cpus} END \n' + '\n* xyz {CHARGE} {SPIN}\n'
            else:
                out_strings[
                    'OptString'] = '! {FNAL} {BAS_SET} OptTS FREQ' + f' CPCM({pcm_solvent})' + '\n\n%TDDFT\nNROOTS 5 \n IRoot {STATE}\n END \n' + f'\n\n%maxcore {round_down_to_1000(memory / cpus)} \n%PAL NPROCS {cpus} END \n' + '\n* xyz {CHARGE} {SPIN}\n'
    return out_strings


def write_orca_input(filename: str, strings_orca: dict, program_keyword_dict, geometry: Molecule):
    with open(filename, 'w') as orcainp:
        orcainp.write(strings_orca['OptString'].format(**program_keyword_dict))
        orcainp.write(geometry.gen_xyz_string().replace(';', '\n'))
        orcainp.write('*\n\n')
        if program_keyword_dict['STATE'] == 0 and program_keyword_dict.get('TDDFT', 'ON') != 'OFF':
            orcainp.write(strings_orca['TDDFT_String'].format(**program_keyword_dict))
    return
