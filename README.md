# Generate_Subs
## A protocol for autogeneration of input files for investigation of the impact of substitution on molecular properties

This gitlab contains the necessary scripts to run an automated protocol for the generation of input files for quantum chemical calculations at the DFT level for substituted molecules, as outlined in [ARXIV LINK]. 

This protocol takes as input a set of key geometries (minima, transition states) for the molecule of interest, along with an input file defining relevant parameters for both the protol itself as well as the quantum chemical calculations to setup. It then generates input files for substituted versions of the molecule, through an automated attaching procedure. This allows calcualation of substitutued geometries to start from pre-optimised best guess structures, while also avoiding the need for human input and thus avoiding human error.

## Installation

The installation of this protocol is straightforward:

Download the files:

```git clone https://github.com/issymerritt/generate_subs```

(Optional) Create and activate a python environment

```conda create -n generate_subs```
```conda activate generate_subs```

Install required python packages (periodictable, numpy, pyscf)

```pip install -r requirements.txt```


## Usage

The main script used to run the program is the ```sub_generator.py```, which can be run with three options:

1. ```python sub_generator.py -S``` will list all substitutions available to the program in the substitution library
2. ```python sub_generator.py -C [INPUT_FILE]``` will check if the given inputfile is present and valid for the protocol
3. ```python sub_generator.py -G [INPUT_FILE]``` will run the protocol for the given inputfile, creating the requested input files for quantum chemical calculations.

## QC Compatability

The protocol is currently compatible with generation of input files for Gaussian16 and Orca6. Additional support for alternative QC programs can be simply added, through modification of the ```program_dependencies.py``` file. 

## Input file

This file defines the parameters of both how the automated protocol should be run, as well as the options to write into the generated QC input files. The file has a modular format, with 4 different block types. Examples of input files can be found in the example_usage/example_input folder.

### ```[GENERAL]``` block

This section defines parameters for the substitution protocol scripts.

Required arguments are:

- ```PATH_TO_CORE``` : informs the program where the input non-substitutued geometry xyz files can be located
- ```PROGRAM``` : which QC program the input files should be generated for (Gaussian16, Orca)

Optional arguments are:

- ```PARALLEL```: whether the input generator should be allowed to run in parallel or not (default = ON), the scripts will automatically determine how many cpus are available and the optimal way to parallelise.

### ```[PROG_PARAMS]``` block

This section defines universal parameters for the quantum chemistry input files to be written.

Required arguments are:

- ```FUNCTIONAL``` : the DFT functional (name convention should match that used by the QC program selected)
- ```BASIS_SET``` : the basis set (name convention should match that used by the QC program selected)
- ```CHARGE``` : the overall charge of the system
- ```SPIN``` : the total spin of the system (0 = singlet)

To select a solvent for PCM, either use;

- ```SOLVENT``` : the PCM solvent (for gas phase calculations, use SOLVENT = None)

or to define your own PCM solvent, use

- ```SOLV_EPS``` : dielectric constant
- ```SOLV_EPSINF``` : refractive index squared

Optional arguments for this section are:

- ```N_CPU``` : number of cpu for the QC calculations (default = 32)
- ```MEM``` : memory for the QC calculation (default = 63900MB)
- ```TDDFT``` : whether single point TD-DFT calculations should be run after optimisation of ground state geometries. Default = ON.

### ```[CORE_INFO x]``` block

Sequentially numbered sections (```[CORE_INFO 1], [CORE_INFO 2], ...```) are required for each "key" geometry (Minima, TS) to be investigated. Within each block, required arguments are:

- ```CORE_NAME``` : the name of the xyz file (without .xyz extension) containing the coordinates of this geometry for the unsubstitutued molecule
- ```STATE``` : whether this geometry corresponds to a stationary point on the ground state or an excited state (0 = GS, 1 = ES1, 2 = ES2, etc.)
- ```CORE_TYPE``` : whether substituted molecules starting from this geometry should be optimised to a minimum or a transition state

### ```[SUBSTITUTION x]``` block

Finally, sequentially numbered sections (```[SUBSTITUTION 1], [SUBSTITUTION 2], ...``` ) for each substitution position of the molecule are required. Within each block, required arguments are:

- ```SUBTYPES``` : a (comma separated) list of the folders within the substitution_library which contain the substitutions to ve investigated
- ```CORE_SUB_POS``` : the atom number (matching xyz ordering) of the substitution position on the central core geometry
- ```CORE_AT_TO_REM``` : the atom number (matching xyz ordering) of the hydrogen/other atom to be deleted from the central core geometry and replaced by the substitution
- ```FRAGMENT_LIST``` : a (comma separated) list of substitutions to be investigated. To include the option of no substitution at this position, the option NONE can be included in this list. 



