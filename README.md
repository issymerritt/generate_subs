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

- ```PARALLEL```: whether the input generator should be allowed to run in parallel or not (default = ON), the scripts will automatically determine how many cpus are available.




