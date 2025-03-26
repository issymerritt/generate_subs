# Generate_Subs
## A protocol for autogeneration of input files for investigation of the impact of substitution on molecular properties

This gitlab contains the necessary scripts to run an automated protocol for the generation of input files for quantum chemical calculations at the DFT level for substituted molecules, as outlined in [ARXIV LINK]. 

This protocol takes as input a set of key geometries (minima, transition states) for the molecule of interest, along with an input file defining relevant parameters for both the protol itself as well as the quantum chemical calculations to setup. It then generates input files for substituted versions of the molecule, through an automated attaching procedure. This allows calcualation of substitutued geometries to start from pre-optimised best guess structures, while also avoiding the need for human input and thus avoiding human error.

### Installation

The installation of this protocol is straightforward:

Download the files:

```git clone https://github.com/issymerritt/generate_subs```

(Optional) Create and activate a python environment

```conda create -n generate_subs```
```conda activate generate_subs```

Install required python packages (periodictable, numpy, pyscf)

```pip install -r requirements.txt```


### Usage

The main script used to run the program is the ```sub_generator.py```, which can be run with three options:

1. ```python sub_generator.py -S``` will list all substitutions available to the program in the substitution library
2. ```python sub_generator.py -C [INPUT_FILE]``` will check if the given inputfile is present and valid for the protocol
3. ```python sub_generator.py -G [INPUT_FILE]``` will run the protocol for the given inputfile, creating the requested input files for quantum chemical calculations.

### QC Compatability

The protocol is currently compatible with generation of input files for Gaussian16 and Orca6. Additional support for alternative QC programs can be simply added, through modification of the ```program_dependencies.py``` file. 

### Input file

This file defines the 
