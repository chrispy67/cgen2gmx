# cgen2charmm

##### Table of Contents  
1. [Introduction](#headers)
2. [Flags and Inputs](#headers)
3. [Examples](#headers)
4. [Contributing](#headers)

# Introduction
The raw output from [CGenFF](https://cgenff.com/) is incompatible with [CHARMM forcefields for GROMACS simulations](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5199616/) in many ways, primarily units and force constants being in different columns. Organizing and curating these parameters by hand can be very tedious and error-prone.<br />

This is a commandline tool written in Python that makes the process of parameterizing multiple small molecules in one CHARMM forcefield directory much easier. This is meant to work with raw output files from CgenFF and existing `ffbonded.itp` files from CHARMM forcefield direcotries. 
When generating and storing forcefield parameters for several small molecules, duplicate entries of parameterized interactions can cause issues with molecular dynamics production runs. This tool reads raw outputs from [CGenFF](https://cgenff.com/) and existing forcefield files, searches for redundant (duplicate) entries  --bonds, angles, dihedrals, and improper dihedrals-- converts from kcal to kJ when needed, and formats unique entries to be directly copy/pasted into a forcefield file

# Flags 
 #### `--cgen`:
 **Path to raw output from CGenFF output file to be read and parsed** 
 #### `--itp`:
 **Path to existing `ffbonded.itp` to be read and parsed.** 
 <br />
NOTE: There are unique functions for reading both forcefield files and CGenFF outputs. There is some error handling built in to account for missing improper dihedral parameters (common) and small formatting changes, but A POORLY ORGANIZED `ffbonded.itp` FILE WILL CAUSE ISSUES! The output format of this script will fit in well with standard CHARMM forcefields. 
 #### `--output`:
 **Path to desired output file** 
  <br />
Output file format closely matches standard CHARMM forcefields. If there is already a file of that name, the option will be given to continue by overwriting or exiting. Header columns are always written for each parameter in `ffbonded.itp`, regardless if any unique parameters are printed out. Units will be indicated where necessary in brackets. 
 #### `--kJ`:
# Examples 
# Contributing



