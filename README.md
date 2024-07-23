# cgen2charmm

##### Table of Contents  
1. [Introduction](#headers)
2. [Flags and Inputs](#headers)
3. [Examples](#headers)
4. [Contributing](#headers)

# Introduction
This is a commandline tool written in Python that makes the process of parameterizing multiple small molecules in one CHARMM forcefield directory much easier. This is meant to work with outputs from the [CGenFF website](https://cgenff.com/) and [CHARMM forcefields for GROMACS simulations](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5199616/). 
When generating and storing forcefield parameters in `ffbonded.itp` for multiple similar small molecules, duplicate entries of parameterized interactions can cause issues molecular dynamics production runs. This tool reads raw outputs from [CgenFF](https://cgenff.com/) and existing forcefield files, searches for redundant entries, converts from kcal to kJ when needed, and formats unique entries to be copy pasted into parameter files. 


# Flags and Inputs
# Examples 
# Contributing



