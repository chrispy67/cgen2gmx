# cgen2charmm

##### Table of Contents  
1. [Introduction](https://github.com/chrispy67/cgen2gmx#introduction)
2. [Installation](https://github.com/chrispy67/cgen2gmx#installation)
3. [Flags and Inputs](https://github.com/chrispy67/cgen2gmx#flags-and-inputs)
4. [Examples](https://github.com/chrispy67/cgen2gmx#examples)
5. [Contributing](https://github.com/chrispy67/cgen2gmx#contributing)

# Introduction
The raw output from [CGenFF](https://cgenff.com/) is incompatible with [CHARMM forcefields for GROMACS simulations](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5199616/) in many ways, primarily units and force constants being in different columns. Organizing and curating these parameters by hand can be very tedious and error-prone.<br />

This is a commandline tool written in Python that makes the process of parameterizing multiple small molecules in one CHARMM forcefield directory much easier. This is meant to work with raw output files from CgenFF and existing `ffbonded.itp` files from CHARMM forcefield direcotries. 
When generating and storing forcefield parameters for several small molecules, duplicate entries of parameterized interactions can cause issues with molecular dynamics production runs. This tool reads raw outputs from [CGenFF](https://cgenff.com/) and existing forcefield files, searches for redundant (duplicate) entries--bonds, angles, dihedrals, and improper dihedrals--converts from kcal to kJ when needed, and formats unique entries to be directly copy/pasted into an existing forcefield file.

# Installation

There is a pip package installable via `pip install cgen2gmx` that will enable the use of the package anywhere by using `cgen2gmx` anywhere in the CLI. 

If you wish to make changes to the source code for your specific use case, cloning or branching the repo is reccomended. 


# Flags and Inputs 
 #### `--cgen`: **Path to raw output from CGenFF output file to be read and parsed** 
 #### `--itp`: **Path to existing `ffbonded.itp` to be read and parsed.** 
NOTE: There are unique functions for reading both forcefield files and CGenFF outputs. There is some error handling built in to account for missing improper dihedral parameters (common) and small formatting changes, but A POORLY ORGANIZED `ffbonded.itp` FILE WILL CAUSE ISSUES! The output format of this script will fit in well with standard CHARMM forcefields. 
 #### `--output`: **Path to desired output file** 
Output file format closely matches standard CHARMM forcefields. If there is already a file of that name, the option will be given to continue by overwriting or exiting. Header columns are always written for each parameter in `ffbonded.itp`, regardless if any unique parameters are printed out. Units will be indicated where necessary in brackets. 
  <br />
 #### `--kJ`: (ON/OFF) Optional flag that specifies a unit conversion. Units will remain kcal and Å by default, but `--kJ` will convert to kJ and nm units used by GROMACS.
 This adds extra functionality if you want to just check for duplicate entries in your `ffbonded.itp` file and keep the units the same. Units are indicated in headers. 
# Demo: 
  Inside `demo/` there are example `ffbonded.itp` files from different CHARMM builds, a clean, protonated .mol2 file that is compatible with CGenFF, and corresponding outputs from CGenFF. Below is a demonstration of using `cgen2gmx` with sample data, along with a short tutorial on using CGenFF. 

------------------------------------------------------------------------------------------------------------------------------
1. **Generate CGenFF output**
 
   Input a properly formatted and protonated .mol2 file to [CGenFF](https://cgenff.com/) (recently moved to SilcsBio) and retrieve the .str file with parameters. Check this file carefully for poor estimations and high penalty approximations! This script is flexible enough to recognize most commented-out and unnecessary lines, but the output of the .str file can vary. If there are any issues with parse_cgen(), make sure all lines that aren't parameters are commented out.

   **Be sure to include ALL parameters, not just the parameters that aren't already in CHARMM!** SilcsBio has no knowledge of what parameters exist in your CHARMM build, especially if an older version is being used. 

3. **Familiarize yourself with `ffbonded.itp` file and prepare to add new parameters**

    Understand how `ffbonded.itp` is formatted and which entries go where. These are large files with standard formatting rules and any unexpected lines or inconsistencies may cause issues with parse_ff(). Sample files are provided to show original up to date .itp files, as well as a slightly modified .itp to demonstrate the formatting rules.

4. **Use the cgen2gmx.py script**
   
   Generate an output file that contains unique entries-new entries to be added to the existing .itp file that will not clash with existing entries-with specified units and column order for CHARMM forcefields for GROMACS.
   From the directory containing `cgen2gmx.py` if repo is cloned:
   > python cgen2gmx.py --itp demo/ffbonded-36m_jul2022.itp --cgen demo/CGEN_OUTPUT/CRO_ex.str --output cgen2gmx_DEMO.dat --kJ

   
   **if installed via pip, `cgen2gmx` allows you to interact with the module anywhere. Enter `cgen2gmx --help` for more information.**

   This command takes in a current, unmodified `ffbonded.itp` file, searches for duplicates found in a raw CGenFF output file, and outputs unique entries in the proper order and units for GROMACS simulations (kJ and nm).


   > python cgen2gmx.py --itp demo/ffbonded-36m_jul2022.itp --cgen demo/CGEN_OUTPUT/CRO_ex.str --output cgen2gmx_DEMO.dat

   This command does the exact same thing as the previous one, but leaves the units as-is from CGenFF (kcal and Å). This should add flexibility in case this script needs to be adapted for other MD engines. However, the header column order is hardcoded and fixed to follow the format present in CHARMM forcefields.

5. **Add unique entries to `ffbonded.itp` with a text editor**

   `cgen2gmx_DEMO.dat` will contain the unique entries, complete with consistent formatting to blend in with existing entires and headers with units. Headers will be printed regardless if there are any unique entries. Use your text editor of choice to add the unique entries under the proper bracketed heading, [ bondtypes ], [ angletypes ], and so forth.  


# Contributing
 If this was useful to you in any way, whether parameterizing small molecules for a research project or just to see how someone else accomplished this, please let me know! I used a clunkier, less flexible version of this code to deal with MD parameters for chromophores in fluorescent proteins for my MS thesis. Polishing and publishing this script was a great exercise in building reusable code, as well as an opprotunity to showcase my Python abilities. 

 While `parse_cgen()` and `classes.py` handle the atom types and topolgy information found in the CGenFF output, that information is not used by the script. There are issues translating the atom types used by different MD engines, and the topologies found in CGenFF are GENERALLY compatible. Future releases may target this information for greater functionality.

 If you think this this repository would fit in well with an existing codebase, please let me know and I would be happy to contribute to more estabished software packages for visibility!
