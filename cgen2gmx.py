import pandas as pd
import os
import click
import sys
import argparse
import config

# to-do:
# - tie units to appropriate keys for header columns
# - error handling for the new format_string() function.
# - what about an option where I can append all of the new values
#   to the existing ffbonded.itp in a recognizable way. must make backups
# - NOT RECOGNIZING BACKWARDS-DEFINED PARAMETERS!
#   - unsure if there is a convention, but a,b,c,d == d,c,b,a are equivalent interactions, 
#     yet these backwards matches are not captured by get_uniques() 

###LONG GAME
# - i/o for files up HERE, will make it easy for end user. 
# - put in some debug lines with the logger? putting time in here might save time later.
# - Something that can detect large errors or poor predictions?
#   - would need to parse this as an optional flag in parse_files.py


#may need to change, but this is how I get files
pwd = '{}/'.format(sys.path[0])

#general help string
parser = argparse.ArgumentParser(description="""
A script that takes in raw outputs from CHARMM CgenFF forcefield parameter generator,
seaches for duplicate entries in a provided ffbonded.itp file, and neatly prints
out all of the needed parameters with headers and correct units
""")

#input charmm file
parser.add_argument('--itp', type=str, default='ffbonded.itp', 
    help='Path to ffbonded.itp file to be read from CHARMM forcefield')

#input cgen parameter file
parser.add_argument('--cgen', type=str,
    help='?Raw output file from CgenFF by SilcsBio. PLEASE CHECK CONTENTS OF FILE CAREFULLY')

parser.add_argument('--output', type=str, 
    help="""File for unique forcefield parameters to be printed to. If file is present, 
    an additional prompt will be presented to exit""")

#boolean flag for units, need another for kcals!
#this would be good to handle in the classes.py, changing the output of the df units
#currently not functional
unit_inputs = ['--kj', '--kJ', '--kJ/mol', '--kJ*mol^-1']
parser.add_argument(*unit_inputs, dest='convert_to_kj', action='store_true',
    help='Specifying units (CgenFF uses kcal/mol, but GROMACS needs kJ/mol.')


args = parser.parse_args()
config.input_file_CHARMM = args.itp
config.input_file_CGEN = args.cgen 
config.output_file = args.output 
config.convert_to_kj = args.convert_to_kj
print('from cgen2gmx', config.convert_to_kj, config.input_file_CGEN)


if config.convert_to_kj:
    print('Converting to kJ/mol with values:', config.convert_to_kj)
else:
    print('No unit conversion specified.')


#importing functions and files AFTER parsing the arguments is VERY important here
from parse_files import parse_cgen, parse_ff
from get_uniques import get_uniques, update_charmm

ff = parse_ff(config.input_file_CHARMM) #class
cgen = parse_cgen(config.input_file_CGEN) #class

unique_bonds = get_uniques(ff.get_bonds(), 
    cgen.get_bonds(), 
    list(ff.get_bonds().keys()),
    param='bonds')

unique_angles = get_uniques(ff.get_angles(), 
    cgen.get_angles(),
    list(ff.get_angles().keys()),
    param='angles')

unique_dihedrals = get_uniques(ff.get_dihedrals(),
    cgen.get_dihedrals(),
    list(ff.get_dihedrals().keys()),
    param='dihedrals')

unique_impropers = get_uniques(ff.get_impropers(), 
    cgen.get_impropers(),
    list(ff.get_impropers().keys()),
    param='impropers')

#
if os.path.isfile(config.output_file):
    if click.confirm(f'{config.output_file} exists, would you like to overwrite the previous file?', default=True):
        with open(config.output_file, 'w') as fp:
            pass
    else:  
        sys.exit('canceling')


update_charmm(unique_bonds, config.output_file, cgen.bond_header())
update_charmm(unique_angles, config.output_file, cgen.angle_header())
update_charmm(unique_dihedrals, config.output_file, cgen.dihedral_header())
update_charmm(unique_impropers, config.output_file, cgen.improper_header())

