import pandas as pd
import os
import click
import sys
import argparse
from src import config
# to-do:l
# - update_charmm() can't hang with the empty set()
# - error handling for the new format_string() function.
# - NOT RECOGNIZING BACKWARDS-DEFINED PARAMETERS!
#   - unsure if there is a convention, but a,b,c,d and d,c,b,a are equivalent interactions, 
#     yet these backwards matches are not captured by get_uniques()?

###LONG GAME
# - put in some debug lines with the logger? putting time in here might save time later.
# - Something that can detect large errors or poor predictions?
#   - would need to parse this as an optional flag in parse_files.py


#may need to change, but this is how I get files
pwd = '{}/'.format(sys.path[0])

def main():
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
    from src.parse_files import parse_cgen, parse_ff
    from src.get_uniques import get_uniques, update_charmm

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

    # #
    if os.path.isfile(config.output_file):
        if click.confirm(f'{config.output_file} exists, would you like to overwrite the previous file?', default=True):
            with open(config.output_file, 'w') as fp:
                pass
        else:  
            sys.exit('canceling')

    update_charmm(unique_bonds, config.output_file, cgen.bond_header(), param='bonds')
    update_charmm(unique_angles, config.output_file, cgen.angle_header(), param='angles')
    update_charmm(unique_dihedrals, config.output_file, cgen.dihedral_header(), param='dihderals')
    update_charmm(unique_impropers, config.output_file, cgen.improper_header(), param='imrpropers')


if __name__ == '__main__':
    main() #run the script 