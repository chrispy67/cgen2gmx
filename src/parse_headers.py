import pandas as pd
import os
import click
import sys
import argparse
import subprocess
from src import config, parse_headers

params = ['bonds', 'angles', 'dihedrals, impropers'] #all 4 variables

# ffbonded.itp is target and these are default values.

# I need to figure out where this fits in the order of the scripts. Rememeber, these values are getting REORDERED in get_uniques(). 
#   can I adjust the header values directly in the dataframe that is printed?
#   
#   Types here are really important.
bonds_dict = {
    'row0': str, # i
    'row1': str, # j
    'row2': int, # b0
    'row3': int, # kb
}

angle_dict = {
    'row0': str, # i
    'row1': str, # j
    'row2': str, # k
    'row3': int, # func
    'row4': int, # theta0 DOUBLE CHECK THE ORDER OF THESE
    'row5': int, # ktheta DOUBLE CHECK THE ORDER OF THESE
    'row6': int, # kub DOUBLE CHECK THE ORDER OF THESE
    'row7': int, # ub0 DOUBLE CHECK THE ORDER OF THESE
}

dihedral_dict = {
    'row0': str, # i 
    'row1': str, # j
    'row2': str, # k
    'row3': str, # l
    'row4': int, # func,
    'row5': int, # phi0
    'row6': int, # kphi
    'row7': int  # mult
}

improper_dict = {
    'row0': str, # i 
    'row1': str, # j
    'row2': str, # k
    'row3': str, # l
    'row4': int, # func,
    'row5': int, # phi0
    'row6': int, # kphi
    'row7': int  # mult
}


# function that loops through these dicts and assigns values to be passed at the end of the script. 
#   should it be a string?
def headers_interactive():
    header = 'testing 1 2 3'
    return header