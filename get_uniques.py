"""
Tl;dr:
    I am trying to simulate fluorescent proteins with molecular dynamics simualations, including their autocatalyzed fluorophore.
    This is an issue because the fluorophores are often NOT parameterized out of the box from Charmm36 and will not be recognized during grompp
    HOWEVER, many of the bonds and atom types are similar enough to other systems that have been paramterized, so this is how they are calculated. 
    
    This comes with THREE issues:
        1. CGenFF does not know what atoms types are or are not paramertarized already, leading to duplicates.
        2. Including duplicates in the forcefield paramters causes an error with grompp that I don't want to bypass with -maxwarn flag
            I would much rather keep all the settings standard than overwrite with approximations. 
        3. The CGenFF program outputs all the parameters needed for Charmm36 in different units AND a different order. 
"""

""" EXAMPLES FROM charmm36/merged.rtp
[ EDO ]
    [ atoms ]
    C1     CG321   0.041    0
    O1     OG311  -0.640    1
    C2     CG321   0.041    2
    O2     OG311  -0.640    3
    H01    HGA2    0.090    4
    H02    HGA2    0.090    5
    H03    HGP1    0.419    6
    H04    HGA2    0.090    7
    H05    HGA2    0.090    8
    H06    HGP1    0.419    9
    [ bonds ]
    C1    O1
    C1    C2
    C1    H01
    C1    H02...
"""

"""EXAMPLES FROM charmm36/ffbonded.itp
[ bondtypes ]
; i    j       func     b0(A)        kb**
; CR2 parameters
CG2R52 CG2R53    1      0.15080     220915.20
CG2R52 CG321     1      0.15000     246856.00
CG2R52 NG2O1     1      0.14020     192464.00
CG2R53 CG321     1      0.15000     355640.00
CG2R53 NG2O1     1      0.14020     192464.00 ...

[ angletypes ]
; Original Charmm
;i        j        k    func     theta0       ktheta          ub0          kub
H        NH2      CT1     5   111.000000   418.400000   0.00000000         0.00
H        NH2      CT2     5   111.000000   418.400000   0.00000000         0.00
NH2      CT1      CT1     5   110.000000   566.513600   0.00000000         0.00
NH2      CT1      CT2     5   110.000000   566.513600   0.00000000         0.00
NH2      CT1      CT3     5   110.000000   566.513600   0.00000000         0.00
CT1       CD      OH1     5   110.500000   460.240000   0.00000000         0.00
CT3      CT1       CD     5   108.000000   435.136000   0.00000000         0.00
NH2      CT1      HB1     5   109.500000   317.984000   0.21400000     41840.00 ...

[ dihedrals ]
;   i       j       k       l       func    phi0 [deg]    kphi [kJ/mol]          mult

H      NH2      CT1      CT1     9     0.000000     0.000000     1
H      NH2      CT1        C     9     0.000000     0.000000     1
H      NH2      CT2        C     9     0.000000     0.000000     1
H      NH2      CT1      HB1     9     0.000000     0.460240     3
H      NH2      CT2      HB2     9     0.000000     0.460240     3
H      NH2      CT1      CT2     9     0.000000     0.460240     3
H      NH2      CT1      CT3     9     0.000000     0.460240     3 ...

[ dihedral types] ;impropers
;      i        j        k        l  func         eta0 [deg]         keta [kJ/mol/rad^2]
HE2      HE2      CE2      CE2     2     0.000000    25.104000
HR1      NR1      NR2     CPH2     2     0.000000     4.184000
HR1      NR2      NR1     CPH2     2     0.000000     4.184000
HR3     CPH1      NR1     CPH1     2     0.000000     4.184000
HR3     CPH1      NR2     CPH1     2     0.000000     4.184000
HR3     CPH1      NR3     CPH1     2     0.000000     8.368000
HR3      NR1     CPH1     CPH1     2     0.000000     4.184000

"""

import logging as log
import pandas as pd
from classes import MolecularData, ForceFieldInfo
from parse_files import parse_cgen, parse_ff
import re
import os
import click
import sys

#to-do:
# - put units in all classes
# - i/o for files up HERE, will make it easy for end user. 
# - put in some debug lines with the logger? putting time in here might save time later.
# - I want to make it print out stuff REALLY consistently, whitespace and everything
# - Something that can detect large errors or poor predictions??


pwd = '{}/'.format(sys.path[0])
cgen = parse_cgen('CGFF-CR2_output.dat') #get data by treating cgen as the class.
ff = parse_ff('ffbonded.itp')

#how to loop through the dictionary and where the keys are addressable.
def iterate_nested_dict(nested_dict, parent_key=''):
    for key, value in nested_dict.items():
        full_key = f"{parent_key}.{key}" if parent_key else key
        if isinstance(value, dict):
            iterate_nested_dict(value, full_key)
        else:
            print(f"{full_key}")


#this is how I can access data from the nested list with some over complicated files
##EXAMPLE
# print(cgen.data['ffbonded.itp'].keys())
# print(cgen.data['ffbonded.itp']['bonds'])
# print(cgen.data['ffbonded.itp']['bonds'][0])

#same as FF data!
# print(ff.data['ffbonded.itp']['bonds'][2])

def get_uniques(ff, cgen, index_columns, param):
    try:
        ff_df = pd.DataFrame(ff)
        cgen_df = pd.DataFrame(cgen)
        merged = cgen_df.merge(ff_df, on=index_columns, how='left', indicator=True)
        
        # Unique entries 
        unique_entries = merged[merged['_merge'] == 'left_only']

        # Drop columns from ff_df that are not in index_columns
        columns_to_drop = ff_df.columns.difference(index_columns)
        columns_to_drop = [col for col in columns_to_drop if col in unique_entries.columns]
        
        # Curate dataset to be easier to work with
        unique_entries = unique_entries.drop(columns=columns_to_drop)
        unique_entries = unique_entries.drop(columns=['_merge', 'mult_y'], errors='ignore')  # Drop unnecessary columns

        # Remove columns and reset index
        unique_entries = unique_entries.rename(columns=lambda x: x.rstrip('_x'))
        unique_entries.reset_index(drop=True, inplace=True)

        return pd.DataFrame(unique_entries)
        

    except AttributeError as e:
        print(f'AttributeError: {e}')
        print(f'cannot merge empty array, possibly there are no entires for {param}?')
    except KeyError as e:
        if len(ff_df.columns) == 0:
            print(f'KeyError: {e}')
            print(f'length of merged array for {param} is 0, no unique entries detected?')
        else:
            print(f'KeyError: {e}')
            print(f'ff_df columns: {ff_df.columns}')
            print(f'index_columns: {index_columns}')

unique_bonds = get_uniques(ff.get_bonds(), 
    cgen.get_bonds(), 
    ff.headers['bonds'].split()[0:2],
    param='bonds')

unique_angles = get_uniques(ff.get_angles(), 
    cgen.get_angles(),
    ff.headers['angles'].split()[0:3],
    param='angles')

unique_dihedrals = get_uniques(ff.get_dihedrals(),
    cgen.get_dihedrals(),
    ff.headers['dihedrals'].split()[0:4],
    param='dihedrals')

unique_impropers = get_uniques(ff.get_impropers(), 
    cgen.get_impropers(),
    ff.headers['impropers'].split()[0:4],
    param='impropers')

    # unit_regex = r'\[.*?\]'  # remove anything between brackets
    # columns = re.sub(unit_regex, '', headers).strip().split()


def update_charmm(df, headers, param):
    columns = ff.headers[param].split()
    unit_regex = r'\[.*?\]'
    columns_regex = re.sub(unit_regex, '',headers).strip().split()
    # print(columns_regex)

#this corrects for the issue that the code has with empty arrays and no parsed files. 
    if df is not None and hasattr(df, 'columns'):
        header_columns = columns_regex
    else:
        header_columns = []

    with open('ffbonded-edits.dat', 'w') as f:
        f.write(';' + str(header_columns) + '\n')
        if df is not None and hasattr(df, 'iterrows'):
            for index, row in df.iterrows():
                line = '   '.join(str(item) if not isinstance(item, list) else ' '.join(map(str, item)) for item in row[header_columns])
                f.write(line + '\n')

##i       j       k       l       func    phi0 [deg]    kphi [kJ/mol]          mult
update_charmm(unique_bonds, str(ff.headers['bonds']), 'bonds')
update_charmm(unique_angles, str(ff.headers['angles']), 'angles')
update_charmm(unique_dihedrals, str(ff.headers['dihedrals']), 'dihedrals')
#update_charmm(unique_impropers, str(ff.headers['impropers']), 'impropers')