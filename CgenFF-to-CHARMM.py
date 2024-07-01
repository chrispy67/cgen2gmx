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
import re
from classes import MolecularData, ForceFieldInfo

##to do:
#re-figure out how to handle the data, try to think about future use cases. 
#environment or .yml for for script
#demo

##to fix:
#There is a key  or some attribute issue here, but I like the separation of grouping
# self.data{} by filetype! This will allow for good loops and more comprehensive 
# and faster data handling. These could be longer files, so speed is good. 



def parse_CGenFF(file_path):
    molecular_data = MolecularData()
    current_section = None
    
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            print(current_section)
            if not line or line.startswith(('*', '#', '!')):
                continue
            elif not line.strip():
                current_section = None

            ####----merged.itp----####
            if line.startswith('ATOM'):
                current_section = 'ATOM'
            elif line.startswith('BOND '): #SINGULAR-merged.rtp
                current_section = 'BOND'  #SINGULAR-merged.rtp
            elif line.startswith('IMPR'):
                current_section = 'IMPR'
            
            ####----ffbonded.itp----####
            elif line.startswith('BONDS'): #PLURAL-ffbonded.itp
                current_section = 'BONDS'  #PLURAL-ffbonded.itp
            elif line.startswith('ANGLES'):
                current_section = 'ANGLES'
            elif line.startswith('DIHEDRALS'):
                current_section = 'DIHEDRAL'
            elif line.startswith('IMPROPERS'):
                current_section = 'IMPROPER'
            # Add more sections as needed

            #FOR SELECTIONS
            
            if current_section == 'ATOM':
                atom_data = line.split()
                #['ATOM', 'H11', 'HGR61', '0.115', '!', '0.000']
                row = {
                    'i': atom_data[1],
                    'j': atom_data[2],
                    'charge': atom_data[3]
                }
                molecular_data.add_atom(row)


#HANDLES THE PLURAL AND SINGLUAR VERSIONS THAT CgenFF likes to use.
#similar treatment to 
            elif current_section == 'BOND ' and len(line.split()) != 1:
                bond_data = line.split()
                if len(bond_data) == 3:
                    row = {
                    'i': bond_data[1],
                    'j': bond_data[2]
                    }
                    molecular_data.add_connectivity(row)
                else:
                    row = {
                    'i': bond_data[1],
                    'j': bond_data[2],
                    'func': int(1),
                    'kb': bond_data[3], ###UNITS spring constant
                    'b0': bond_data[4]  ###UNITS A -> nm
                    }
                    molecular_data.add_bonds(row)

#angles: ;      i        j        k  func       theta0       ktheta      ub0       kub

            elif current_section == 'ANGLES' and len(line.split()) != 1 :
                angle_data = line.split()
                if angle_data[6] == '!':
                    row = {
                        'i': angle_data[1],
                        'j': angle_data[2],
                        'k': angle_data[3],
                        'func': 5,
                        'ktheta': angle_data[4],
                        'theta0': angle_data[5],
                        'kub': float(0.00),
                        'ub0': float(0.00)
                    }
                else:
                        row = {
                        'i': angle_data[1],
                        'j': angle_data[2],
                        'k': angle_data[3],
                        'ktheta': angle_data[4],
                        'theta0': angle_data[5],
                        'kub': angle_data[6],
                        'ub0': angle_data[7]
                    }
                molecular_data.add_angles(row)
            
            elif current_section == 'DIHEDRAL' and len(line.split()) != 1:
                dihedral_data = line.split()
                row = {
                    'i': dihedral_data[1],  #STANDARD DIHEDRALS HAVE FUNCTION 9
                    'j': dihedral_data[2],
                    'k': dihedral_data[3],
                    'l': dihedral_data[4],
                    'func': int(9),
                    'kphi': dihedral_data[5],
                    'mult': dihedral_data[6],
                    'phi0': dihedral_data[7]
                }
                molecular_data.add_dihedrals(row)

# ;   i       j       k       l       func    phi0 [deg]    kphi [kJ/mol]          mult


            elif current_section == 'IMPROPER':
                #OFTEN EMPTY
                #BUT CONSIDER ERROR HANDLING THE REMAINDER OF THE SECTIONS FOR NO ENTRIES
                
                improper_data = line.split()
                row = {
                    'i': improper_data[1],    #IMPROPERS HAVE FUNCTION 2
                    'j': improper_data[2],
                    'k': improper_data[3],
                    'l': improper_data[4],
                    'func': int(2),
                    'mult': improper_data[5],
                    'keta': improper_data[6],
                    'eta0': improper_data[7]
                }
                
                molecular_data.add_impropers(row)

    return molecular_data

parse_CGenFF('CGFF-CR2_output.dat')

# i       j       k       l       func    phi0 [deg]    kphi [kJ/mol]          mult

def parse_FF(ff_file): ##WIP
    #this is going to be the tricky part: searching the ffbonded.itp file 
    #for matching strings.
    selection = False
    ff_data = ForceFieldInfo()

    with open(ff_file, 'r') as file:
        for line in file:
            try:
                data = line.split()
            except IndexError:
                continue

            #checks for comments
            if not line.strip() or line.startswith(('*', '#', ';', '!')):
                continue
            #MORE ELEGANT SOLUTION FOR DEALING WITH SELECTIONS
            if line.startswith('[ '):
                selection = line.split('[')[-1].split(']')[0].strip().lower()
                continue

            if selection == 'bondtypes':
                try:
                    bond_data = line.split()
                    existing_data = {
                        'i': bond_data[0],
                        'j': bond_data[1]
                    }
                    ff_data.add_bonds(existing_data)
                except IndexError:
                    continue
            
            elif selection == 'angletypes':
                try:
                    angle_data = line.split()
                    existing_data = {
                        'i': angle_data[0],
                        'j': angle_data[1],
                        'k': angle_data[2]
                    }
                    ff_data.add_angles(existing_data)
                except IndexError:
                    continue

#LESS THAN PERFECT SOLUTION FOR BYPASSING HEADERS, MIGHT GIVE SOME WEIRD DATAPOINTS
            elif selection == 'dihedraltypes'  and len(data) >=4 and data[4] == '9': #ACTUAL DIHEDRALS
                try:
                    dihedral_data = line.split()
                    existing_data = {
                        'i': dihedral_data[0],
                        'j': dihedral_data[1],
                        'k': dihedral_data[2],
                        'l': dihedral_data[3],
                        'mult': dihedral_data[7] 
                    }
                    ff_data.add_dihedrals(existing_data)
                except IndexError:

                    continue
            
            elif selection == 'dihedraltypes' and len(data) >= 4 and data[4] == '2': #IMPROPERS
                try:
                    improper_data = line.split()
                    existing_data = {
                        'i': improper_data[0],
                        'j': improper_data[1],
                        'k': improper_data[2],
                        'l': improper_data[3],
                        'mult':improper_data[7]
                    }
                    ff_data.add_impropers(existing_data)
                except IndexError:
                    continue

                #ONLY NEED TO CHECK i, j, k, l indexes
    return ff_data


# ff_dihedrals = pd.DataFrame(parsed_FF.get_dihedrals())
# CGEN_dihedrals = pd.DataFrame(parsed_CGen.get_dihedrals())


# def get_unique_entries(ff, cgen, index_columns):
#     try:
#         merged = cgen.merge(ff, on=index_columns, how='left', indicator=True)

#         #unique entries 
#         unique_entries = merged[merged['_merge'] == 'left_only']
#         #curating dataset to be easier to work with
#         unique_entries = unique_entries.drop(columns=ff.columns.difference  (index_columns))
#         unique_entries = unique_entries.drop(columns=['_merge'])

#         #reset index
#         unique_entries.reset_index(drop=True, inplace=True)
#         return unique_entries
#     except AttributeError:
#         print('cannot merge empty array')

# parsed_FF = parse_FF('ffbonded.itp')
# parsed_CGen = parse_CGenFF('CGFF-CR2_output.dat')

# cgen_impropers = parsed_CGen.get_impropers()

# unique_bonds = get_unique_entries(parsed_FF.get_bonds(), 
#     parsed_CGen.get_bonds(),
#     ['i', 'j'])

# unique_angles = get_unique_entries(parsed_FF.get_angles(), 
#     parsed_CGen.get_angles(), 
#     ['i', 'j', 'k'])

# unique_dihedrals = get_unique_entries(parsed_FF.get_dihedrals(),
#     parsed_CGen.get_dihedrals(),
#     ['i', 'j', 'k', 'l', 'mult'])

# unique_impropers = get_unique_entries(parsed_FF.get_impropers(), 
#     parsed_CGen.get_impropers(),
#     ['i', 'j', 'k', 'l', 'mult'])

# def update_charmm(df, headers):

#     unit_regex = r'\[.*?\]' #remove anything between brackets

#     #DF.keys() MUST MATCH THE HEADERS
#     #I changed some of this stuff quickly so it may need some tweaking
#     columns = re.sub(unit_regex, '', headers).strip().split()
#     header_columns = [col for col in columns if col in df.columns]

#     with open('ffbonded-edits.dat', 'a') as f:
#         f.write(';'+str(header_columns)+'\n')
#         for index, row in df.iterrows():
#             line = '   '.join(str(item) if not isinstance(item, list) else ' '.join(map(str, item)) for item in row[header_columns])
#             f.write(line+'\n')

# # i       j       k       l       func    phi0 [deg]    kphi [kJ/mol]          mult

# update_charmm(unique_bonds, ' ; i    j  func       b0 [nm]          kb [kJ/mol/nm^2]')
# update_charmm(unique_angles, ' ; i    j  k  func  theta0  ktheta  ub0 kub')
# update_charmm(unique_dihedrals,' i j k l  func  phi0   kphi  mult' )





