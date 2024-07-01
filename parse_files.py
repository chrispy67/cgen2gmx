import logging as log
import pandas as pd
import re
from classes import MolecularData, ForceFieldInfo

##to do:
#re-figure out how to handle the data, try to think about future use cases. 
#environment or .yml for for script
#demo




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