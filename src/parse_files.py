from src.classes import MolecularData, ForceFieldInfo

def parse_cgen(file_path):
    molecular_data = MolecularData()
    
    #selection indicates where to store data in ForceFieldInfo()
    current_selection = None

    #CGEN has some lines to skip that aren't parameters OR comments...
    #might need to add things to this list if you are having issues with reading files. 
    to_skip = ('*', '#', '!', ';', 'END ', 'RETURN', ' ', 'RESI')

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if not line or line.startswith(to_skip):
                continue
            elif not line.strip():
                current_selection = None

            ####----merged.itp----####
            if line.startswith('ATOM'):
                current_selection = 'ATOM'
            elif line.startswith('BOND '): #SINGULAR-merged.rtp
                current_selection = 'BOND'  #SINGULAR-merged.rtp
            elif line.startswith('IMPR '):
                current_selection = 'IMPR'
            
            ####----ffbonded.itp----####
            elif line.startswith('BONDS'): #PLURAL-ffbonded.itp
                current_selection = 'BONDS'  #PLURAL-ffbonded.itp
            elif line.startswith('ANGLES'):
                current_selection = 'ANGLES'
            elif line.startswith('DIHEDRALS'):
                current_selection = 'DIHEDRAL'
            elif line.startswith('IMPROPERS'):
                current_selection = 'IMPROPER'
            # Add more sections as needed

            #FOR SELECTIONS
            
            if current_selection == 'ATOM':
                atom_data = line.split()
                row = {
                    'i': atom_data[1],
                    'j': atom_data[2],
                    'charge': atom_data[3]
                }
                molecular_data.add_atom(row)


#HANDLES THE PLURAL AND SINGLUAR VERSIONS THAT CgenFF likes to use.
#similar treatment to 
            elif current_selection == 'BOND ' and len(line.split()) != 1: #avoids header
                bond_data = line.split()
                row = {
                'i': bond_data[1],
                'j': bond_data[2]
                }
                molecular_data.add_connectivity(row)
            elif current_selection == 'BONDS' and len(line.split()) != 1:
                connectivity_data = line.split()
                row = {
                'i': connectivity_data[0],
                'j': connectivity_data[1],
                'func': int(1),
                'kb': connectivity_data[2], ###UNITS spring constant
                'b0': connectivity_data[3]  ###UNITS A -> nm
                }
                molecular_data.add_bonds(row)



#angles: ;      i        j        k  func       theta0       ktheta      ub0       kub

            elif current_selection == 'ANGLES' and len(line.split()) != 1 :
                angle_data = line.split()
                # Only some entries have Urey-Bradley potentials. Otherwise, the entry for
                #   these columns is zero.
                if angle_data[6] == '!':
                    row = {
                        'i': angle_data[0],
                        'j': angle_data[1],
                        'k': angle_data[2],
                        'func': 5,
                        'ktheta': angle_data[3],
                        'theta0': angle_data[4],
                        'kub': float(0.00),
                        'ub0': float(0.00)
                    }
                else:
                        row = {
                        'i': angle_data[0],
                        'j': angle_data[1],
                        'k': angle_data[2],
                        'func': 5, 
                        'ktheta': angle_data[3],
                        'theta0': angle_data[4],
                        'kub': angle_data[5],
                        'ub0': angle_data[6]
                    }
                molecular_data.add_angles(row)
            
            elif current_selection == 'DIHEDRAL' and len(line.split()) != 1:
                dihedral_data = line.split()
                row = {
                    'i': dihedral_data[0],  #STANDARD DIHEDRALS HAVE FUNCTION 9
                    'j': dihedral_data[1],
                    'k': dihedral_data[2],
                    'l': dihedral_data[3],
                    'func': int(9),
                    'kphi': dihedral_data[4],
                    'multi': dihedral_data[5],
                    'phi0': dihedral_data[6]
                }
                molecular_data.add_dihedrals(row)


# ;   i       j       k       l       func    phi0 [deg]    kphi [kJ/mol]          mult


            elif current_selection == 'IMPROPER' and len(line.split()) != 1:
                
                improper_data = line.split()
                row = {
                    'i': improper_data[0],    #IMPROPERS HAVE FUNCTION 2
                    'j': improper_data[1],
                    'k': improper_data[2],
                    'l': improper_data[3],
                    'kphi': improper_data[4],
                    #'multi': improper_data[5], # there are NO multiplicity functions here
                    'func': str(2),
                    'phi0': improper_data[6]
                }
                
                molecular_data.add_impropers(row)

    return molecular_data

def parse_ff(ff_file): ##WIP
    selection = False
    ff_data = ForceFieldInfo()

    with open(ff_file, 'r') as file:
        for line in file:
            try:
                data = line.split()
            except IndexError:
                continue

            # Checks for comments
            if not line.strip() or line.startswith(('*', '#', ';', '!')):
                continue
            if line.startswith('[ '):
                selection = line.split('[')[-1].split(']')[0].strip().lower()
                continue

            # Adheres to standard CHARMM forcefield headers, but subject to change.
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
            
            # Adheres to standard CHARMM forcefield headers, but subject to change.
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

            # Adheres to standard CHARMM forcefield headers, but subject to change.
            # 'func' being 9 is the only thing that makes this an improper.  
            elif selection == 'dihedraltypes' and data[4] == '9': #ACTUAL DIHEDRALS
                try:
                    dihedral_data = line.split()
                    existing_data = {
                        'i': dihedral_data[0],
                        'j': dihedral_data[1],
                        'k': dihedral_data[2],
                        'l': dihedral_data[3],
                        'multi': dihedral_data[7] 
                    }
                    ff_data.add_dihedrals(existing_data)
                except IndexError:
                    continue

            # Adheres to standard CHARMM forcefield headers, but subject to change.
            elif selection == 'dihedraltypes' and data[4] == '2': #IMPROPERS
                try:
                    improper_data = line.split()
                    existing_data = {
                        'i': improper_data[0],
                        'j': improper_data[1],
                        'k': improper_data[2],
                        'l': improper_data[3],
                        'func':improper_data[4]
                    }
                    ff_data.add_impropers(existing_data)
                except IndexError:
                    continue

    return ff_data