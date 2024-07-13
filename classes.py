import logging as log
import pandas as pd
import re

class MolecularData:
    def __init__(self):
        self.data = {
            'merged.rtp': {
                'atom': [],
                'connectivity': []
            },
            'ffbonded.itp': {
                'bonds': [],
                'angles': [],
                'dihedrals': [],
                'impropers': []
            }
        }

    def add_atom(self, atom_data):
        self.data['merged.rtp']['atom'].append(atom_data)

    def add_connectivity(self, connectivity_data):
        self.data['merged.rtp']['connectivity'].append(connectivity_data)

    def add_bonds(self, bond_data):
        self.data['ffbonded.itp']['bonds'].append(bond_data)
    def add_angles(self, angle_data):
        self.data['ffbonded.itp']['angles'].append(angle_data)

    def add_dihedrals(self, dihedral_data):
        self.data['ffbonded.itp']['dihedrals'].append(dihedral_data)

    def add_impropers(self, improper_data):
        self.data['ffbonded.itp']['impropers'].append(improper_data)



#retrieve data from the class w/ converted units
    def get_atoms(self):
        return pd.DataFrame(self.data['merged.rtp']['atom'])

    def get_connectivity(self):
        return pd.DataFrame(self.data['merged.rtp']['connectivity'])

    def get_bonds(self):
        df = pd.DataFrame(self.data['ffbonded.itp']['bonds'])

        b0_factor = 0.1 #Å -> nm
        kb_factor = 836.8 # kcal/mol -> kJ/mol

        df['b0'] = pd.to_numeric(df['b0'], errors='coerce')
        df['b0'] *= b0_factor

        df['kb'] = pd.to_numeric(df['kb'], errors='coerce')
        df['kb'] *= kb_factor
        return pd.DataFrame(df)

#angles: ;      i        j        k  func       theta0       ktheta      ub0       kub

    def get_angles(self):
        df = pd.DataFrame(self.data['ffbonded.itp']['angles'])
        theta_factor = 1
        ktheta_factor = 8.368 
        r0_factor = 0.1 #Å -> nm
        kub_factor = 836.8 #kcal/mol -> kJ/mol

        #another way to convert all values from strings to floats
        df[['ktheta', 'theta0', 'kub', 'ub0']] = df[['ktheta', 'theta0', 'kub', 'ub0']].apply(pd.to_numeric, errors='coerce')
        df['theta0'] *= theta_factor #for redundancy, units are same (deg)
        df['ktheta'] *= ktheta_factor
        df['ub0'] *= r0_factor
        df['kub'] *= kub_factor
        return pd.DataFrame(df)

    def get_dihedrals(self):
        df = pd.DataFrame(self.data['ffbonded.itp']['dihedrals'])
        theta_0_factor = 1
        kphi_factor = 4.184 #kcal -> kJ

        df[['kphi', 'phi0']] = df[['kphi', 'phi0']].apply(pd.to_numeric, errors='coerce')
        df['kphi'] *= kphi_factor
        df['phi0'] *= theta_0_factor
        return pd.DataFrame(df)
# ;   i       j       k       l       func    phi0 [deg]    kphi [kJ/mol]          mult


    def get_impropers(self):
        df = pd.DataFrame(self.data['ffbonded.itp']['impropers'])
        keta_factor = 8.368
        eta0_factor = 1
        try:
            df[['keta', 'eta0']] = df[['keta', 'eta0']].apply(pd.to_numeric, errors='coerce')

            df['keta'] *= keta_factor
            df['eta0'] *= eta0_factor
            return pd.DataFrame(df)
        
        except KeyError:
            print('no improper dihedrals parameterized!')

#impropers ;      i        j        k        l  func         eta0 [deg]         keta [kJ/mol/rad^2]


#nearly identical class to Molecular Data

class ForceFieldInfo:
    def __init__(self):
        self.data = {
            'merged.rtp': {
                'atom': [],
                'connectivity': []
            },
            'ffbonded.itp': {
                'bonds': [],
                'angles': [],
                'dihedrals': [],
                'impropers': []
            }
        }
        self.headers = {
            'bonds': str('i    j     func  b0[nm]  kb[kJ/mol/nm^2]'),
            'angles': str('i   j  k  func  theta0  ktheta  ub0 kub'),
            'dihedrals': str('i j k l  func  phi0   kphi  mult'),
            'impropers': str('i j k l  func  phi0   kphi  mult')
        }

    def add_atom(self, atom_data):
        self.data['merged.rtp']['atom'].append(atom_data)

    def add_connectivity(self, connectivity_data):
        self.data['merged.rtp']['connectivity'].append(connectivity_data)

    def add_bonds(self, bond_data):
        self.data['ffbonded.itp']['bonds'].append(bond_data)

    def add_angles(self, angle_data):
        self.data['ffbonded.itp']['angles'].append(angle_data)

    def add_dihedrals(self, dihedral_data):
        self.data['ffbonded.itp']['dihedrals'].append(dihedral_data)

    def add_impropers(self, improper_data):
        self.data['ffbonded.itp']['impropers'].append(improper_data)

    def get_atoms(self):
        return pd.DataFrame(self.data['merged.rtp']['atom'])

    def get_connectivity(self):
        return pd.DataFrame(self.data['merged.rtp']['connectivity'])

    def get_bonds(self):
        return pd.DataFrame(self.data['ffbonded.itp']['bonds'])

    def get_angles(self):
        return pd.DataFrame(self.data['ffbonded.itp']['angles'])

    def get_dihedrals(self):
        return pd.DataFrame(self.data['ffbonded.itp']['dihedrals'])

    def get_impropers(self):
        return pd.DataFrame(self.data['ffbonded.itp']['impropers'])
