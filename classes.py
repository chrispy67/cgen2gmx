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

#store data in the class and empty arrays above
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
        kb_factor = 836.8 #kcal/molÅ^-2 -> kJ/mol*nm^-2 #CORRECT FOR FUNC_1

        df['b0'] = pd.to_numeric(df['b0'], errors='coerce')
        df['b0'] *= b0_factor

        df['kb'] = pd.to_numeric(df['kb'], errors='coerce')
        df['kb'] *= kb_factor

        unit_dict = {
            'b0': 'nm',
            'kb': 'kJ/mol*nm^-2'}


        units_df = pd.DataFrame({key: [value] for key, value in unit_dict.items()}) # (1, 4)
        units_df = pd.concat([units_df] * len(df), ignore_index=True) # (len(parsed_entries), 4) each column is the same value
        
        # Rename the unit columns to distinguish from the original columns
        units_df = units_df.add_suffix('_unit')
        
        df_with_units = pd.concat([df, units_df], axis=1)# (len(parsed_entries), 11)

        return df_with_units

#angles: ;      i        j        k  func       theta0       ktheta      ub0       kub

    def get_angles(self):
        df = pd.DataFrame(self.data['ffbonded.itp']['angles'])
        theta_factor = 1 #deg -> deg
        ktheta_factor = 8.368 #kcal -> kJ * func_5_factor
        r0_factor = 0.1 #Å -> nm
        kub_factor = 836.8 #kcal/molÅ^-2 -> kJ/mol*nm^-2 #CORRECT FOR FUNC_5

        #another way to convert all values from strings to floats
        df[['ktheta', 'theta0', 'kub', 'r0']] = df[['ktheta', 'theta0', 'kub', 'r0']].apply(pd.to_numeric, errors='coerce')
        df['theta0'] *= theta_factor #for redundancy, units are same (deg)
        df['ktheta'] *= ktheta_factor
        df['kub'] *= kub_factor
        df['r0'] *= r0_factor
        
        
        #UNITS
        unit_dict = {
            'ktheta': 'kJ/mol',
            'theta0': 'deg',
            'kub': 'kJ/mol*nm^-2',
            'r0': 'nm'} 
        # Create a DataFrame to hold the units
        units_df = pd.DataFrame({key: [value] for key, value in unit_dict.items()}) # (1, 4)
        units_df = pd.concat([units_df] * len(df), ignore_index=True) # (len(parsed_entries), 4) each column is the same value
        
        # Rename the unit columns to distinguish from the original columns
        units_df = units_df.add_suffix('_unit')
        
        df_with_units = pd.concat([df, units_df], axis=1)# (len(parsed_entries), 11)

        return df_with_units

    def get_dihedrals(self):
        df = pd.DataFrame(self.data['ffbonded.itp']['dihedrals'])
        phi_0_factor = 1 #deg -> deg
        kphi_factor = 4.184 #kcal -> kJ

        df[['kphi', 'phi0', ]] = df[['kphi', 'phi0']].apply(pd.to_numeric, errors='coerce')
        df['kphi'] *= kphi_factor #TO BE MULTIPLIED BY MULTI FOR THIS FUNCTION
        df['phi0'] *= phi_0_factor

        unit_dict = {
            'kphi': 'kJ/mol',
            'phi0': 'deg'}
        units_df = pd.DataFrame({key: [value] for key, value in unit_dict.items()}) # (1, 4)
        units_df = pd.concat([units_df] * len(df), ignore_index=True) # (len(parsed_entries), 4) each column is the same value
        
        # Rename the unit columns to distinguish from the original columns
        units_df = units_df.add_suffix('_unit')
        
        df_with_units = pd.concat([df, units_df], axis=1)# (len(parsed_entries), 11)

        return df_with_units


    def get_impropers(self):
        df = pd.DataFrame(self.data['ffbonded.itp']['impropers'])
        keta_factor = 8.368 #kcal/mol -> kJ/mol *2
        eta0_factor = 1 #deg -> deg 
        try:
            df[['keta', 'eta0']] = df[['keta', 'eta0']].apply(pd.to_numeric, errors='coerce')

            df['keta'] *= keta_factor
            df['eta0'] *= eta0_factor
            return pd.DataFrame(df)
        
        except KeyError:
            print('no improper dihedrals parameterized!')
        unit_dict = {
            'kphi': 'kJ/mol',
            'phi0': 'deg'}
        units_df = pd.DataFrame({key: [value] for key, value in unit_dict.items()}) # (1, 2)
        units_df = pd.concat([units_df] * len(df), ignore_index=True) # (len(parsed_entries), 2) each column is the same value
        
        # Rename the unit columns to distinguish from the original columns
        units_df = units_df.add_suffix('_unit')
        
        df_with_units = pd.concat([df, units_df], axis=1)# (len(parsed_entries), 11)

        return df_with_units


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

