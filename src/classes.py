import logging as log
import pandas as pd
import re
import os
from src import config


# import user flags stored in config file
input_file_CHARMM = config.input_file_CHARMM
input_file_CGEN = config.input_file_CGEN
output_file = config.output_file
unit_convert = config.convert_to_kj


class MolecularData:
    def __init__(self):
        
        #broken down by file type for future functionality
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

#retrieve data from the class

    def get_atoms(self):
        #not used yet; the issue of unique entries is not present in merged.rtp
        return pd.DataFrame(self.data['merged.rtp']['atom'])

    def get_connectivity(self):
        #not used yet; the issue of unique entries is not present in merged.rtp
        return pd.DataFrame(self.data['merged.rtp']['connectivity'])

    def get_bonds(self):
        df = pd.DataFrame(self.data['ffbonded.itp']['bonds'])
        df['b0'] = pd.to_numeric(df['b0'], errors='coerce')
        df['kb'] = pd.to_numeric(df['kb'], errors='coerce')

        #CGEN uses kcal by default, but GROAMCS uses kJ
        if unit_convert:
            b0_factor = 0.1 #Å -> nm
            kb_factor = 836.8 #kcal/molÅ^-2 -> kJ/mol*nm^-2 #CORRECT FOR FUNC_1
            df['b0'] *= b0_factor
            df['kb'] *= kb_factor

            unit_dict = {
                'b0': 'nm',
                'kb': 'kJ/mol*nm^-2'}
        
        else:
            unit_dict = {
                'b0': 'Å',
                'kb': 'kcal/mol*Å^2'
            }

        # Rename the unit columns to distinguish from the original columns
        units_df = pd.DataFrame({key: [value] for key, value in unit_dict.items()}) # (1, 4)
        units_df = pd.concat([units_df] * len(df), ignore_index=True) # (len(parsed_entries), 4) each column is the same value
        units_df = units_df.add_suffix('_unit')
        
        df_with_units = pd.concat([df, units_df], axis=1)# (len(parsed_entries), 11)

        # df containing all entries from cgen file with units in each row of *_unit keys
        return df_with_units

#angles: ;      i        j        k  func       theta0       ktheta      ub0       kub

    def get_angles(self):
        df = pd.DataFrame(self.data['ffbonded.itp']['angles'])
        df[['ktheta', 'theta0', 'kub', 'ub0']] = df[['ktheta', 'theta0', 'kub', 'ub0']].apply(pd.to_numeric, errors='coerce')

        #CGEN uses kcal by default, but GROAMCS uses kJ
        if unit_convert:
            theta_factor = 1 #deg -> deg
            ktheta_factor = 8.368 #kcal -> kJ * func_5_factor
            ub0_factor = 0.1 #Å -> nm
            kub_factor = 836.8 #kcal/molÅ^-2 -> kJ/mol*nm^-2 #CORRECT FOR FUNC_5
            df['theta0'] *= theta_factor #for redundancy, units are same (deg)
            df['ktheta'] *= ktheta_factor
            df['kub'] *= kub_factor
            df['ub0'] *= ub0_factor

            #UNITS
            unit_dict = {
                'ktheta': 'kJ/mol',
                'theta0': 'deg',
                'kub': 'kJ/mol*nm^-2',
                'ub0': 'nm'} 
        else:
            #UNITS
            unit_dict = {
                'ktheta': 'kcal/mol',
                'theta0': 'deg',
                'kub': 'kcal/mol*Å^-2',
                'ub0': 'Å'} 

        # Rename the unit columns to distinguish from the original columns
        units_df = pd.DataFrame({key: [value] for key, value in unit_dict.items()}) # (1, 4)
        units_df = pd.concat([units_df] * len(df), ignore_index=True) # (len(parsed_entries), 4) each column is the same value
        units_df = units_df.add_suffix('_unit')
        df_with_units = pd.concat([df, units_df], axis=1)# (len(parsed_entries), 11)

        # df containing all entries from cgen file with units in each row of *_unit keys
        return df_with_units

    def get_dihedrals(self):
        df = pd.DataFrame(self.data['ffbonded.itp']['dihedrals'])
        df[['kphi', 'phi0', ]] = df[['kphi', 'phi0']].apply(pd.to_numeric, errors='coerce')

        #CGEN uses kcal by default, but GROAMCS uses kJ
        if unit_convert:
            phi_0_factor = 1 #deg -> deg
            kphi_factor = 4.184 #kcal -> kJ

            df['kphi'] *= kphi_factor #TO BE MULTIPLIED BY MULTI FOR THIS FUNCTION
            df['phi0'] *= phi_0_factor

            unit_dict = {
                'kphi': 'kJ/mol',
                'phi0': 'deg'}
        else:
            unit_dict = {
                'kphi': 'kcal/mol',
                'phi0': 'deg'}

        # Rename the unit columns to distinguish from the original columns
        units_df = pd.DataFrame({key: [value] for key, value in unit_dict.items()}) # (1, 4)
        units_df = units_df.add_suffix('_unit')
        units_df = pd.concat([units_df] * len(df), ignore_index=True) # (len(parsed_entries), 4) each column is the same value
        
        df_with_units = pd.concat([df, units_df], axis=1)# (len(parsed_entries), 11)

        # df containing all entries from cgen file with units in each row of *_unit keys
        return df_with_units


    def get_impropers(self):
        df = pd.DataFrame(self.data['ffbonded.itp']['impropers'])
        
        try: #often empty
            df[['kphi', 'phi0']] = df[['kphi', 'phi0']].apply(pd.to_numeric, errors='coerce')

        except KeyError as e:
            print(f'Key Error! {e}', 'is the parameter field entry?')
        
        #CGEN uses kcal by default, but GROAMCS uses kJ
        if unit_convert:
            keta_factor = 8.368 #kcal/mol -> kJ/mol *2
            eta0_factor = 1 #deg -> deg 
            df['kphi'] *= keta_factor
            df['phi0'] *= eta0_factor

            unit_dict = {
                'kphi': 'kJ/mol',
                'phi0': 'deg'}
        else:
            unit_dict = {
                'kphi': 'kcal/mol',
                'phi0': 'deg'}

        units_df = pd.DataFrame({key: [value] for key, value in unit_dict.items()}) # (1, 2)
        units_df = pd.concat([units_df] * len(df), ignore_index=True) # (len(parsed_entries), 2) each column is the same value
        units_df = units_df.add_suffix('_unit')        

        df_with_units = pd.concat([df, units_df], axis=1)# (len(parsed_entries), 11)

        # df containing all entries from cgen file with units in each row of *_unit keys
        return df_with_units

    #HARDCODED headers that are standard for modern CHARMM forcefields
    def bond_header(self):
        header = '[ bondtypes ]' +'\n'+'; i        j          func      b0[b0_unit]         kb[kb_unit]'
        return str(header)

    def angle_header(self):
        header = '[ angletypes ]' + '\n'+'; i        j        k           func  theta0[theta0_unit]  ktheta[ktheta_unit]  ub0[ub0_unit]    kub[kub_unit]'
        return str(header)

    def dihedral_header(self):
        header = '[ dihedraltypes ]'+'\n'+'; i       j       k       l         func   phi0[phi0_unit]   kphi[kphi_unit] mult'
        return str(header)

    def improper_header(self):
        header = '[ dihedraltypes ] ;IMPROPERS'+'\n'+'; i       j       k       l         func   phi0[phi0_unit]   kphi[kphi_unit]'
        return(str(header))


#nearly identical class to MolecularData, but less functionality

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

