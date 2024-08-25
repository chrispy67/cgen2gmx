import logging as log
import pandas as pd
from src.parse_files import parse_cgen, parse_ff
import os
import click
import sys
from src import config

# How to loop through the dictionary and where the keys are addressable.
def iterate_nested_dict(nested_dict, parent_key=''):
    for key, value in nested_dict.items():
        full_key = f"{parent_key}.{key}" if parent_key else key
        if isinstance(value, dict):
            iterate_nested_dict(value, full_key)
        else:
            print(f"{full_key}")


##EXAMPLE
# print(cgen.data['ffbonded.itp']['bonds'][0]) #DOES NOT CONTAIN UNITS 
# print(cgen.get_angles().keys()) #CONTAINS UNITS 


# WIP, checking for backwards entries is off by default. 
def create_entries_set(df, index_columns):
    entries_set = set()
    for row in df[index_columns].values:
        direct_entry = tuple(row)
        reversed_entry = tuple(row[::-1])
        entries_set.add(direct_entry)
        entries_set.add(reversed_entry)
    return entries_set


def get_uniques(ff, cgen, index_columns, param):

    try:
        if config.reverse_entries: #WIP
            print('WIP--not checking for backwards entries')
            pass
        ##JUMBLES THE ORDER???
            ff_entries = create_entries_set(ff, index_columns) 
            cgen_entries = create_entries_set(cgen, index_columns)
            print('checking reverse entries...')
    
            # # screen for unique entries into df
            unique_entries = [entry for entry in cgen_entries if entry not in ff_entries]
            unique_entries_df = pd.DataFrame(unique_entries, columns=index_columns)
    
            # this merge checks for duplicates in 4 columns and merges them with EXISTING
            # parameters found in CGEN file
            unique_entries_df = unique_entries_df.merge(cgen, on=index_columns, how='left')
    
            columns_to_drop = ff.columns.difference(index_columns)
            unique_entries_df = unique_entries_df.drop(columns=columns_to_drop, errors='ignore')
            unique_entries_df.reset_index(drop=True, inplace=True)
            unique_entries_df = unique_entries_df.dropna()
    
            unique_entries_df.reset_index(drop=True, inplace=True)  
            return unique_entries_df
        
        # Traditional matches
        # merging two dfs
        merged = cgen.merge(ff, on=index_columns, how='left', indicator=True)
        
        # merge both arrays
        unique_entries = merged[merged['_merge'] == 'left_only']

        # Drop columns from ff_df that are not in index_columns
        columns_to_drop = ff.columns.difference(index_columns)
        columns_to_drop = [col for col in columns_to_drop if col in unique_entries.columns]
        unique_entries_df = unique_entries.drop(columns=columns_to_drop)
        unique_entries_df = unique_entries.drop(columns=['_merge', 'mult_y'], errors='ignore')  # Drop unnecessary columns

        # Remove columns and reset index
        unique_entries_df.reset_index(drop=True, inplace=True)
        return unique_entries_df
        
    except AttributeError as e:
        print(f'AttributeError: {e}')
        print(f'cannot merge empty array, possibly there are no entires for {param}?')
    except KeyError as e:
        if len(ff.columns) == 0:
            print(f'KeyError: {e}')
            print(f'length of merged array for {param} is 0, no unique entries detected?')
        else:
            print(f'KeyError: {e}')
            print(f'ff_df columns: {ff.columns}')
            print(f'index_columns: {index_columns}')


# print(unique_bonds['kb_unit'].iloc[0])

def format_string(df, unformatted_string):

    # Screen for type of parameter based on headers from incoming dataframe 
    # There may be a more flexible way to do this
    if all(col in df.columns for col in ['i', 'j', 'func', 'kb', 'b0']): #bonds
        format_template = "{:>8}  {:>8}  {:>8}   {:>8.5f}        {:>8.2f}" 
        
        #REORDER AND RECAST
        split_line = unformatted_string.split()
        func = int(split_line[2])
        kb = float(split_line[3])
        b0 = float(split_line[4])
        
        reordered_line = split_line[:2] + [func, b0, kb]
        line  = format_template.format(*reordered_line)

    elif all(col in df.columns for col in ['i', 'j', 'k', 'func', 'ktheta', 'theta0','kub', 'ub0']): #angles
        format_template = "{:>8}{:>8}{:>8}   {:>8}    {:>8.3f}   {:>8.5f}  {:>6.4f}  {:>8.4f}"
        
        #REORDER AND RECAST
        split_line = unformatted_string.split() 
        ktheta = float(split_line[4])
        theta0 = float(split_line[5])
        r0 = float(split_line[7])
        kub = float(split_line[6])

        reordered_line = split_line[:4] + [theta0, ktheta, r0, kub]
        line = format_template.format(*reordered_line)

    elif all(col in df.columns for col in ['i', 'j', 'k', 'l', 'func', 'kphi', 'multi', 'phi0']): #dihedrals
        format_template = "{:>8}{:>8}{:>8}{:>8}{:>8}    {:>8.4f}   {:>8.5f}     {:>1}"

        #REORDER AND RECAST
        split_line = unformatted_string.split()
        kphi = float(split_line[5])
        multi = int(split_line[6])
        phi0 = float(split_line[7])

        reordered_line = split_line[:5] + [phi0, kphi, multi]
        line = format_template.format(*reordered_line)

    elif all(col in df.columns for col in ['i', 'j', 'k', 'l', 'func', 'kphi', 'phi0']): #impropers
        format_template = "{:>8}{:>8}{:>8}{:>8}{:>8}    {:>8.4f}   {:>8.5f}"

        #REORDER AND RECAST
        split_line = unformatted_string.split()
        kphi = float(split_line[4])
        func = split_line[5]
        phi0 = float(split_line[6])

        reordered_line = split_line[:4] + [func, phi0, kphi]
        line = format_template.format(*reordered_line)
    
    return line


def update_charmm(df, outfile, header, param):

    # Create a mapping between unit columns and their respective units
    # Necessary when reordering of the columns in format_string()
    unit_mapping = {}
    unit_cols = [col for col in df.columns if '_unit' in col]
    try:
        for col in unit_cols:
            unit_mapping[col.replace('_unit', '')] = df[col].iloc[0] 
    
    # to circumvent an empty array from get_uniques() which doesn't have df[col].iloc[0]
    except IndexError as e:
        unit_mapping = {}
        print(e, ' *from update_charmm() from get_uniques.py*')
        print(f' - It is likely that no unique {param} were found in --itp input')
    
    # This redefine is crucial here
    formatted_header = header

    #match each '*_unit' key to its matching '*' key
    for key, unit in unit_mapping.items():
        formatted_header = formatted_header.replace(f'[{key}_unit]', f'[{unit}]')

    file_exists = os.path.exists(outfile)
    with open(outfile, 'a' if file_exists else 'w') as f:
        f.write(formatted_header + '\n')

        if df is not None:
            df = df.fillna(float(0.0))
            header_columns = [col for col in df.columns if '_unit' not in col]
            for index, row in df.iterrows():
                unformatted_line = '   '.join(str(row[col]) if not isinstance(row[col], list) else ' '.join(map(str, row[col])) for col in header_columns)

                formatted_line = format_string(df, unformatted_line)
                f.write(formatted_line + '\n')


