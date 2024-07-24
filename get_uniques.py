import logging as log
import pandas as pd
from parse_files import parse_cgen, parse_ff
import os
import click
import sys

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
# print(cgen.data['ffbonded.itp']['bonds'][0]) #DOES NOT CONTAIN UNITS 
# print(cgen.get_angles().keys()) #CONTAINS UNITS 


def get_uniques(ff, cgen, index_columns, param):
    try:
        merged = cgen.merge(ff, on=index_columns, how='left', indicator=True)
        
        # merge both arrays
        unique_entries = merged[merged['_merge'] == 'left_only']

        # Drop columns from ff_df that are not in index_columns
        columns_to_drop = ff.columns.difference(index_columns)
        columns_to_drop = [col for col in columns_to_drop if col in unique_entries.columns]

        unique_entries = unique_entries.drop(columns=columns_to_drop)
        unique_entries = unique_entries.drop(columns=['_merge', 'mult_y'], errors='ignore')  # Drop unnecessary columns

        # Remove columns and reset index
        unique_entries.reset_index(drop=True, inplace=True)
        return pd.DataFrame(unique_entries)
        
    #are there any other errors I can forsee?
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

    #screen for type of parameter based on headers from incoming dataframe 
    #there is probably a better way to do this...
    if all(col in df.columns for col in ['i', 'j', 'func', 'kb', 'b0']): #bonds
        format_template = "{:>8}  {:>8}  {:>8}   {:>8.8f}    {:>8.2f}" #bonds
        #REORDER AND RECAST
        split_line = unformatted_string.split()

        func = int(split_line[2])
        kb = float(split_line[3])
        b0 = float(split_line[4])
        
        reordered_line = split_line[:2] + [func, b0, kb]
        line  = format_template.format(*reordered_line)

    elif all(col in df.columns for col in ['i', 'j', 'k', 'func', 'ktheta', 'theta0','kub', 'r0']): #angles
        format_template = "{:>8}{:>8}{:>8}   {:>8}  {:>8.6}   {:>8.6f}  {:>8.4f}  {:>8.4f}"
        
        #REORDER AND RECAST
        split_line = unformatted_string.split() 
        
        ktheta = float(split_line[4])
        theta0 = float(split_line[5])
        r0 = float(split_line[7])
        kub = float(split_line[6])

        reordered_line = split_line[:4] + [theta0, ktheta, r0, kub]
        line = format_template.format(*reordered_line)

    elif all(col in df.columns for col in ['i', 'j', 'k', 'l', 'func', 'kphi', 'multi', 'phi0']): #dihedrals
        
        #REORDER AND RECAST
        format_template = "{:>8}{:>8}{:>8}{:>8}{:>8}    {:>8.8f}   {:>8.6f}  {:>8}"
        split_line = unformatted_string.split()

        kphi = float(split_line[5])
        multi = int(split_line[6])
        phi0 = float(split_line[7])

        reordered_line = split_line[:5] + [phi0, kphi, multi]
        line = format_template.format(*reordered_line)

    elif all(col in df.columns for col in ['i', 'j', 'k', 'l', 'func', 'kphi', 'phi0']): #impropers
        format_template = "{:>8}{:>8}{:>8}{:>8} {:>8}    {:>8.6f}   {:>8.6f}"
        split_line = unformatted_string.split()

        #REORDER AND RECAST
        kphi = float(split_line[4])
        func = split_line[5]
        phi0 = float(split_line[6])

        reordered_line = split_line[:4] + [func, phi0, kphi]
        # print(reordered_line)
        line = format_template.format(*reordered_line)
    
    return line


def update_charmm(df, outfile, header):
    # units are stored in classes.py

    #OLD WAY, more elegant
    header_columns = [col for col in df.columns if '_unit' not in col]
    units = [col for col in df.columns if '_unit' in col]

    file_exists = os.path.exists(outfile)
    with open(outfile, 'a' if file_exists else 'w') as f:
        f.write(header + '\n')
        # if not os.path.exists(outfile):
        #     f.write(';'+' '.join(header_columns) + '\n')

        if df is not None:
            df = df.fillna(float(0.0))
            for index, row in df.iterrows():
                unformatted_line = '   '.join(str(row[col]) if not isinstance(row[col], list) else ' '.join(map(str, row[col])) for col in header_columns)

                formatted_line = format_string(df, unformatted_line)
                f.write(formatted_line + '\n')

