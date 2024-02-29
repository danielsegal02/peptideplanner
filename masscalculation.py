import pandas as pd

from collections import Counter

def char_occ(pep_string):
    char_count = Counter(pep_string)
    return dict(char_count)

def filter_dataframe_mass(aa_data, pep_string):

    aa_array = [char for char in pep_string]
    
    filtered_df = aa_data[aa_data['Code'].isin(aa_array)]
    
    result_dict = dict(zip(filtered_df['Code'], filtered_df['Residue Mass']))
    
    return result_dict

def multiply_dictionaries(dict1, dict2):

    result_dict = {key: value * dict2[key] for key, value in dict1.items()}

    return result_dict

def calculate_mass():

    aa_data = pd.read_csv('AminoAcidTable.csv')

    aa_chain = input("Enter your amino acid chain: ")

    char_count = char_occ(aa_chain)

    filt_mass = filter_dataframe_mass(aa_data, char_count.keys())

    final_dict = multiply_dictionaries(char_count, filt_mass)

    return final_dict

def sum_mass():
    result = calculate_mass()
    result_mass = sum(result.values()) + 18.01
    return result_mass