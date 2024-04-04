import pandas as pd

def calculate_mass(peptide_string):
    # In case the input box is empty
    if peptide_string == "":
        return 0

    # Load the CSV file
    df = pd.read_csv("AminoAcidTable.csv")
    
    # Map each code to its Residue Mass
    code_to_mass = df.set_index('Code')['Residue Mass'].to_dict()
    
    # Sum the Residue Mass for each code in the string (and account for the extra missing water mass)
    total_mass = sum([code_to_mass[code] for code in peptide_string if code in code_to_mass]) + 18.01528
    
    return total_mass

def calculate_charge(peptide_string):
    # Load the CSV file
    df = pd.read_csv("AminoAcidTable.csv")
    
    # Map each code to its Charge
    code_to_charge = df.set_index('Code')['Charge'].to_dict()
    
    # Sum the Charge for each code in the string
    total_charge = sum([code_to_charge[code] for code in peptide_string if code in code_to_charge]) + 1
    
    return total_charge
