import pandas as pd

aa_data = pd.read_csv('AminoAcidTable.csv')

def calculate_charge(pep_string):

    aa_array = [char for char in pep_string]

    filtered_aa = aa_data[(aa_data["Code"].isin(aa_array))]

    total_charge = filtered_aa["Charge"].sum()
    
    print("The total charge of the peptide is: ", total_charge)

    return(int(total_charge))





