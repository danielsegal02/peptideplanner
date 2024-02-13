import pandas as pd

aa_data = pd.read_csv('AminoAcidTable.csv')

peptide_string = input("Enter your peptide chain: ")

def calculate_charge(peptide_string):

    aa_array = [char for char in pep_string]

    filtered_aa = aa_data[(aa_data["Code"].isin(aa_array))]

    total_charge = filtered_aa["Charge"].sum()
    
    print("The total charge of the peptide is: ", total_charge)

print(aa_data)

