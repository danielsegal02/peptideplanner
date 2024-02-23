import pandas as pd

aa_data = pd.read_csv('AminoAcidTable.csv')

def calculate_mass(pep_string):

    aa_array = [char for char in pep_string]

    filtered_aa = aa_data[(aa_data["Code"].isin(aa_array))]

    peptide_mass = filtered_aa["Residue Mass"].sum()

    peptide_mass = int(peptide_mass + 18.01528)

    print("The total mass of the peptide is: ", peptide_mass)

    return(int(peptide_mass))





