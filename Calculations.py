from pyteomics import mass,electrochem

def format_peptide_sequence(input_string):
    return input_string.replace('-', '')

def calculate_mass_and_charge(peptide):
    # Calculate the monoisotopic mass
    peptide_mass = mass.calculate_mass(sequence=peptide)
    
    # Calculate the charge.
    peptide_charge = electrochem.charge(peptide, 7)

    return peptide_mass, peptide_charge


# DOES NOT WORK FOR AMINO ACIDS LIKE (am, ac, etc...)


# example_peptide_sequence_input = "ACDEFGHIK"  #
