def format_peptide_sequence(input_string):
    return input_string.replace('-', '')

def calculate_mass(peptide_string):
    ###
    # Process to calculate Mass:
        # add the mass of each amino acid in the chain
        # multiply the mass of water by the number of bonds (n-1 where n = the number of amino acids)
        # subtract multiplied value from the summed mass of all the amino acids
    ###
    mass_of_water = 18.01528    # fixed mass of water
    mass = 0    # initialize mass variable
    number_of_AA = 0    # initialize number of amino acids variable
   
    for amino_acid in peptide_string:
        #mass += amino_acid.mass # sum the mass of each amino acid from the "Database"  
        number_of_AA += 1   # increment the number of amino acids with each iteration

    mass += 89 + 174 + 132  # currently hard-coded just to test the function
    bond_mass = mass_of_water * (number_of_AA-1)    # find the mass of the bonds (mass of water - number of bonds)
    mass = mass - bond_mass # calculate total mass
    return mass

# example_peptide_sequence_input = "ACDEFGHIK"
