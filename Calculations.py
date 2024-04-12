import pandas as pd

def calculate_mass(peptide_string, n_terminus, c_terminus):
    # In case the input box is empty
    if peptide_string == "":
        return 0

    # Load the CSV file
    df = pd.read_csv("AminoAcidTable.csv")
    
    # Map each code to its Residue Mass
    code_to_mass = df.set_index('Code')['Residue Mass'].to_dict()
    
    # Sum the Residue Mass for each code in the string (and account for the extra missing water mass)
    total_mass = sum([code_to_mass[code] for code in peptide_string if code in code_to_mass]) + 18.01528

    if n_terminus == "Acetyl":
        total_mass += 42.04
    if c_terminus == "Amide":
        total_mass -= 0.986

    return total_mass

def calculate_charge(peptide_string, n_terminus, c_terminus):
    # Load the CSV file
    df = pd.read_csv("AminoAcidTable.csv")
    
    # Map each code to its Charge
    code_to_charge = df.set_index('Code')['Charge'].to_dict()
    
    # Sum the Charge for each code in the string
    total_charge = sum([code_to_charge[code] for code in peptide_string if code in code_to_charge]) + 1
    
    if n_terminus == "Amine":
        total_charge += 1
    
    return total_charge

def calculate_reagent_mass(pep_scale, percent_resin_used, reagent_MW, reagent_equiv):
    # Perform the calculation based on the formula
    reagent_mass = pep_scale * percent_resin_used * reagent_MW * reagent_equiv

    # Round the result to 1 decimal place
    reagent_mass_rounded = round(reagent_mass, 1)

    return reagent_mass_rounded

def calculate_reagent_volume(pep_scale, percent_resin_used, reagent_MW, reagent_equiv, reagent_density):
    # Perform the calculation based on the formula
    modified_reagent_density = 1/reagent_density 
    reagent_volume = pep_scale * percent_resin_used * reagent_MW * reagent_equiv * modified_reagent_density

    # Round the result to no decimal place
    reagent_volume_rounded = round(reagent_volume, 0)
    
    # Cast to an int to remove the decimal from the final value
    reagent_volume_rounded = int(reagent_volume_rounded)

    return reagent_volume_rounded

def calculate_solvent_volume(pep_scale, percent_resin_used, solvent_factor):
    # Perform the calculation based on the formula
    modified_solvent_factor = solvent_factor * 10
    solvent_volume = pep_scale * percent_resin_used * modified_solvent_factor

    # Round the result to 1 decimal place
    solvent_volume_rounded = round(solvent_volume, 1)

    return solvent_volume_rounded
