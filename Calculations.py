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

def calculate_reagent_mass(mol_equiv, reagent_MW, pepscale, resin_p_used):
    #what do we do if the input boxes are empty? we should account for that.
    resin_percent = resin_p_used * 0.01
    return (mol_equiv*reagent_MW*pepscale*resin_percent)

def calculate_reagent_vol(mol_equiv, reagent_MW, reagent_density, pepscale, resin_p_used, is_dry):
    rg = 1
    if is_dry == False:
        rg = reagent_density** -1 #maybe set density to 1 and/or make it unmodifyable wwhen dry is true
    #if dry is clicked reagent density will be 1 and if not it will be modifhable, so ahve that modified.
    resin_percent = resin_p_used * 0.01
    return (mol_equiv*reagent_MW*rg*pepscale*resin_percent) # might need to round up or down...

def calculate_solvent_vol(volFactor, pepscale, resin_p_used):
    resin_percent = resin_p_used * 0.01
    return (10*volFactor*pepscale*resin_percent) 
