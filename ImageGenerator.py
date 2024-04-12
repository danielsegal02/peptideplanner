import pandas as pd
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Draw, rdChemReactions, rdDepictor

## Functions to retrieve the user input and put it into the correct format

def get_smiles_from_code(code_string):
    # Load the CSV file
    df = pd.read_csv("AminoAcidTable.csv")
    
    # Map each code to its SMILES representation
    code_to_smiles = df.set_index('Code')['SMILES'].to_dict()
    
    # Convert the code string to a list of corresponding SMILES
    smiles_list = [code_to_smiles[code] for code in code_string if code in code_to_smiles]
    
    return smiles_list

def modify_n_terminus(pep_smiles_lst):
    """
    Modifies the N-terminus of the first amino acid in a list of amino acid SMILES strings by reacting it with an acetyl group.
    
    Args:
        pep_smiles_lst (list of str): List of SMILES strings for the amino acids.
        
    Returns:
        list of str: The original list with the N-terminus of the first amino acid modified.
    """
    # Acetyl group SMILES for N-terminus modification
    acetyl_smiles = 'CC(=O)O'
    amine_mol = Chem.MolFromSmiles(acetyl_smiles)

    # Take the first amino acid SMILES from the list
    last_aa_smiles = pep_smiles_lst[0]
    last_aa_mol = Chem.MolFromSmiles(last_aa_smiles)

    # Define a reaction for adding the acetyl group to the N-terminus
    # The reaction targets an amine group and attaches the acetyl group
    amidation_rxn = rdChemReactions.ReactionFromSmarts('[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]')
    product_set = amidation_rxn.RunReactants((amine_mol, last_aa_mol))

    if product_set:
        # Update the first amino acid with the acetylated product
        modified_last_aa_mol = product_set[0][0]
        modified_last_aa_smiles = Chem.MolToSmiles(modified_last_aa_mol)
        # Replace the first amino acid SMILES in the list with the modified one
        pep_smiles_lst[0] = modified_last_aa_smiles

    return pep_smiles_lst

def modify_c_terminus(final_smiles):
    """
    Modifies the C-terminus of the last amino acid in a list of amino acid SMILES strings by replacing the carboxylic acid group with an amide group.
    
    Args:
        pep_smiles_lst (list of str): List of SMILES strings for the amino acids.
        
    Returns:
        list of str: The original list with the C-terminus of the last amino acid modified.
    """
    # Define the carboxylic acid and amide groups
    carboxylic_acid = "C(=O)O"
    amide_group = "C(=O)N"

    # Reverse the strings
    final_smiles_reversed = final_smiles[::-1]
    carboxylic_acid_reversed = carboxylic_acid[::-1]
    amide_group_reversed = amide_group[::-1]

    # Replace the first occurrence of the reversed carboxylic acid with the reversed amide group
    modified_smiles_reversed = final_smiles_reversed.replace(carboxylic_acid_reversed, amide_group_reversed, 1)

    # Reverse the result back to get the original orientation with the last occurrence replaced
    modified_final_smiles = modified_smiles_reversed[::-1]

    return modified_final_smiles



def combine_smiles(amino_acids_smiles):
    """
    Combines individual amino acid SMILES into a single peptide SMILES string using RDKit.
    
    Args:
        amino_acids_smiles (list of str): List of SMILES strings for the amino acids.
    
    Returns:
        str: A single SMILES string representing the peptide.
    """
    # Initialize an empty molecule that will serve as the base for the peptide chain
    peptide_mol = Chem.RWMol()

    # Loop through each amino acid SMILES string
    for i, smiles in enumerate(amino_acids_smiles):
        # Convert the SMILES string to a RDKit molecule
        aa_mol = Chem.MolFromSmiles(smiles)
        
        if i == 0:
            # For the first amino acid, just add it to the peptide molecule
            peptide_mol = Chem.RWMol(aa_mol)
        else:
            # For subsequent amino acids, perform a peptide bond formation reaction
            # Define a generic peptide coupling reaction
            # The reaction removes a water molecule to form the peptide bond
            peptide_rxn = rdChemReactions.ReactionFromSmarts('[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]')
            # Combine the current peptide molecule with the new amino acid molecule
            product_set = peptide_rxn.RunReactants((peptide_mol, aa_mol))
            
            if product_set:
                # Update the peptide molecule with the first product
                peptide_mol = product_set[0][0]
    
    # Convert the RDKit molecule back to a SMILES string
    final_peptide_smiles = Chem.MolToSmiles(peptide_mol)
    
    return final_peptide_smiles


## Functions to generate the actual images for the GUI

def generate_peptide_image(pep_str, n_terminus, c_terminus):
    # Turn input into equivalent SMILES codes
    pep_smiles_lst = get_smiles_from_code(pep_str)
    # Modify the first amino acid if the Acetyl option is chosen
    if n_terminus == "Acetyl":
        pep_smiles_lst = modify_n_terminus(pep_smiles_lst)
    # Combine all the amino acids smiles into one long SMILES code
    final_smiles = combine_smiles(pep_smiles_lst)
    # Modifies the c-terminus after the smiles has been combined
    if c_terminus == "Amide":
        final_smiles = modify_c_terminus(final_smiles)
    # Creates a molecule object from the final, combined, modified SMILES
    mol = Chem.MolFromSmiles(final_smiles)
    # Sets the preferences to get the structure to try and have a straight backbone
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(mol)
    # Produces the image from the molecule object in the designated path (saves within the Images folder)
    image_path = 'Images/peptide_image.png'
    Draw.MolToFile(mol, image_path, size=(800, 300))  # try to customize size to window size if possible
    return image_path

def generate_mass_spec(mass, charge):
    # Increase figure size for better layout of plot and table
    plt.figure(figsize=(12, 4))  # Increase width to accommodate table on the right
    peaks = [0 for i in range(charge)]
    
    while charge > 0:
        peaks[charge - 1] = (mass + charge) / charge
        charge -= 1
    
    # Plot each peak
    for i in range(len(peaks)):
        plt.axvline(peaks[i], ymin=0, ymax=0.8)
        plt.text(peaks[i], 0.81, [round(peaks[i], 1), f'M$^{i+1}$$^+$'], fontsize=11, ha='center', va='center', rotation=30)

    # Adjust the plot area to make room for the table on the right
    plt.subplots_adjust(left=0.1, right=0.6)
    plt.xlim(0, max(peaks)+100)
    plt.xlabel("m/z")
    plt.ylabel("Intensity")
    
    # Data for table
    data = [{'Peak': f'M + {i+1}', 'Mass': peaks[i]} for i in range(len(peaks))]
    dfdata = pd.DataFrame(data)

    # Create the table on the right side of the plot
    table_data = [dfdata.columns.tolist()] + dfdata.values.tolist()
    final_table = plt.table(cellText=table_data, loc='right', cellLoc='center', bbox=[1.1, 0.0, 0.4, 1])
    final_table.auto_set_font_size(False)
    final_table.set_fontsize(7.5)

    # Save the finished image in the designated path (saves within the Images Folder)
    file_path = "Images/MassSpec.png"
    plt.savefig(file_path, bbox_inches="tight")
    plt.close()  # Close the plot to free up memory
    return file_path
