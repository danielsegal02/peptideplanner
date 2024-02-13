# RDKit modules
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdChemReactions

def combine_smiles(amino_acids_smiles):
    """
    Combines individual amino acid SMILES into a single peptide SMILES string using RDKit.
    
    Args:
        amino_acids_smiles (list of str): List of SMILES strings for the amino acids.
    
    Returns:
        str: A single SMILES string representing the peptide.
    """
    amino_acids_smiles = [                    # currently hard-coded, but will call a function that creates this list based on the correlating smiles in the csv
    "N[C@@H](C)C(=O)O",  # Alanine
    "N[C@@H](CCCNC(=N)N)C(=O)O",  # Arginine
    "N[C@@H](CC(=O)N)C(=O)O"  # Asparagine
    ]
    
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


def generate_peptide_image(pep_str):
    final_smiles = combine_smiles(pep_str)
    mol = Chem.MolFromSmiles(final_smiles)    # currently hard-coded, but later make it turn the peptide_string to the mapped smiles, then use it to make the image
    image_path = 'peptide_image.png'
    Draw.MolToFile(mol, image_path, size=(600, 200))  #customize size to window size if possible
    return image_path


# Test Amino Acid input
# biln = "ac-D-T-H-F-E-I-A-am"
