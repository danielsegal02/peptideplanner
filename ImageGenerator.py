# RDKit modules
from rdkit import Chem
from rdkit.Chem import Draw

def generate_peptide_image(peptide_string):
    mol = Chem.MolFromSmiles("N[C@H](Cc1c(F)c(F)c(F)c(F)c1F)C(=O)O")    # currently hard-coded, but later make it turn the peptide_string to the mapped smiles, then use it to make the image
    image_path = 'peptide_image.png'
    Draw.MolToFile(mol, image_path, size=(600, 200))  #customize size to window size if possible
    return image_path


# Test Amino Acid input
# biln = "ac-D-T-H-F-E-I-A-am"
