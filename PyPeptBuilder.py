from pyPept.sequence import Sequence
from pyPept.sequence import correct_pdb_atoms
from pyPept.molecule import Molecule
from pyPept.conformer import Conformer
from pyPept.conformer import SecStructPredictor

# RDKit modules
from rdkit import Chem
from rdkit.Chem import Draw

def generate_peptide_image(peptide_string):
    seq = Sequence(peptide_string)
    seq = correct_pdb_atoms(seq)    # Correct atom names in the sequence object
    mol = Molecule(seq, depiction='rdkit')  # Generate the RDKit object
    romol = mol.get_molecule(fmt='ROMol')   # Uses rdkit to get the molecule and draw it to a pdf
    image_path = 'peptide_image.png'
    Draw.MolToFile(romol, image_path, size=(600, 200))  #customize size to window size if possible
    return image_path


# Test Amino Acid input
# biln = "ac-D-T-H-F-E-I-A-am"

# Code to print the SMILES format
# print("The SMILES of the peptide is: {}".format(Chem.MolToSmiles(romol)))


# Does difference between HELM and BILN matter??
