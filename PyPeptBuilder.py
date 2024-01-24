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
    seq = correct_pdb_atoms(seq)
    mol = Molecule(seq, depiction='rdkit')
    romol = mol.get_molecule(fmt='ROMol')
    image_path = 'peptide_image.png'
    Draw.MolToFile(romol, image_path, size=(400, 400))
    return image_path
























# # Start the Sequence object
# biln = "ac-D-T-H-F-E-I-A-am"


# # difference between HELM and BILN?????

# ################################################
# ## With HELM: Call the converter to change from HELM to BILN
# #helm = "PEPTIDE1{[ac].D.T.H.F.E.I.A.[am]}$$$$V2.0"
# #b = Converter(helm=helm)
# #biln = b.get_biln()
# ################################################

# seq = Sequence(biln)
# # Correct atom names in the sequence object
# seq = correct_pdb_atoms(seq)

# # # Loop wit the included monomers
# # mm_list = seq.s_monomers
# # for i, monomer in enumerate(mm_list):
# #     mon = monomer['m_romol']

# # Generate the RDKit object
# mol = Molecule(seq, depiction='rdkit')
# romol = mol.get_molecule(fmt='ROMol') #uses rdkit to get the molecule and draw it to a pdf
# # print("The SMILES of the peptide is: {}".format(Chem.MolToSmiles(romol)))
# Draw.MolToFile(romol, 'peptide_attempt.png', size=(1200, 1200))
