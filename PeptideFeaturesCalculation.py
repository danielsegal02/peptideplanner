import pyteomics as pyt
from pyteomics import mass
from PeptideFeaturesClass import PeptideFeatures

aminochain = input("Enter amino acid chain: ")

pH = int(input("Enter pH of your solution: "))

pep = PeptideFeatures(aminochain, pH)

molmass = pep.calc_mass()

print(aminochain, "has a mass of ",molmass, " g/mol")

molcharge = pep.calc_charge()

print(aminochain, "has a charge of ",molcharge)
