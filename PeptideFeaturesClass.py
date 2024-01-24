import pyteomics as pyt
from pyteomics import mass
from pyteomics import electrochem

class PeptideFeatures:

    def __init__(self, chain, pH):
        self.chain = chain
        self.pH = pH
        
    def calc_mass(self):
        return mass.calculate_mass(sequence = self.chain)
    
    def calc_charge(self):
        return electrochem.charge(self.chain, self.pH)  



