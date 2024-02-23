from correctedchargecalculation import calculate_charge
from masscalculation  import calculate_mass
import numpy as np
import matplotlib.pyplot as plt

pep_string = input("Enter your peptide chain: ")

def mass_spec(pep_string):

    charge = calculate_charge(pep_string)

    mass = calculate_mass(pep_string)

    print(type(mass))

    peaks = [0 for i in range(charge)]
    while charge>0:
        peaks[charge-1] = mass/charge
        charge = charge - 1

    print(peaks)
    return(peaks)
peaks = mass_spec(pep_string)

for i in range(len(peaks)):
    plt.axvline(peaks[i], ymin = 0, ymax = 0.8)
    plt.text(peaks[i], 0.81, peaks[i], ha='center', va='center')
plt.show()
        