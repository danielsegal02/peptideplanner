from correctedchargecalculation import sum_charge
from masscalculation  import sum_mass
import numpy as np
import matplotlib.pyplot as plt

def mass_spec():

    charge = sum_charge()

    mass = sum_mass()

    peaks = [0 for i in range(charge)]

    while charge>0:

        peaks[charge-1] = mass/charge

        charge = charge - 1

    print(peaks)

    return(peaks)

peaks = mass_spec()

for i in range(len(peaks)):

    plt.axvline(peaks[i], ymin = 0, ymax = 0.8)
    
    plt.text(peaks[i], 0.81, round(peaks[i], 1), ha='center', va='center')

plt.xlabel("m/z")

plt.ylabel("Intensity")

plt.show()
        