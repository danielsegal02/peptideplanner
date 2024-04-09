from correctedchargecalculation import sum_charge
from masscalculation  import sum_mass
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def mass_spec():

    charge = sum_charge()

    mass = sum_mass()

    peaks = [0 for i in range(charge)]

    while charge>0:

        peaks[charge-1] = (mass + charge)/charge

        charge = charge - 1

    print(peaks)

    return(peaks)

peaks = mass_spec()

data = []

for i in range(len(peaks)):

    plt.axvline(peaks[i], ymin = 0, ymax = 0.8)
    
    plt.text(peaks[i], 0.81, [round(peaks[i], 1),f'M$^{i+1}$$^+$'], fontsize = 11, ha='center', va='center', rotation = 30)

    data.append(
            { 'Peak': f'M + {i+1}',
              'Mass': peaks[i]}
    )    

dfdata = pd.DataFrame(data)

plt.xlim(0, max(peaks)+100)
plt.xlabel("m/z")

plt.ylabel("Intensity")

table_data = [dfdata.columns.tolist()] + dfdata.values.tolist()

table = plt.table(cellText=table_data, loc='bottom', cellLoc='center', bbox = [0, -0.3, 1, 0.2])

plt.subplots_adjust(bottom=0.2)

fig = plt.gcf()
fig.set_size_inches(20, 10.5)
fig.savefig('massspec.png')
