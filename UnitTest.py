import pandas as pd 
from ImageGenerator import generate_peptide_image, generate_mass_spec
from Calculations import calculate_mass, calculate_charge
# sources used for verifying correct info:
# FOR MASS:
# - https://pepcalc.com
# - https://web.expasy.org/peptide_mass/
# - The excel provided by the Medina Lab
# CHARGE: calculated manually in the UnitTest.py itself
# MASS SPECTROMETRY:
#
# NOTES:
# - mass function gives "monoisotopic mass"
# - O, K, R, o, k, r amino acids raise charge by 1, all others make no change 
#
# RECOMMENDED PEPTIDES:
#
#
#

# Amino Acids for testing the validity mass and charge calculations
peptides=["ARN","DARN","ARNDE","EEA"]
correct_masses=[359.19,474.22,603.26,347.13] # the monoisotpoic mass before adjusting for termini
c_terminus=["Amide","Acid","Amide","Acid"]
n_terminus=["Acetyl","Amine","Acetyl","Amine"]
# Amino Acids for testing validity of Mass Spec calculations
m_peptides=["KRWHWWRRHWVVW","WKWLKKWIK","KRWWKWWRR","RRWWRWVVW"]
m_c_terminus="Amide"
m_n_terminus="Amine"
expected_peaks=[[1009.6,673.4,505.3],[657.9,438.9,330.0],[743.9,496.3,372.5],[1428.8,714.9,476.9]]#2-4,1-3 for last

def generate_mass_spec(mass, charge):
    peaks = [0 for i in range(charge)]

    while charge > 0:
        peaks[charge - 1] = (mass + charge) / charge
        charge -= 1
    return peaks


for i in range(len(peptides)):
    
    # Manually calculates the charge without using calculations.py according to specific rules for the app to follow
    test_charge=0
    for element in peptides[i]:
        if element=="O" or element == "K" or element == "R":
            test_charge+=1
        if element=="o" or element == "k" or element == "r":
            test_charge+=1
    test_charge+=1
    
    if n_terminus[i] == "Acetyl":
        correct_masses[i] += 42.04
    if c_terminus[i] == "Amide":
        correct_masses[i] -= 0.986
    if n_terminus[i] == "Amine":
        test_charge += 1
    
    #Show what is calculated by the program alongside what expected values are 
    calc_mass=round(calculate_mass(peptides[i],n_terminus[i],c_terminus[i]),2)
    calc_charge=calculate_charge(peptides[i],n_terminus[i],c_terminus[i])
    print("Calculated Mass Of ",peptides[i],": ",calc_mass)
    print("Calculated Charge Of ",peptides[i],": ",calc_charge)
    print("Expected Mass Of ",peptides[i],":",correct_masses[i])
    print("Expected Charge Of ",peptides[i],":",test_charge)
    
    
for i in range(len(m_peptides)):
    mass=calculate_mass(m_peptides[i],m_n_terminus,m_c_terminus)
    charge=calculate_charge(m_peptides[i],m_n_terminus,m_c_terminus)    
    peaks=generate_mass_spec(mass,charge)
    if i==3:
        print("Calculated Mass Spec Peaks for ",m_peptides[i],": ",round(peaks[0],1),round(peaks[1],1),round(peaks[2],1))
    else:
        print("Calculated Mass Spec Peaks for ",m_peptides[i],": ",round(peaks[1],1),round(peaks[2],1),round(peaks[3],1))
        
    print("Expected Mass Spec Peaks for ",m_peptides[i],": ",expected_peaks[i][0],expected_peaks[i][1],expected_peaks[i][2])
    
    