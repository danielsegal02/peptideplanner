import pandas as pd 
from ImageGenerator import generate_peptide_image, generate_mass_spec
from Calculations import calculate_mass, calculate_charge
# sources used for verifying correct info:
# - https://pepcalc.com
# - https://web.expasy.org/peptide_mass/
# - The excel provided by the Medina Lab
# NOTES on thinks I noticed while playing with the functions:
# - mass function gives "monoisotopic mass"
# - charge function may be incorrect for charge going up or down seems to only give 0 or 1
peptides=["ARN","DARN","ARNDE"]
correct_masses1=[359.19,474.22,603.26]#the monoisotpoic mass
correct_masses2=[359.38,474.47,603.58]#the molecular weight/avg mass (is slightly higher generally)
correct_charges=[1, 0, -1]
correct_tests=0
for i in range(len(peptides)):
    calc_mass=round(calculate_mass(peptides[i]),2)
    calc_charge=calculate_charge(peptides[i])
    print("Calculated Mass Of ",peptides[i],": ",calc_mass)
    print("Calculated Charge Of ",peptides[i],": ",calc_charge)
    print("Expected Mass Of ",peptides[i],":",correct_masses1[i])
    print("Expected Charge Of ",peptides[i],":",correct_charges[i])
    if calc_mass==correct_masses1[i] and calc_charge==correct_charges[i]:
        print("Calculations Correct")
        correct_tests+=1
    else:
        print("Calcualtions did not match actual within margin of error.")
    print("Tests Passesd: ",correct_tests)
    if correct_tests==len(peptides):
        print("All tests ran correctly")

