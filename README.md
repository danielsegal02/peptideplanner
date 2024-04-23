# Peptide Planner Application


## File Descriptions:
1. **MainWindow.py**
    - This is where all the code for the visual GUI is. The button, labels, and stlying and function calls all happen here.
2. **ImageGenerator.py** 
    - This is where the chemical structure png image, and the mass spec png image for the peptide is generated. It pull the SMILES from the csv file and manipulates them to then use in RdKit functions which creates the images. All the functions to generate the image are in this file, but they are called in the MainWindow.py file.
3. **Calculations.py**
    - This is where the mass, net charge, and all the reagent tab values of the peptide are calculated. All the necessary functions to generate the calculations are in this file, but the functions are called in the MainWindow.py file.
4. **AminoAcidTable.csv**
    - This is the "database" of all the amino acids that can be used as inputs for the application.
5. **Images Folder**
    - This is the folder that stores all the generated images, as well as the group picture in the credit tab.
6. **UnitTest.py**
    - This is a unit tester file that checks the accuracy of our functions. When run, it prints the outputs of our function calls and compares them to accurate numbers of what they should be.
7. **UnitTester/Sources**
    - Extra sources to reference for the unit tests / checking the accuracy of the results.

## Setup:
1. Install the following dependencies:
    - Python (3.9 but better versions MIGHT work)
    - pip install:
        - pandas
        - matplotlib
        - rdkit
2. Clone the repo
    - link: https://github.com/danielsegal02/peptideplanner.git
 

## Running the application:
1. Open a terminal (i.e. powershell on windows)
2. cd into the directory that has the repo cloned into it.
3. run this command in the terminal to start the application:
    python MainWindow.pya
