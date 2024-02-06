# Peptide Planner


Currently, there are four files:
1. MainWindow.py
    - This is where all the code for the visual GUI is. The button, labels, and stlying and function calls all happen here.
2. ImageGenerator.py 
    - This is where the chemical structure png image for the peptide is generated. It applies the rdKit modules to make the image by using the SMILES format of the given amino acid input. All the functions to generate the image are in this file, but they are called in the MainWindow.py file.
3. Calculations.py
    - This is where the mass of the peptide is calculated without any outside libraries. Soon the Net charge will be calculated here as well. All the necessary functions to generate the mass and net charge calculations are in this file, but those are called in the MainWindow.py file.
4. AminoAcids.py
    - This was a start of a "Database" of Amino Acids. It was not needed for this code currently, but it might be useful to build off of later.


Setup:
- Install python if you don't already have it downloaded. I believe it should be python 3.9 or better. (It might not matter anymore)
- Must of rdKit installed
    - I have it since I pip installed PyPept, but rdKit can be installed on it's own as well.
    - I will look into uninstalling PyPept and Pyteomics and checking what exactly is needed for it to run.
 

To run the application:
1. Open a terminal (i.e. powershell on windows)
2. cd into the directory that has the repo cloned into it.
3. run this command in the terminal to start the application:
    python MainWindow.py