# peptideplanner


Currently, there are four files:
1. MainWindow.py
    - This is where all the code for the visual GUI is. The button, labels, and stlying and function calls all happen here.
2. PyPeptBuilder.py 
    - This is where the PyPept library is accessed. All the necessary functions to generate the chemical structure png are in this file, but they are called in the MainWindow.py file.
3. Calculations.py
    - This is where the Pyteomics library is accessed. All the necessary functions to generate the mass and net charge calculations are in this file, but those are also called in the MainWindow.py file. There is also a function in here that formats the peptide string to be accepted by the pyteomics library.
4. AminoAcids.py
    - This was a start of a "Database" of Amino Acids. It was not needed for this code currently, but it might be useful to build off of later.


Setup:
- Install python if you don't already have it downloaded. I believe it should be python 3.9.
- Install PyPept:
    - pip install git+https://github.com/Boehringer-Ingelheim/pyPept.git
- Install Pyteomics:
    - pip install pyteomics
 

To run the application:
1. Open a terminal (i.e. powershell on windows)
2. cd into the directory that has the repo cloned into it.
3. run this command in the terminal to start the application:
    python MainWindow.py