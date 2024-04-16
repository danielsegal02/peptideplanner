import sys,os

import tkinter as tk

IPYNB_FILENAME = 'secondarystructure.ipynb'
CONFIG_FILENAME = '.config_ipynb'

def generate_sec_struct(amino):
    arg = "--str_param "+amino
    with open(CONFIG_FILENAME,'w') as f:
        f.write(arg)
   # os.system('jupyter nbconvert --to notebook --execute test_argv.ipynb')
   # os.system('jupyter nbconvert --execute {:s} --to notebook'.format(IPYNB_FILENAME))
    os.system('jupyter notebook secondarystructure.ipynb')
    return None



def retrieve_input(entry):
    input_value = entry.get()
    print(input_value)
    #generate_sec_struct(input_value)




    