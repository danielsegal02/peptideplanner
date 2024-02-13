import pandas as pd 

aadata = pd.read_csv('AminoAcidTable.csv')

aadict = aadata.to_dict('records')

print(aadata)