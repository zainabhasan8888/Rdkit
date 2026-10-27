import rdkit
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors

# Draw the structure of Palbociclib.

mol = Chem.MolFromSmiles('CC(=O)C1=C(C)C2=CN=C(NC3=NC=C(C=C3)N3CCNCC3)N=C2N(C2CCCC2)C1=O')
img = Draw.MolToImage(mol)
img

# Calculate Lipinski parameters. 
mw = Descriptors.MolWt(mol)
logp = Descriptors.MolLogP(mol)
h_donors = Descriptors.NumHDonors(mol)
h_acceptors = Descriptors.NumHAcceptors(mol)

print(f"Molecular weight: {mw} g/mol")
print (f"LogP: {logp}")
print(f"Number of hydrogen donors: {h_donors}")
print (f"H-Bond acceptors: {h_acceptors}")

# Check if Palbociclib passes the Lipinski's rule of five successfully.
pass_lipinski = (mw <= 500) and (logp <= 5) and (h_donors <= 5) and (h_acceptors <= 10)
print (f" The drug passes the Lipinski's rule of Five; {pass_lipinski}")
