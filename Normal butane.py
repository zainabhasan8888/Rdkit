from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit import DataStructs
import numpy as np

# Draw Normal Butane.
mol = Chem.MolFromSmiles('CCCC')
img = Draw.MolToImage(mol)
img

#Convert Image of the molecule back to smiles.
smiles = Chem.MolToSmiles(mol)
smiles

# Count the molecular weight of the molecule.
mw = Descriptors.MolWt(mol)
mw
