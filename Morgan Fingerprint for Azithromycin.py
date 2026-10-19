import pandas as pd
import numpy as np
import warnings

import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, PandasTools, MACCSkeys, AtomPairs, rdFingerprintGenerator
from rdkit.Chem import DataStructs
from rdkit.Chem.rdmolops import PatternFingerprint
from rdkit.Avalon import pyAvalonTools
from rdkit.Chem.AtomPairs.Pairs import GetAtomPairFingerprintAsBitVect

mol = Chem.MolFromSmiles('CC[C@H]1OC(=O)[C@H](C)[C@@H](O[C@H]2C[C@@](C)(OC)[C@@H](O)[C@H](C)O2)[C@H](C)[C@@H](O[C@@H]2O[C@H](C)C[C@@H]([C@H]2O)N(C)C)[C@](C)(O)C[C@@H](C)CN(C)[C@H](C)[C@@H](O)[C@]1(C)O')
img = Draw.MolToImage(mol)
img

# Generate Morgan Fingerprint for Azithromycin.

mf_bitvector = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits= 2048)


# Convert the Fingerprint into a NumPy Array.

arr= np.zeros((0,), dtype= np.int8)
DataStructs.ConvertToNumpyArray (mf_bitvector, arr)

print("Morgan Fingerprint for Azithromycin is:")
print(arr)
