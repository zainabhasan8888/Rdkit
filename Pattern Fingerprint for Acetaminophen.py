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

mol = Chem.MolFromSmiles('CC(=O)NC1=CC=C(O)C=C1')
img = Draw.MolToImage(mol)
img

# Generate Pattern Fingerprint for acetaminophen.
pf_bitvector = PatternFingerprint(mol)

# Convert the pattern fingerprint to NumPy array.
arr = np.zeros((0,), dtype=np.int8)
DataStructs.ConvertToNumpyArray(pf_bitvector, arr)

print("PatternFingerPrint for Acetaminophen is:")
print(arr)
