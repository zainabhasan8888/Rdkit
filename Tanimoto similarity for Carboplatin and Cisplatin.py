import pandas as pd
import numpy as np
import rdkit
import warnings
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, PandasTools, MACCSkeys, AtomPairs, rdFingerprintGenerator
from rdkit.Chem import DataStructs
from rdkit.Chem.rdmolops import PatternFingerprint
from rdkit.Avalon import pyAvalonTools
from rdkit.Chem.AtomPairs.Pairs import GetAtomPairFingerprintAsBitVect

# Draw the structure of Carboplatin

mol_1 = Chem.MolFromSmiles('[H][N]([H])([H])[Pt]1(OC(=O)C2(CCC2)C(=O)O1)[N]([H])([H])[H]')
img = Draw.MolToImage(mol_1)
img

# Draw the structure of Cisplatin.
mol_2 = Chem.MolFromSmiles('[H][N]([H])([H])[Pt](Cl)(Cl)[N]([H])([H])[H]')
img_2 = Draw.MolToImage(mol_2)
img_2

# Generate Morgan Fingerprint for both structures.
fp1= AllChem.GetMorganFingerprintAsBitVect(mol_1, 2, nBits=2048)
fp2= AllChem.GetMorganFingerprintAsBitVect(mol_2, 2, nBits=2048)

# Calculate Tanimoto Similarity metric.

similarity = DataStructs.TanimotoSimilarity(fp1,fp2)
print(f"Tanimoto similarity for the two compound is {similarity}")
