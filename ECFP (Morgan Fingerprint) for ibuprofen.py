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

# Draw the structure of Ibuprofen from SMILES.
mol = Chem.MolFromSmiles('CC(C)CC1=CC=C(C=C1)C(C)C(=O)O')
img = Draw.MolToImage(mol)
img

mol_fp = rdFingerprintGenerator.GetMorganGenerator(radius= 2, fpSize=2048)

ecfp_mol= mol_fp.GetFingerprint(mol)

ecfp_mol_Numpy = np.zeros((0,), dtype=np.int8)
DataStructs.ConvertToNumpyArray(ecfp_mol, ecfp_mol_Numpy)
print(ecfp_mol_Numpy)
