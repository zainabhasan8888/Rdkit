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

mol = Chem.MolFromSmiles('[H][N]([H])([H])[Pt](Cl)(Cl)[N]([H])([H])[H]')
img = Draw.MolToImage(mol)
img

maccs_bivector = MACCSkeys.GenMACCSKeys(mol)
print(maccs_bivector)
maccs_list= list(maccs_bivector)
print(maccs_list)
