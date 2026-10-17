import rdkit
from rdkit import Chem, RDLogger
from rdkit.Chem import Draw, AllChem, PandasTools, BRICS, MACCSkeys, Descriptors
from rdkit.Chem.Draw import rdMolDraw2D, SimilarityMaps
from rdkit.DataStructs.cDataStructs import TanimotoSimilarity, DiceSimilarity
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem import Descriptors3D

import pandas as pd
import numpy as np

mol = Chem.MolFromSmiles('O=CN1CCCC1')
img = Draw.MolToImage(mol)
img

# save mol as PNG.
Draw.MolToFile(mol, 'mol.png')
