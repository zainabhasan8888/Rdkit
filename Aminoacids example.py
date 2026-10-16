from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit import DataStructs
import numpy as np

# Convert list of smiles (Glycine, Phenylalanin, Histidine, Cysteine) to list of mols.
smiles_list = ['N[C@@H](CS)C(O)=O','C(C(=O)O)N', 'N[C@@H](CC1=CC=CC=C1)C(O)=O', 'N[C@@H](CC1=CNC=N1)C(O)=O']
mol_list = []
for smiles in smiles_list:
  mol= Chem.MolFromSmiles(smiles)
  mol_list.append(mol)

img = Draw.MolsToGridImage(mol_list)
img

# Check which one of these aminoacids has sulfur atom in it.
pattern = Chem.MolFromSmiles('S')
for mol in mol_list:
  print(mol.HasSubstructMatch(pattern))

# Check which one of these aminoacids has carboxyl group in it.
pattern = Chem.MolFromSmiles('C(=O)O')
for mol in mol_list:
  print(mol.HasSubstructMatch(pattern))

# Check which one of these aminoacids has a ring.
pattern = Chem.MolFromSmarts('[r]')
for mol in mol_list:
  print(mol.HasSubstructMatch(pattern))

# Check which one of these aminoacids has a five-member ring.
pattern = Chem.MolFromSmarts('[r5]')
for mol in mol_list:
  print(mol.HasSubstructMatch(pattern))
