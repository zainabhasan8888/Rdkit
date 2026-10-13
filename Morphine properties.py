from rdkit import Chem
from rdkit.Chem import Draw
m = Chem.MolFromSmiles('[H][C@@]12OC3=C(O)C=CC4=C3[C@@]11CCN(C)[C@]([H])(C4)[C@]1([H])C=C[C@@H]2O')
img = Draw.MolToImage(m)
img
# counting the number of atoms in a compound (Non-Hydrogen Atoms)
m.GetNumAtoms()
# Get symbol of atoms in a compound (Non_hydrogen atoms)
[atom.GetSymbol() for atom in m.GetAtoms()]
# show compound with hydrogen atoms.
morphine_with_hydrogen = Chem.AddHs(m)
img_2 = Draw.MolToImage(morphine_with_hydrogen)
img_2
# calculate molecular weight of a compound.
from rdkit.Chem import Descriptors
Descriptors.MolWt(m)
# Calculating Log P for a compound.
Descriptors.MolLogP(m)
# count the number of rotatable bonds in a compound.
Descriptors.NumRotatableBonds(m)
