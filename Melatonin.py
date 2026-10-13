from rdkit import Chem
from rdkit.Chem import Draw
m = Chem.MolFromSmiles('CC(=O)NCCC1=CNc2c1cc(OC)cc2')
img = Draw.MolToImage(m)
img
