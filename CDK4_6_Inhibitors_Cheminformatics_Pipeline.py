import numpy as np
import pandas as pd
import rdkit
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, AllChem

drugs_data = {
    "Drug_name": ["Palbociclib", "Ribociclib", "Abemaciclib"],
    "SMILES": [
        "CC(=O)C1=C(C)C2=CN=C(NC3=NC=C(C=C3)N3CCNCC3)N=C2N(C2CCCC2)C1=O", # Palbociclib
        "CN(C)C(=O)C1=CC2=CN=C(NC3=CC=C(C=N3)N3CCNCC3)N=C2N1C1CCCC1", # Ribociclib
        "CCN1CCN(CC2=CC=C(NC3=NC=C(F)C(=N3)C3=CC(F)=C4N=C(C)N(C(C)C)C4=C3)N=C2)CC1" # Abemaciclib
    ]
}

df = pd.DataFrame(drugs_data)

def evaluate_lipinski_descriptors(smiles):
  mol = Chem.MolFromSmiles(smiles)
  mw = Descriptors.MolWt(mol)
  logp = Descriptors.MolLogP(mol)
  h_donors = Descriptors.NumHDonors(mol)
  h_acceptors = Descriptors.NumHAcceptors(mol)
  passes = (mw <= 500) and (logp <= 5) and (h_donors <= 5) and (h_acceptors <= 10)
  return pd.Series([mw, logp, h_donors, h_acceptors, passes])
  index= ['Molecular_Weight', 'LogP', 'H_donors', 'H_Acceptors', 'Passes_Lipinski']





def evaluate_lipinski_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    h_donors = Descriptors.NumHDonors(mol)
    h_acceptors = Descriptors.NumHAcceptors(mol)
    passes = (mw <= 500) and (logp <= 5) and (h_donors <= 5) and (h_acceptors <= 10)
    return pd.Series([mw, logp, h_donors, h_acceptors, passes], 
                     index=['Molecular_Weight', 'LogP', 'H_donors', 'H_Acceptors', 'Passes_Lipinski'])


lipinski_results = df["SMILES"].apply(evaluate_lipinski_descriptors)
df_combined = pd.concat([df, lipinski_results], axis=1)

print("---- Lipinski's rule screening results ----")
print(df_combined[['Drug_name', 'Molecular_Weight', 'Passes_Lipinski']])

print("---- Tanimoto similarity Matrix (Morgan Fingerprint)----")
mols = [Chem.MolFromSmiles(s) for s in df['SMILES']]
fps = [AllChem.GetMorganFingerprintAsBitVect(m, radius=2, nBits=2048) for m in mols]

similarity_matrix = []
for i in range(len(fps)):
    row = []
    for j in range(len(fps)):
        sim = DataStructs.TanimotoSimilarity(fps[i], fps[j])
        row.append(round(sim, 3))
    similarity_matrix.append(row)

df_sim = pd.DataFrame(similarity_matrix, index=df['Drug_name'], columns=df['Drug_name'])
print(df_sim)

