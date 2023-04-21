import pandas as pd

df = pd.read_csv("data/Chembl_drug_targets_20220517.tsv", sep = '\t')
print(df.head())

df = df[df.Organism == "Homo sapiens"]
df = df[(df.Type == "SINGLE PROTEIN") | (df.Type == "PROTEIN FAMILY")]
df = df[["ChEMBL ID", "UniProt Accessions"]]
df = df.rename(columns = {"ChEMBL ID": "ChEMBLID", "UniProt Accessions": "UniProt"})

# explode protein column
df = (df.set_index(df.columns.drop('UniProt',1).tolist())
                .UniProt.str.split('|', expand=True)
                .stack()
                .reset_index()
                .rename(columns={0:'UniProt'})
                .loc[:, df.columns]
                )

# extract list of protein IDs
prot_df = df[["UniProt"]].drop_duplicates()
prot_df.to_csv("data/Chembl_drug_targets_UniProtIDs.csv", index = False, header = False)

# mapping to EnsemblId: https://www.uniprot.org/uploadlists/
# 2. Select options: From UniProtKB AC/ID to Ensembl

# Mapping message: 4,164 out of 4,200 identifiers from UniProtKB AC/ID were successfully mapped to 4,583 Ensembl IDs.

prot_df = pd.read_csv("data/Chembl_drug_targets_UniProtIDs_EnsemblIDs.tsv", sep = '\t')
prot_df = prot_df.rename(columns = {"From": "UniProt", "To": "EnsemblID"})
df = df.merge(prot_df, on = "UniProt")

df.shape # 5491
df.ChEMBLID.nunique() # 4290

df.to_csv("data/Chembl_drug_targets_ChEMBLID_UniProtIDs_EnsemblIDs.csv", index = False)