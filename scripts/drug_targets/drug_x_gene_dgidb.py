import pandas as pd

"""
DGIDB:
Data downloaded from DGIDB database (Integration of the Drug–Gene Interaction Database (DGIdb 4.0) 
with open crowdsource efforts. Freshour S, Kiwala S, Cotto KC, Coffman AC, McMichael JF, Song J, Griffith M, 
Griffith OL, Wagner AH. Nucleic Acids Research. 2020 Nov 25; doi: https://doi.org/10.1093/nar/gkaa1084. PMID: 33237278): 
- https://www.dgidb.org/downloads

Target selection:
- All DGIDB relationships
- gene_name has to be present

Gene ID:
- genes are matched to EnsemblID using the provided gene vocabulary file

Drug ID:
- if drugbank ID (DBXXX) is not provided, add mappings from PubChem (ChemblID -> DrugBank identifiers)
"""

drug_df = pd.read_csv(snakemake.input["drugs"])

# Map DGIDB identifiers to drugbank IDs
drug_voc_df = pd.read_csv(snakemake.input["drug_voc"], sep = '\t')
drug_voc_df = drug_voc_df.dropna(subset = ["drug_name"])
drug_voc_df['Chembl_id'] = drug_voc_df["concept_id"].apply(lambda x: x[7:])
#drug_voc_df = drug_voc_df[drug_voc_df.drug_claim_source == "DrugBank"]
#drug_dgidb_df = drug_df.merge(drug_voc_df, left_on = "drugbank_id", right_on = "drug_claim_name")
drug_chembl_db_df = pd.read_csv(snakemake.input["drug_voc_chembl_db"], sep = '\t', names = ["Chembl_id", "drugbank_id"])
drug_chembl_db_df = drug_chembl_db_df.dropna()
drug_voc_df = drug_voc_df.merge(drug_chembl_db_df, on = "Chembl_id", how = "left")
drug_voc_df["drugbank_id"][drug_voc_df.drug_claim_source == "DrugBank"] = drug_voc_df["drug_claim_name"][drug_voc_df.drug_claim_source == "DrugBank"]
drug_voc_df = drug_voc_df.dropna(subset = ["drugbank_id"])

# Select required drugs
drug_dgidb_df = drug_df.merge(drug_voc_df, on = "drugbank_id")

# Add EnsemblId information
gene_df = pd.read_csv(snakemake.input["gene_voc"], sep = '\t')
gene_df = gene_df[gene_df.gene_claim_source == "Ensembl"]
gene_df = gene_df.dropna(subset = ['gene_name'])

# get drug targets from DGIDB
df = pd.read_csv(snakemake.input["dgidb_data"], sep = '\t')
df = df.dropna(subset = ['gene_name', 'gene_claim_name', 'drug_concept_id'])

df = df[["gene_name", "interaction_claim_source", "drug_claim_name", "drug_name", "drug_concept_id", "interaction_group_score"]]
df = df.merge(drug_dgidb_df, left_on = 'drug_concept_id', right_on = 'concept_id')
df = df.merge(gene_df, on = 'gene_name')

df = df[['drugbank_id', 'db_name', 'gene_name', 'gene_claim_name']]
df = df.rename(columns = {'gene_claim_name': 'EnsemblId'})

df = df.drop_duplicates()

df.to_csv(snakemake.output[0], index = False)