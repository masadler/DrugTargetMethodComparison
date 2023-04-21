import pandas as pd

df = pd.read_csv(snakemake.input['drugs'])
gene_df = pd.read_csv(snakemake.input['target_map'])

target_df = pd.read_csv(snakemake.input['drug_target'], sep = '\t')

target_df = target_df[["Parent Molecule ChEMBL ID", "Parent Molecule Name", "Parent Molecule Type",
                       "Mechanism of Action", "Target ChEMBL ID", "Target Name"]]

target_df = target_df.rename(columns = {"Parent Molecule ChEMBL ID": "Chembl_id", "Parent Molecule Name": "drug_name", 
                                        "Parent Molecule Type": "drug_type", "Mechanism of Action": "MoA", 
                                        "Target ChEMBL ID": "ChEMBLID_target", "Target Name": "Target_name"})

# select drugs of interest
target_df = target_df[target_df.Chembl_id.isin(df.Chembl_id)]

# merge UniProt and EnsemblId information
gene_df = gene_df.rename(columns = {"ChEMBLID": "ChEMBLID_target"})
target_df = target_df.merge(gene_df, on = "ChEMBLID_target")
target_df = target_df.rename(columns = {"EnsemblID": "EnsemblId"})

target_df = target_df[["Chembl_id", "drug_name", "EnsemblId", "drug_type", "MoA", "Target_name"]]

target_df.to_csv(snakemake.output[0], index = False)