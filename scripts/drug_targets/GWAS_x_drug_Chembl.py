import pandas as pd

df = pd.read_csv(snakemake.input['drug_indication_data'], sep = '\t')
trait_df = pd.read_csv(snakemake.input['gwas_indication'])

trait_df = trait_df[trait_df.Acronym == snakemake.wildcards["trait"]]

trait_df = (trait_df.set_index(trait_df.columns.drop('ChEMBL_ID',1).tolist())
                .ChEMBL_ID.str.split('; ', expand=True)
                .stack()
                .reset_index()
                .rename(columns={0:'ChEMBL_ID'})
                .loc[:, trait_df.columns]
                )

df = df[["Parent Molecule ChEMBL ID", "Parent Molecule Name", "Parent Molecule Type", "EFO IDs", "Max Phase for Indication"]]
df = df.rename(columns = {"Parent Molecule ChEMBL ID": "Chembl_id", "Parent Molecule Name": "drug_name", 
                          "Parent Molecule Type": "drug_type", "EFO IDs": "EFO_ID"})

# explode EFO_ID column
df = (df.set_index(df.columns.drop('EFO_ID',1).tolist())
                .EFO_ID.str.split('|', expand=True)
                .stack()
                .reset_index()
                .rename(columns={0:'EFO_ID'})
                .loc[:, df.columns]
                )

# make sure all indications where found!
missing_df = trait_df[~trait_df.ChEMBL_ID.isin(df.EFO_ID)]
if missing_df.shape[0] > 0:
    print("These traits were missing:")
    print(missing_df[["ChEMBL_ID"]])
    exit()

# only drugs that passed clinical phase 4
df = df[df["Max Phase for Indication"] == 4]

df = df[df.EFO_ID.isin(trait_df.ChEMBL_ID)]

df = df[["Chembl_id", "drug_name", "drug_type"]]
df = df.drop_duplicates()

df.to_csv(snakemake.output[0], index = False)