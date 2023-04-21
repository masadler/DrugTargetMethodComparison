import pandas as pd

dd_df = pd.read_csv(snakemake.input["drug_indication_data"])
gwas_indication_df = pd.read_csv(snakemake.input["gwas_indication"])

trait_df = gwas_indication_df[gwas_indication_df.Acronym == snakemake.wildcards["trait"]]

trait_df = (trait_df.set_index(trait_df.columns.drop('Ruiz_ID',1).tolist())
                .Ruiz_ID.str.split('; ', expand=True)
                .stack()
                .reset_index()
                .rename(columns={0:'Ruiz_ID'})
                .loc[:, trait_df.columns]
                )

# make sure all indications where found!
missing_df = trait_df[~trait_df.Ruiz_ID.isin(dd_df.indication)]
if missing_df.shape[0] > 0:
    print("These traits were missing:")
    print(missing_df[["Ruiz_ID"]])
    exit()

dd_df_trait = dd_df[dd_df.indication.isin(trait_df.Ruiz_ID)]

dd_df_trait = dd_df_trait[["drug", "drug_name"]]
dd_df_trait = dd_df_trait.drop_duplicates()
dd_df_trait = dd_df_trait.rename(columns = {"drug": "drugbank_id", "drug_name": "db_name"})

dd_df_trait.to_csv(snakemake.output[0], index = False)
