import pandas as pd

path = snakemake.params["data_path"]

trait_df = pd.read_csv(snakemake.input["trait_info"])
traits = list(trait_df.Acronym_Exome)

for trait in traits:
    print(trait)
    print(path)
    print(str(trait_df["file_name"][trait_df.Acronym_Exome == trait].item()))
    df = pd.read_csv(path + str(trait_df["file_name"][trait_df.Acronym_Exome == trait].item()), sep = '\t')
    df = df[df.effect_allele == "M3.1"]
    df["GeneName"] = df['Name'].apply(lambda x: x.split('(')[0])
    df["EnsemblID"] = df['Name'].apply(lambda x: x.split('(')[1].split(')')[0])
    if trait in ["BMD", "DBP", "SBP", "LDL", "TC"]:
        df = df[["GeneName", "EnsemblID", "chromosome", "base_pair_location", "beta", "p_value", "effect_allele_frequency", "standard_error"]]
    else:
        df = df[["GeneName", "EnsemblID", "chromosome", "base_pair_location", "odds_ratio", "p_value", "effect_allele_frequency", "standard_error"]]
    df.to_csv("output/exome_genes/" + trait + "_gene_burden_M3_1_associations.tsv", sep = '\t', index = False)