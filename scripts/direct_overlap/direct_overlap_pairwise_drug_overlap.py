import numpy as np
import pandas as pd

#### prepare result matrices ####

results = np.zeros((4,30))

mr_files = snakemake.input["mr"]
traits = [file.split("/")[2].split("_")[0] for file in mr_files]

#### get input files ####
drug_files = snakemake.input["drug_drugbank_dgidb"]
pascal_files = snakemake.input["pascal"]
exome_files = snakemake.input["exome"]
protein_files = snakemake.input["mr_protein"]
gene_info_df = pd.read_csv(snakemake.input["gene_info"], names = ["EnsemblId", "Chr", "GeneStart", "GeneEnd", "strand", "GeneName", "band"], sep = '\t')
protein_info_df = pd.read_csv(snakemake.input['protein_info'])
protein_info_df = protein_info_df[["SeqId", "EnsemblId"]].drop_duplicates()

#### calculate pairwise Jaccard indices ####

def Jaccard_calculation(df_1, df_2, p_1 = 0.01, p_2 = 0.01):

    # get prioritized genes
    pval_1 = np.quantile(df_1.p_value, p_1)
    pval_2 = np.quantile(df_2.p_value, p_2)

    df_1["prioritized"] = 0
    df_1["prioritized"][df_1.p_value <= pval_1] = 1

    df_2["prioritized"] = 0
    df_2["prioritized"][df_2.p_value <= pval_2] = 1

    # merge datasets to be in same universe
    df_merge = pd.merge(df_1, df_2, on = "EnsemblId", suffixes = ("_1", "_2"))
    df_merge = df_merge.drop_duplicates()

    targets_1 = df_merge["EnsemblId"][(df_merge.drug_target_1 == 1) & (df_merge.prioritized_1 == 1)]
    targets_2 = df_merge["EnsemblId"][(df_merge.drug_target_2 == 1) & (df_merge.prioritized_2 == 1)]

    intersection = len(list(set(targets_1).intersection(targets_2)))
    union = (len(targets_1) + len(targets_2)) - intersection
    if union == 0:
        jaccard_index = 0
    else:
        jaccard_index = float(intersection) / union

    return jaccard_index

for i, trait in enumerate(traits):

    #### Drug targets ####
    drug_df = pd.read_csv(drug_files[i])
    drug_genes = drug_df.EnsemblId

    #### Pascal ####
    pascal_df = pd.read_csv(pascal_files[i], sep = '\t', names = ["EnsemblId", "p_value"]) # analysis was only done on protein-coding genes

    # add drug target information
    pascal_df["drug_target"] = 0
    pascal_df["drug_target"][pascal_df.EnsemblId.isin(drug_genes)] = 1

    #### MR ####
    mr_df = pd.read_csv(mr_files[i], sep = '\t')
    mr_df = mr_df.rename(columns = {"ProbeID": "EnsemblId"})
    mr_df = mr_df[mr_df.EnsemblId.isin(gene_info_df.EnsemblId)]

    mr_df["drug_target"] = 0
    mr_df["drug_target"][mr_df.EnsemblId.isin(drug_genes)] = 1

    mr_eqtlgen_df = mr_df[["EnsemblId", "p_eqtlgen", "drug_target"]]
    mr_eqtlgen_df = mr_eqtlgen_df[~(mr_eqtlgen_df.p_eqtlgen).isna()]
    mr_eqtlgen_df = mr_eqtlgen_df.rename(columns = {"p_eqtlgen": "p_value"})

    mr_df = mr_df[["EnsemblId", "pmin", "drug_target"]]
    mr_df = mr_df.rename(columns = {"pmin": "p_value"})

    #### MR - protein
    protein_df = pd.read_csv(protein_files[i], sep = '\t')
    protein_df['SeqId'] = protein_df.ProbeID.apply(lambda x: x.split("_")[0] + '_' + x.split("_")[1])
    protein_df = protein_df.merge(protein_info_df, on = "SeqId")
    protein_df = protein_df[["EnsemblId", "p_ivw"]]
    protein_df = protein_df.rename(columns = {"p_ivw": "p_value"})
    protein_df["drug_target"] = 0
    protein_df["drug_target"][protein_df.EnsemblId.isin(drug_genes)] = 1

    # subset eQTLGen to decode universe
    universe_genes = list(set(mr_eqtlgen_df.EnsemblId) &  set(protein_df.EnsemblId))
    mr_eqtlgen_decode_df = mr_eqtlgen_df[mr_eqtlgen_df.EnsemblId.isin(universe_genes)]

    #### Exome ####
    exome_df = pd.read_csv(exome_files[i], sep = '\t')
    exome_df = exome_df[["EnsemblID", "p_value"]]
    exome_df = exome_df.rename(columns = {"EnsemblID": "EnsemblId"})
    exome_df["drug_target"] = 0
    exome_df["drug_target"][exome_df.EnsemblId.isin(drug_genes)] = 1

    #### Calculate pairwise Jaccard indices
    results[0,i] = Jaccard_calculation(pascal_df, mr_df)
    results[1,i] = Jaccard_calculation(pascal_df, exome_df)
    results[2,i] = Jaccard_calculation(mr_df, exome_df)
    results[3,i] = Jaccard_calculation(mr_eqtlgen_decode_df, protein_df, p_1 = 0.05, p_2 = 0.05)


result_df = pd.DataFrame(results, index = ["pascal_mr", "pascal_exome", "mr_exome", "transcript_protein"])
result_df.columns = traits

result_df.to_csv(snakemake.output[0], sep = '\t')