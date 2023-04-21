import numpy as np
import pandas as pd

gene_info_df = pd.read_csv(snakemake.input["gene_info"], sep = '\t', names = ["EnsemblId", "Chr", "GeneStart", "GeneEnd", "strand", "GeneName", "band"])
tissue_df = pd.read_csv(snakemake.input["tissue"])

mr_files = snakemake.input["mr"]

# get tissue for each gene
def get_col_name(row):    
    b = (df.ix[row.name] == row['pmin'])
    return b.idxmax()

# create result dataframe 
result_df = pd.DataFrame()

for file in mr_files:

    trait = file.split("/")[2].split("_")[0]
    print(trait)

    df = pd.read_csv(file, sep = '\t')

    df_t = df[["ProbeID", "pmin"]]
    df = df.drop('pmin', axis = 1)
    df_t['p_tissue'] = df_t.apply(get_col_name, axis=1)
    df_t['Tissue'] = df_t.apply(lambda x: x['p_tissue'].split("p_")[1], axis = 1)

    # subset to protein-coding genes
    df_t = df_t[df_t.ProbeID.isin(gene_info_df.EnsemblId)]

    # select 1%
    mr_p = np.quantile(df_t.pmin, 0.01)
    df_t = df_t[df_t.pmin <= mr_p]
    df_t = df_t.merge(tissue_df, on = "Tissue")

    count_df = df_t["Tissue_category"].value_counts().rename_axis('Tissue').reset_index(name='counts')
    count_df['trait'] = trait

    result_df = pd.concat([result_df, count_df])

result_df.to_csv(snakemake.output[0], sep = '\t', index = False)





