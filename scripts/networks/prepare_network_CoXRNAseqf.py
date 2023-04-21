import numpy as np
import pandas as pd

df = pd.read_csv(snakemake.input["network"], sep = '\t', index_col=0)
df[df < 4] = 0
df = df.replace(0, np.nan)

edge_df = df.stack().reset_index()
edge_df = edge_df.rename(columns = {"-": "gene1", "level_1": "gene2", 0: "weight"})
del df

# reduce network to protein-coding genes
gene_info_df = pd.read_csv(snakemake.input["gene_info"], sep = '\t', names = ["EnsemblId", "Chr", "GeneStart", "GeneEnd", "strand", "GeneName", "band"])
gene_info_df = gene_info_df[["EnsemblId", "GeneName"]]

edge_df = edge_df[(edge_df.gene1.isin(gene_info_df.EnsemblId)) & (edge_df.gene2.isin(gene_info_df.EnsemblId))]
edge_df = edge_df.dropna()
edge_df = edge_df[["gene1", "gene2", "weight"]]
edge_df = edge_df.drop_duplicates()

edge_df.to_csv(snakemake.output[0], sep = '\t', index = False)
