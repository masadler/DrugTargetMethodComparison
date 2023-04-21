import pandas as pd

ndf = pd.read_csv(snakemake.input["network"], delim_whitespace=True)
print(ndf.head())
ndf = ndf[ndf.combined_score >= 150]
ndf["protein1"] = ndf["protein1"].apply(lambda x: x[5:])
ndf["protein2"] = ndf["protein2"].apply(lambda x: x[5:])

# merge Ensembl Identifier
gene_df = pd.read_csv(snakemake.input["prot_gene_ensembl_voc"], sep = '\t')
gene_df = gene_df[["ensembl_gene_id", "ensembl_peptide_id"]]
ndf = ndf.merge(gene_df, left_on = "protein1", right_on = "ensembl_peptide_id")
ndf = ndf.rename(columns = {"ensembl_gene_id": "gene1"})
ndf = ndf.merge(gene_df, left_on = "protein2", right_on = "ensembl_peptide_id")
ndf = ndf.rename(columns = {"ensembl_gene_id": "gene2"})

ndf = ndf[["gene1", "gene2", "combined_score"]]
ndf = ndf.drop_duplicates()
ndf["weight"] = ndf.combined_score
ndf = ndf[["gene1", "gene2", "weight"]]

ndf.to_csv(snakemake.output[0], sep = '\t', index = False)