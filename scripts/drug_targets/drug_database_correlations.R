gene_info_df = read.table(snakemake@input[["gene_info"]])
names(gene_info_df) = c("EnsemblId", "Chr", "GeneStart", "GeneEnd", "strand", "GeneName", "band")

#### drug targets ####
n_databases = 5
drug_targets = vector(mode = "list", length = n_databases)
drug_targets_db = c("Ruiz_dgidb", "Ruiz_stitch", "drugbank_dgidb", "drugbank_stitch", "chembl")

#print(snakemake@input[["drug_Ruiz_dgidb"]])

for (i in 1:length(snakemake@input[["drug_Ruiz_dgidb"]])){

    #### Drug targets ####
    drug_Ruiz_dgidb_df = read.csv(snakemake@input[["drug_Ruiz_dgidb"]][[i]], header = T)
    drug_Ruiz_stitch_df = read.csv(snakemake@input[["drug_Ruiz_stitch"]][[i]], header = T)
    drug_drugbank_dgidb_df = read.csv(snakemake@input[["drug_drugbank_dgidb"]][[i]], header = T)
    drug_drugbank_stitch_df = read.csv(snakemake@input[["drug_drugbank_stitch"]][[i]], header = T)
    drug_chembl_df = read.csv(snakemake@input[["drug_chembl"]][[i]], header = T)

    drug_targets[[1]] = c(drug_targets[[1]], as.vector(drug_Ruiz_dgidb_df$EnsemblId))
    drug_targets[[2]] = c(drug_targets[[2]], as.vector(drug_Ruiz_stitch_df$EnsemblId))
    drug_targets[[3]] = c(drug_targets[[3]], as.vector(drug_drugbank_dgidb_df$EnsemblId))
    drug_targets[[4]] = c(drug_targets[[4]], as.vector(drug_drugbank_stitch_df$EnsemblId))
    drug_targets[[5]] = c(drug_targets[[5]], as.vector(drug_chembl_df$EnsemblId))

}

d_df = gene_info_df[, c("EnsemblId", "GeneName")]

d_df$Ruiz_dgidb = 0
d_df[d_df$EnsemblId %in% drug_targets[[1]], "Ruiz_dgidb"] = 1

d_df$Ruiz_stitch = 0
d_df[d_df$EnsemblId %in% drug_targets[[2]], "Ruiz_stitch"] = 1

d_df$drugbank_dgidb = 0
d_df[d_df$EnsemblId %in% drug_targets[[3]], "drugbank_dgidb"] = 1

d_df$drugbank_stitch = 0
d_df[d_df$EnsemblId %in% drug_targets[[4]], "drugbank_stitch"] = 1

d_df$chembl = 0
d_df[d_df$EnsemblId %in% drug_targets[[5]], "chembl"] = 1

#### Correlation ####

d_df = d_df[, c("Ruiz_dgidb", "Ruiz_stitch", "drugbank_dgidb", "drugbank_stitch", "chembl")]

R = cor(d_df)

write.table(R, snakemake@output[[1]], quote = F)
