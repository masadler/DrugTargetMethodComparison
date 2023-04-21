
seed = as.integer(snakemake@params[["seed"]])

gene_info_df = read.table(snakemake@input[["gene_info"]])
names(gene_info_df) = c("EnsemblId", "Chr", "GeneStart", "GeneEnd", "strand", "GeneName", "band")
gene_order = unique(gene_info_df$EnsemblId)

set.seed(seed)
random_scores <- rnorm(length(gene_order))
random_df <- data.frame(EnsemblId = gene_order,
                         z = random_scores)

write.table(random_df, snakemake@output[[1]], row.names = F, quote = F, sep = '\t')