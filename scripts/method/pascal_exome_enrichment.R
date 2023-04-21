

pascal_files = snakemake@input[["pascal"]]
pascal_UKBB_files = snakemake@input[["pascal_UKBB"]]
exome_files = snakemake@input[["exome"]]

num_files = length(pascal_files)

enrichment <- function(genes_1, genes_2, universe, results_df, TP_col, OR_col, pval_col, i){
        TP = length(intersect(genes_1, genes_2))
        FP = length(setdiff(genes_1, genes_2))
        FN = length(setdiff(genes_2, genes_1))
        TN = length(setdiff(universe, union(genes_1, genes_2)))
        cont.matrix = matrix(c(TP, FP, FN, TN), nrow = 2)
        fish.test = fisher.test(cont.matrix, alternative = "two.sided")
        pval = fish.test$p.val
        OR = fish.test$estimate
        results_df[i, TP_col] = TP
        results_df[i, OR_col] = OR
        results_df[i, pval_col] = pval

        return(results_df)
}

percentiles = c(0.1, 0.2, 0.5, 1, 2, 3, 5, 7.5, 10) # in percentage
num_perc = length(percentiles)
results_df = data.frame(matrix(NA, num_perc*num_files, 11))
names(results_df) = c("trait", "n_universe", "n_universe_UKBB", "n_exome", "perc", "TP", "OR", "pval", "UKBB_TP", "UKBB_OR", "UKBB_pval")

index = 1
for (file in 1:num_files){

    trait = strsplit(strsplit(pascal_files[file], "_")[[1]][5], ".tsv")[[1]][1]
    print(trait)
    results_df[index:(index+num_perc-1), 'trait'] = trait

    pascal_df = read.table(pascal_files[file])
    names(pascal_df) = c("gene", "pval")
    pascal_UKBB_df = read.table(pascal_UKBB_files[file])
    names(pascal_UKBB_df) = c("gene", "pval")
    print(head(pascal_df))
    exome_df = read.table(exome_files[file], header = T)
    print(head(exome_df))
    results_df[index:(index+num_perc-1), 'n_exome'] = nrow(exome_df)

    # Consortia
    universe = intersect(pascal_df$gene, exome_df$EnsemblID)
    print(paste0("Length of universe: ", length(universe)))
    pascal_df = pascal_df[pascal_df$gene %in% universe,]
    exome_df_cons = exome_df[exome_df$EnsemblID %in% universe,]
    results_df[index:(index+num_perc-1), 'n_universe'] = length(universe)

    # UKBB
    universe_UKBB = intersect(pascal_UKBB_df$gene, exome_df$EnsemblID)
    print(paste0("Length of universe: ", length(universe_UKBB)))
    pascal_UKBB_df = pascal_UKBB_df[pascal_UKBB_df$gene %in% universe_UKBB,]
    exome_df = exome_df[exome_df$EnsemblID %in% universe_UKBB,]
    results_df[index:(index+num_perc-1), 'n_universe_UKBB'] = length(universe_UKBB)

    for (j in 1:length(percentiles)){

        # Consortia data
        pascal_p = quantile(pascal_df$pval, percentiles[j]/100)
        print(paste0("Percentile: ", percentiles[j], "; pascal p-value:", pascal_p))
        pascal_genes = pascal_df[pascal_df$pval <= pascal_p, 'gene']
        exome_p = quantile(exome_df_cons$p_value, percentiles[j]/100)
        exome_genes = exome_df_cons[exome_df_cons$p_value <= exome_p, 'EnsemblID']
        results_df = enrichment(pascal_genes, exome_genes, universe, results_df, "TP", "OR", "pval", index + j -1)

        # UKBB data
        pascal_p = quantile(pascal_UKBB_df$pval, percentiles[j]/100)
        print(paste0("Percentile: ", percentiles[j], "; pascal p-value:", pascal_p))
        pascal_genes = pascal_UKBB_df[pascal_UKBB_df$pval <= pascal_p, 'gene']
        exome_p = quantile(exome_df$p_value, percentiles[j]/100)
        exome_genes = exome_df[exome_df$p_value <= exome_p, 'EnsemblID']
        results_df = enrichment(pascal_genes, exome_genes, universe_UKBB, results_df, "UKBB_TP", "UKBB_OR", "UKBB_pval", index + j -1)

        results_df[index + j -1, 'perc'] = percentiles[j]
    }
    index = index + num_perc
}

write.table(results_df, snakemake@output[[1]], row.names = F, quote = F, sep = '\t')
