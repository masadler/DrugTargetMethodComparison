

pascal_files = snakemake@input[["pascal"]]
mr_files = snakemake@input[["mr"]]
gtex_files = snakemake@input[["mr_gtex"]]
gene_info_df = read.table(snakemake@input[["gene_info"]])
names(gene_info_df) = c("EnsemblId", "Chr", "GeneStart", "GeneEnd", "strand", "GeneName", "band")

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
colnames = c("trait", "n_universe_eqtlgen", "n_universe_eqtlgen_GTEX", "n_pascal", "n_universe_GTEX", "n_eqtlgen", "n_eqtlgen_GTEX", "n_GTEX", "perc", "eqtlgen_TP", "eqtlgen_OR", "eqtlgen_pval", "eqtlgen_GTEX_TP", "eqtlgen_GTEX_OR", "eqtlgen_GTEX_pval", "GTEX_TP", "GTEX_OR", "GTEX_pval")
results_df = data.frame(matrix(NA, num_perc*num_files, length(colnames)))
names(results_df) = colnames

index = 1
for (file in 1:num_files){

    trait = strsplit(strsplit(pascal_files[file], "_")[[1]][5], ".tsv")[[1]][1]
    print(trait)
    results_df[index:(index+num_perc-1), 'trait'] = trait

    pascal_df = read.table(pascal_files[file])
    names(pascal_df) = c("gene", "pval")
    results_df[index:(index+num_perc-1), 'n_pascal'] = nrow(pascal_df)
    print(head(pascal_df))
    mr_df = read.csv(mr_files[file], header = T, sep = '\t')
    mr_df = mr_df[mr_df$ProbeID %in% gene_info_df$EnsemblId,]

    gtex_df = read.csv(gtex_files[file], header = T, sep = '\t')
    gtex_df = gtex_df[gtex_df$ProbeID %in% gene_info_df$EnsemblId,]

    # distinguish between eQTLGen & eQTLGen + GTEx & GTEx
    eqtlgen_df = mr_df[!is.na(mr_df$p_eqtlgen),]
    results_df[index:(index+num_perc-1), 'n_eqtlgen'] = nrow(eqtlgen_df)
    results_df[index:(index+num_perc-1), 'n_eqtlgen_GTEX'] = nrow(mr_df)
    results_df[index:(index+num_perc-1), 'n_GTEX'] = nrow(gtex_df)

    universe_eqtlgen = intersect(pascal_df$gene, eqtlgen_df$ProbeID)
    results_df[index:(index+num_perc-1), 'n_universe_eqtlgen'] = length(universe_eqtlgen)
    pascal_eqtlgen_df = pascal_df[pascal_df$gene %in% universe_eqtlgen,]
    eqtlgen_df = eqtlgen_df[eqtlgen_df$ProbeID %in% universe_eqtlgen,]

    universe_GTEX = intersect(pascal_df$gene, gtex_df$ProbeID)
    results_df[index:(index+num_perc-1), 'n_universe_GTEX'] = length(universe_GTEX)
    pascal_gtex_df = pascal_df[pascal_df$gene %in% universe_GTEX,]
    gtex_df = gtex_df[gtex_df$ProbeID %in% universe_GTEX,]

    universe_eqtlgen_GTEX = intersect(pascal_df$gene, mr_df$ProbeID)
    results_df[index:(index+num_perc-1), 'n_universe_eqtlgen_GTEX'] = length(universe_eqtlgen_GTEX)
    pascal_df = pascal_df[pascal_df$gene %in% universe_eqtlgen_GTEX,]
    mr_df = mr_df[mr_df$ProbeID %in% universe_eqtlgen_GTEX,]
    
    for (j in 1:length(percentiles)){
        pascal_p = quantile(pascal_eqtlgen_df$pval, percentiles[j]/100)
        pascal_genes = pascal_eqtlgen_df[pascal_eqtlgen_df$pval <= pascal_p, 'gene']
        mr_p = quantile(eqtlgen_df$p_eqtlgen, percentiles[j]/100)
        mr_genes = eqtlgen_df[eqtlgen_df$p_eqtlgen <= mr_p, 'ProbeID']
        results_df = enrichment(pascal_genes, mr_genes, universe_eqtlgen, results_df, "eqtlgen_TP", "eqtlgen_OR", "eqtlgen_pval", index + j -1)

        pascal_p = quantile(pascal_gtex_df$pval, percentiles[j]/100)
        pascal_genes = pascal_gtex_df[pascal_gtex_df$pval <= pascal_p, 'gene']
        mr_p = quantile(gtex_df$pmin, percentiles[j]/100)
        mr_genes = gtex_df[gtex_df$pmin <= mr_p, 'ProbeID']
        results_df = enrichment(pascal_genes, mr_genes, universe_GTEX, results_df, "GTEX_TP", "GTEX_OR", "GTEX_pval", index + j -1)

        pascal_p = quantile(pascal_df$pval, percentiles[j]/100)
        pascal_genes = pascal_df[pascal_df$pval <= pascal_p, 'gene']
        mr_p = quantile(mr_df$pmin, percentiles[j]/100)
        mr_genes = mr_df[mr_df$pmin <= mr_p, 'ProbeID']
        results_df = enrichment(pascal_genes, mr_genes, universe_eqtlgen_GTEX, results_df, "eqtlgen_GTEX_TP", "eqtlgen_GTEX_OR", "eqtlgen_GTEX_pval", index + j -1)
        results_df[index + j -1, 'perc'] = percentiles[j]
    }
    index = index + num_perc
}

write.table(results_df, snakemake@output[[1]], row.names = F, quote = F, sep = '\t')
