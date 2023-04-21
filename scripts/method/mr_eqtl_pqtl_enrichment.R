library(stringr)

protein_info_df <- read.csv(snakemake@input[['protein_info']], header = T)
protein_info_df <- protein_info_df[, c("SeqId", "EnsemblId")]
protein_info_df <- protein_info_df[!duplicated(protein_info_df),]

transcript_files <- snakemake@input[['mr_e']]
protein_files <- snakemake@input[['mr_p']]

num_files = length(transcript_files)

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
results_df = data.frame(matrix(NA, num_perc*num_files, 7))
names(results_df) = c("trait", "n_universe", "n_protein", "perc", "TP", "OR", "pval")

index = 1
for (file in 1:num_files){

    trait <- strsplit(transcript_files[file], "_")[[1]][4]
    results_df[index:(index+num_perc-1), 'trait'] = trait

    t_df <- read.table(transcript_files[file], header = T)

    # add Ensembl information to protein dataframe
    p_df <- read.table(protein_files[file], header = T)
    id_df <- str_split_fixed(p_df$ProbeID, "_", 5)
    p_df$SeqId <- str_c(id_df[, 1], "_", id_df[, 2])
    p_df <- merge(p_df, protein_info_df, by = "SeqId")
    results_df[index:(index+num_perc-1), 'n_protein'] = length(unique(p_df$EnsemblId))

    df <- merge(t_df, p_df, by.x = "ProbeID", by.y = "EnsemblId", 
                suffixes = c("_t", "_p"))

    results_df[index:(index+num_perc-1), 'n_universe'] = nrow(df)
    universe = df$ProbeID

    for (j in 1:length(percentiles)){

        # Consortia data
        transcript_p = quantile(df$p_ivw_t, percentiles[j]/100)
        transcript_genes = df[df$p_ivw_t <= transcript_p, 'ProbeID']
        protein_p = quantile(df$p_ivw_p, percentiles[j]/100)
        protein_genes = df[df$p_ivw_p <= protein_p, 'ProbeID']
        results_df = enrichment(transcript_genes, protein_genes, universe, results_df, "TP", "OR", "pval", index + j -1)

        results_df[index + j -1, 'perc'] = percentiles[j]
    }
    index = index + num_perc
}

write.table(results_df, snakemake@output[[1]], row.names = F, quote = F, sep = '\t')
