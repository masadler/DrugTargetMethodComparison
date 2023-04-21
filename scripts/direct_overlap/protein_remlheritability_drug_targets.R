library(stringr)

#### Drug targets ####
n_databases = 5
drug_Ruiz_dgidb_df = read.csv(snakemake@input[["drug_Ruiz_dgidb"]], header = T)
drug_Ruiz_stitch_df = read.csv(snakemake@input[["drug_Ruiz_stitch"]], header = T)
drug_drugbank_dgidb_df = read.csv(snakemake@input[["drug_drugbank_dgidb"]], header = T)
drug_drugbank_stitch_df = read.csv(snakemake@input[["drug_drugbank_stitch"]], header = T)
drug_chembl_df = read.csv(snakemake@input[["drug_chembl"]], header = T)

# store drug target genes in a list
drug_targets = list()
drug_targets[[1]] = drug_Ruiz_dgidb_df$EnsemblId
drug_targets[[2]] = drug_Ruiz_stitch_df$EnsemblId
drug_targets[[3]] = drug_drugbank_dgidb_df$EnsemblId
drug_targets[[4]] = drug_drugbank_stitch_df$EnsemblId
drug_targets[[5]] = drug_chembl_df$EnsemblId
drug_targets_db = c("Ruiz_dgidb", "Ruiz_stitch", "drugbank_dgidb", "drugbank_stitch", "chembl")

#### MR - deCODE proteins ####
protein_df = read.table(snakemake@input[["mr"]], header = T)

gene_df = read.table(snakemake@input[["decode"]])
names(gene_df) = c("Chr", "ProbeID", "dummy", "bp", "ProbeName", "strand")

protein_df = merge(gene_df, protein_df[, c("ProbeID", "n_iv", "topSNP")], by = "ProbeID", all.x = TRUE)
protein_df[is.na(protein_df$n_iv), "n_iv"] = 0

# add heritability
her_df = read.csv(snakemake@input[["her"]], sep = '\t')
protein_df = merge(protein_df, her_df, by = "ProbeID")

# add EnsemblId
id_df <- str_split_fixed(protein_df$ProbeID, "_", 5)
protein_df$SeqId <- str_c(id_df[, 1], "_", id_df[, 2])
protein_info_df <- read.csv(snakemake@input[['protein_info']], header = T)
protein_info_df <- protein_info_df[, c("SeqId", "EnsemblId")]
protein_info_df <- protein_info_df[!duplicated(protein_info_df),]
protein_df <- merge(protein_df, protein_info_df, by = "SeqId")
protein_df = subset(protein_df, select=-c(ProbeID)) # remove ProbeID -> replace name with EnsemblId
protein_df = protein_df[,c("EnsemblId", "h2cis", "n_iv", "topSNP")]
names(protein_df) = c("ProbeID", "h2cis", "n_iv", "topSNP")

#### Result dataframes & Enrichment function ####
colnames = c("Universe", "Database", "n_targets", "n_genes", "h2cis_mean_target", "h2cis_mean_nottarget", 
"h2cis_stderr", "h2cis_tstat", "h2cis_pval", "poly_diff_loc", "poly_tstat", "poly_pval")
result_df = data.frame(matrix(NA, n_databases, length(colnames)))
names(result_df) = colnames

enrichment <- function(df, drug_genes, universe, database, result_df, i){
        
    # split genes into drug targets and not drug targets
    df$target = 0
    df[df$ProbeID %in% drug_genes, "target"] = 1

    n_targets = sum(df$target)
    result_df[i, "n_targets"] = n_targets
    result_df[i, "n_genes"] = nrow(df)
    print(n_targets)

    if (n_targets > 3){
        # T-test for enrichment - h2cis
        res = t.test(h2cis ~ target, data = df)
        result_df[i, "h2cis_mean_target"] = res$estimate[2]
        result_df[i, "h2cis_mean_nottarget"] = res$estimate[1]
        result_df[i, "h2cis_stderr"] = res$stderr
        result_df[i, "h2cis_tstat"] = res$statistic
        result_df[i, "h2cis_pval"] = res$p.value
        print("h2cis done")

        # Wilcoxon rank test for polygenicity

        if (n_targets > 3){
            res = wilcox.test(n_iv ~ target, data = df, conf.int = TRUE)
            result_df[i, "poly_diff_loc"] = res$estimate
            result_df[i, "poly_tstat"] = res$statistic
            result_df[i, "poly_pval"] = res$p.value
            print("n_iv done")
        }
    }

    # fill in context
    result_df[i, "Universe"] = universe
    result_df[i, "Database"] = database

    return(result_df)
}

#### Calculate enrichments ####

j = 1

for (i in 1:5){
    print(i)
    result_df = enrichment(protein_df, drug_targets[[i]], "decode", drug_targets_db[i], result_df, j)
    j = j+1
}

write.table(result_df, snakemake@output[[1]], sep = '\t', row.names = FALSE, quote = FALSE)