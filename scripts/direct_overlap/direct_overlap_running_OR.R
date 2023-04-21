library(stringr)
library(dplyr) # df = rename(df, new_col01 = old_col01, new_col02 = old_col02, ...)

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

#### Pascal ####
pascal_df = read.table(snakemake@input[["pascal"]], header = F) # analysis was only done on protein-coding genes
names(pascal_df) = c("EnsemblId", "p_value")
pascal_UKBB_df = read.table(snakemake@input[["pascal_ukbb"]], header = F)
names(pascal_UKBB_df) = c("EnsemblId", "p_value")

#### MR ####
mr_df = read.csv(snakemake@input[["mr"]], header = T, sep = '\t')
gene_info_df = read.table(snakemake@input[["gene_info"]])
names(gene_info_df) = c("EnsemblId", "Chr", "GeneStart", "GeneEnd", "strand", "GeneName", "band")

mr_UKBB_df = read.csv(snakemake@input[["mr_UKBB"]], header = T, sep = '\t')
mr_UKBB_df = mr_UKBB_df[mr_UKBB_df$ProbeID %in% gene_info_df$EnsemblId,]
mr_UKBB_df = mr_UKBB_df[, c("ProbeID", "pmin")]
names(mr_UKBB_df) = c("EnsemblId", "p_value")

mr_df = mr_df[mr_df$ProbeID %in% gene_info_df$EnsemblId,]
mr_eqtlgen_df = mr_df[, c("ProbeID", "p_eqtlgen")]
mr_eqtlgen_df = mr_eqtlgen_df[!is.na(mr_eqtlgen_df$p_eqtlgen),]
names(mr_eqtlgen_df) = c("EnsemblId", "p_value")

mr_df = mr_df[, c("ProbeID", "pmin")]
names(mr_df) = c("EnsemblId", "p_value")

mr_gtex_df = read.csv(snakemake@input[["mr_gtex"]], header = T, sep = '\t')
mr_gtex_df = mr_gtex_df[, c("ProbeID", "pmin")]
mr_gtex_df = mr_gtex_df[mr_gtex_df$ProbeID %in% gene_info_df$EnsemblId,]
names(mr_gtex_df) = c("EnsemblId", "p_value")

#### MR - protein
protein_df = read.table(snakemake@input[["mr_protein"]], header = T)
id_df <- str_split_fixed(protein_df$ProbeID, "_", 5)
protein_df$SeqId <- str_c(id_df[, 1], "_", id_df[, 2])
protein_info_df <- read.csv(snakemake@input[['protein_info']], header = T)
protein_info_df <- protein_info_df[, c("SeqId", "EnsemblId")]
protein_info_df <- protein_info_df[!duplicated(protein_info_df),]
protein_df <- merge(protein_df, protein_info_df, by = "SeqId")
protein_df = protein_df[,c("EnsemblId", "p_ivw")]
names(protein_df) = c("EnsemblId", "p_value")

#### Exome ####
exome_df = read.table(snakemake@input[["exome"]], header = T)
exome_df = exome_df[,c("EnsemblID", "p_value")]
names(exome_df) = c("EnsemblId", "p_value")

#### Result dataframe & Enrichment function ####
percentiles = c(0.1, 0.2, 0.5, 1, 2, 3, 5, 7.5, 10)
colnames = c("Percentile", "Method", "Universe", "Database", "TP", "OR", "Pval", "diagonal_padding", "ci_low", "ci_up")
results_df = data.frame(matrix(NA, n_databases*12*length(percentiles), length(colnames)))
names(results_df) = colnames

enrichment <- function(method_genes, drug_genes, universe_genes, perc, method, universe, database, results_df, i){
        
        # intersect with universe
        method_genes = intersect(method_genes, universe_genes)
        drug_genes = intersect(drug_genes, universe_genes)

        # calculate enrichment  
        diagonal_padding = FALSE

        TP = length(intersect(method_genes, drug_genes))
        FP = length(setdiff(method_genes, drug_genes))
        FN = length(setdiff(drug_genes, method_genes))
        TN = length(setdiff(universe_genes, union(drug_genes, method_genes)))

        if (FP == 0){
            diagonal_padding = TRUE
            FP = 1
        }

        if (FN == 0){
            diagonal_padding = TRUE
            FN = 1
        }

        cont.matrix = matrix(c(TP, FP, FN, TN), nrow = 2)
        fish.test = fisher.test(cont.matrix, conf.int = TRUE, conf.level = 0.95)
        pval = fish.test$p.val
        OR = fish.test$estimate
        results_df[i, "TP"] = TP
        results_df[i, "OR"] = OR
        results_df[i, "Pval"] = pval
        results_df[i, "ci_low"] = fish.test$conf.int[1]
        results_df[i, "ci_up"] = fish.test$conf.int[2]

        # fill in context
        results_df[i, "Percentile"] = perc
        results_df[i, "Method"] = method
        results_df[i, "Universe"] = universe
        results_df[i, "Database"] = database
        results_df[i, "diagonal_padding"] = diagonal_padding

        return(results_df)
}

#### Calculate enrichments ####

j = 1

## 1) Pascal
universe_genes = unique(pascal_df$EnsemblId)

for (i in 1:5){
    for (perc in percentiles){
        pascal_genes = unique(pascal_df[pascal_df$p_value <= quantile(pascal_df$p_value, perc/100), "EnsemblId"])
        results_df = enrichment(pascal_genes, drug_targets[[i]], universe_genes, perc, "Pascal", "Pascal", drug_targets_db[i], results_df, j)
        j = j+1
    }
}

## 2) Pascal UKBB
universe_genes = unique(pascal_UKBB_df$EnsemblId)

for (i in 1:5){
    for (perc in percentiles){
        pascal_genes = unique(pascal_UKBB_df[pascal_UKBB_df$p_value <= quantile(pascal_UKBB_df$p_value, perc/100), "EnsemblId"])
        results_df = enrichment(pascal_genes, drug_targets[[i]], universe_genes, perc, "Pascal_UKBB", "Pascal_UKBB", drug_targets_db[i], results_df, j)
        j = j+1
    }
}

## 3) Pascal - MR universe (GTEx & eQTLGen)
universe_genes = intersect(pascal_df$EnsemblId, mr_df$EnsemblId)
pascal_mr_df = pascal_df[pascal_df$EnsemblId %in% universe_genes,]

for (i in 1:5){
    for (perc in percentiles){
        pascal_genes = unique(pascal_mr_df[pascal_mr_df$p_value <= quantile(pascal_mr_df$p_value, perc/100), "EnsemblId"])
        results_df = enrichment(pascal_genes, drug_targets[[i]], universe_genes, perc, "Pascal", "Pascal_eqtlgen_gtex", drug_targets_db[i], results_df, j)
        j = j+1
    }
}

## 4) Pascal - MR universe (deCODE proteins)
universe_genes = intersect(pascal_df$EnsemblId, protein_df$EnsemblId)
pascal_mr_df = pascal_df[pascal_df$EnsemblId %in% universe_genes,]

for (i in 1:5){
    for (perc in percentiles){
        pascal_genes = unique(pascal_mr_df[pascal_mr_df$p_value <= quantile(pascal_mr_df$p_value, perc/100), "EnsemblId"])
        results_df = enrichment(pascal_genes, drug_targets[[i]], universe_genes, perc, "Pascal", "Pascal_decode", drug_targets_db[i], results_df, j)
        j = j+1
    }
}

## 5) MR - eQTLGen & GTEx
universe_genes = unique(mr_df$EnsemblId)

for (i in 1:5){
    for (perc in percentiles){
        mr_genes = unique(mr_df[mr_df$p_value <= quantile(mr_df$p_value, perc/100), "EnsemblId"])
        results_df = enrichment(mr_genes, drug_targets[[i]], universe_genes, perc, "MR_eqtlgen_gtex", "eqtlgen_gtex", drug_targets_db[i], results_df, j)
        j = j+1
    }
}

## 6) MR - eQTLGen & GTEx - UKBB
universe_genes = unique(mr_UKBB_df$EnsemblId)

for (i in 1:5){
    for (perc in percentiles){
        mr_genes = unique(mr_UKBB_df[mr_UKBB_df$p_value <= quantile(mr_UKBB_df$p_value, perc/100), "EnsemblId"])
        results_df = enrichment(mr_genes, drug_targets[[i]], universe_genes, perc, "MR_eqtlgen_gtex_UKBB", "eqtlgen_gtex_UKBB", drug_targets_db[i], results_df, j)
        j = j+1
    }
}

## 7) MR - eQTLGen
universe_genes = unique(mr_eqtlgen_df$EnsemblId)

for (i in 1:5){
    for (perc in percentiles){
        mr_genes = unique(mr_eqtlgen_df[mr_eqtlgen_df$p_value <= quantile(mr_eqtlgen_df$p_value, perc/100), "EnsemblId"])
        results_df = enrichment(mr_genes, drug_targets[[i]], universe_genes, perc, "MR_eqtlgen", "eqtlgen", drug_targets_db[i], results_df, j)
        j = j+1
    }
}

## 8) MR - GTEx
universe_genes = unique(mr_gtex_df$EnsemblId)

for (i in 1:5){
    for (perc in percentiles){
        mr_genes = unique(mr_gtex_df[mr_gtex_df$p_value <= quantile(mr_gtex_df$p_value, perc/100), "EnsemblId"])
        results_df = enrichment(mr_genes, drug_targets[[i]], universe_genes, perc, "MR_gtex", "gtex", drug_targets_db[i], results_df, j)
        j = j+1
    }
}

## 9) MR - eQTLGen (decode universe)
universe_genes = intersect(mr_eqtlgen_df$EnsemblId, protein_df$EnsemblId)
mr_eqtlgen_decode_df = mr_eqtlgen_df[mr_eqtlgen_df$EnsemblId %in% universe_genes,]

for (i in 1:5){
    for (perc in percentiles){
        mr_genes = unique(mr_eqtlgen_decode_df[mr_eqtlgen_decode_df$p_value <= quantile(mr_eqtlgen_decode_df$p_value, perc/100), "EnsemblId"])
        results_df = enrichment(mr_genes, drug_targets[[i]], universe_genes, perc, "MR_eqtlgen", "eqtlgen_decode", drug_targets_db[i], results_df, j)
        j = j+1
    }
}

## 10) MR - decode proteins
universe_genes = unique(protein_df$EnsemblId)

for (i in 1:5){
    for (perc in percentiles){
        mr_genes = unique(protein_df[protein_df$p_value <= quantile(protein_df$p_value, perc/100), "EnsemblId"])
        results_df = enrichment(mr_genes, drug_targets[[i]], universe_genes, perc, "MR_decode", "decode", drug_targets_db[i], results_df, j)
        j = j+1
    }
}

## 11) MR - decode proteins (eQTLGen universe)
universe_genes = intersect(mr_eqtlgen_df$EnsemblId, protein_df$EnsemblId)
protein_eqtlgen_df = protein_df[protein_df$EnsemblId %in% universe_genes,]

for (i in 1:5){
    for (perc in percentiles){
        mr_genes = unique(protein_eqtlgen_df[protein_eqtlgen_df$p_value <= quantile(protein_eqtlgen_df$p_value, perc/100), "EnsemblId"])
        results_df = enrichment(mr_genes, drug_targets[[i]], universe_genes, perc, "MR_decode", "eqtlgen_decode", drug_targets_db[i], results_df, j)
        j = j+1
    }
}

## 12) Exome
universe_genes = unique(exome_df$EnsemblId)

for (i in 1:5){
    for (perc in percentiles){
        exome_genes = unique(exome_df[exome_df$p_value <= quantile(exome_df$p_value, perc/100), "EnsemblId"])
        results_df = enrichment(exome_genes, drug_targets[[i]], universe_genes, perc, "Exome", "Exome", drug_targets_db[i], results_df, j)
        j = j+1
    }
}

write.table(results_df, snakemake@output[[1]], row.names = F, quote = F, sep = '\t')