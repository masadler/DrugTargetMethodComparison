library(stringr)
library(dplyr) # df = rename(df, new_col01 = old_col01, new_col02 = old_col02, ...)
library(pROC)

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

# store all drug/target information in a list
drug_Ruiz_dgidb_df = drug_Ruiz_dgidb_df[, c("drugbank_id", "db_name", "EnsemblId", "gene_name")]
names(drug_Ruiz_dgidb_df) = c("drug_id", "drug_name", "EnsemblId", "Gene")
drug_Ruiz_stitch_df = drug_Ruiz_stitch_df[, c("drugbank_id", "db_name", "EnsemblId", "gene_name")]
names(drug_Ruiz_stitch_df) = c("drug_id", "drug_name", "EnsemblId", "Gene")
drug_drugbank_dgidb_df = drug_drugbank_dgidb_df[, c("drugbank_id", "db_name", "EnsemblId", "gene_name")]
names(drug_drugbank_dgidb_df) = c("drug_id", "drug_name", "EnsemblId", "Gene")
drug_drugbank_stitch_df = drug_drugbank_stitch_df[, c("drugbank_id", "db_name", "EnsemblId", "gene_name")]
names(drug_drugbank_stitch_df) = c("drug_id", "drug_name", "EnsemblId", "Gene")

gene_voc_df = read.csv(snakemake@input[["prot_gene_ensembl_voc"]], sep = '\t')
gene_voc_df = gene_voc_df[, c("ensembl_gene_id", "hgnc_symbol")]
names(gene_voc_df) = c("EnsemblId", "Gene")
gene_voc_df = gene_voc_df[!duplicated(gene_voc_df),]
drug_chembl_df = merge(drug_chembl_df, gene_voc_df)
drug_chembl_df = drug_chembl_df[, c("Chembl_id", "drug_name", "EnsemblId", "Gene")]
names(drug_chembl_df) = c("drug_id", "drug_name", "EnsemblId", "Gene")
drug_chembl_df = drug_chembl_df[!duplicated(drug_chembl_df),]

drug_info = list()
drug_info[[1]] = drug_Ruiz_dgidb_df
drug_info[[2]] = drug_Ruiz_stitch_df
drug_info[[3]] = drug_drugbank_dgidb_df
drug_info[[4]] = drug_drugbank_stitch_df
drug_info[[5]] = drug_chembl_df

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

#### MR - Protein ####

protein_df = read.table(snakemake@input[["mr_protein"]], header = T)
id_df <- str_split_fixed(protein_df$ProbeID, "_", 5)
protein_df$SeqId <- str_c(id_df[, 1], "_", id_df[, 2])
protein_info_df <- read.csv(snakemake@input[['protein_info']], header = T)
protein_info_df <- protein_info_df[, c("SeqId", "EnsemblId")]
protein_info_df <- protein_info_df[!duplicated(protein_info_df),]
protein_df <- merge(protein_df, protein_info_df, by = "SeqId")
protein_df = protein_df[, c("EnsemblId", "p_ivw")]
names(protein_df) = c("EnsemblId", "p_value")

# UKBB
protein_ukbb_df = read.table(snakemake@input[["mr_protein_UKBB"]], header = T)
id_df <- str_split_fixed(protein_ukbb_df$ProbeID, "_", 5)
protein_ukbb_df$SeqId <- str_c(id_df[, 1], "_", id_df[, 2])
protein_ukbb_df <- merge(protein_ukbb_df, protein_info_df, by = "SeqId")
protein_ukbb_df = protein_ukbb_df[,c("EnsemblId", "p_ivw")]
names(protein_ukbb_df) = c("EnsemblId", "p_value")

#### Exome ####
exome_df = read.table(snakemake@input[["exome"]], header = T)
exome_df = exome_df[, c("EnsemblID", "p_value")]
names(exome_df) = c("EnsemblId", "p_value")

#### Random ####
random_df = read.table(snakemake@input[["random"]], header = T)
random_df$p_value = 2*pnorm(abs(random_df$z), lower.tail = F)
random_df = random_df[, c("EnsemblId", "p_value")]

#### Result dataframes & Enrichment function & drug info function ####
colnames = c("Percentile", "Method", "Universe", "Database", "TP", "OR", "Pval", "diagonal_padding", "ci_low", "ci_up", "AUC", "ci_low_AUC", "ci_up_AUC")
results_df = data.frame(matrix(NA, n_databases*15, length(colnames)))
names(results_df) = colnames

colnames = c("Percentile", "Method", "Universe", "Database", "TP", "EnsemblId", "Gene", "drug_id", "drug_name", "rank", "gene_perc")
target_df = data.frame(matrix(NA, 1, length(colnames)))
names(target_df) = colnames

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

AUC <- function(method_df, drug_genes, results_df, i){

    method_df$rank = rank(method_df$p_value, ties.method = "min")
    method_df$target = 0
    method_df[method_df$EnsemblId %in% drug_genes, "target"] = 1
    if (sum(method_df$target) == 0){
        results_df[i, "AUC"] = 0.5
        results_df[i, "ci_low_AUC"] = 0.975
        results_df[i, "ci_up_AUC"] = 0.025
    } else {
        roc1 <- roc(method_df$target,
                method_df$rank, direction = ">",
                ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE)
        results_df[i, "AUC"] = auc(roc1)[1]
        results_df[i, "ci_low_AUC"] = ci(roc1)[1]
        results_df[i, "ci_up_AUC"] = ci(roc1)[3]
    }

    return(results_df)
}

drug_target_info <- function(method_genes, method_df, drug_genes, drug_info, perc, method, universe, database, target_df){
    
    targets = intersect(method_genes, drug_genes)
    num_targets = length(targets)
    
    if (num_targets == 0) return (target_df)

    target_info = drug_info[drug_info$EnsemblId %in% targets, ]

    # get gene rank
    method_df$rank = rank(method_df$p_value, ties.method = "min")
    method_df$gene_perc = method_df$rank/nrow(method_df)*100
    target_info = merge(target_info, method_df[, c("EnsemblId", "rank", "gene_perc")])

    target_info$TP = num_targets
    target_info$Percentile = perc
    target_info$Method = method
    target_info$Universe = universe
    target_info$Database = database

    target_info = target_info[, c("Percentile", "Method", "Universe", "Database", "TP", "EnsemblId", "Gene", "drug_id", "drug_name", "rank", "gene_perc")]

    target_df = rbind(target_df, target_info)
    
    return (target_df)
}

#### Calculate enrichments ####

j = 1

## 1) Pascal
perc = 0.01
universe_genes = unique(pascal_df$EnsemblId)
pascal_genes = unique(pascal_df[pascal_df$p_value <= quantile(pascal_df$p_value, perc), "EnsemblId"])

for (i in 1:5){
    results_df = enrichment(pascal_genes, drug_targets[[i]], universe_genes, perc, "Pascal", "Pascal", drug_targets_db[i], results_df, j)
    results_df = AUC(pascal_df, drug_targets[[i]], results_df, j)
    target_df = drug_target_info(pascal_genes, pascal_df, drug_targets[[i]], drug_info[[i]], perc, "Pascal", "Pascal", drug_targets_db[i], target_df)
    j = j+1
}

## 2) Pascal UKBB
perc = 0.01
universe_genes = unique(pascal_UKBB_df$EnsemblId)
pascal_genes = unique(pascal_UKBB_df[pascal_UKBB_df$p_value <= quantile(pascal_UKBB_df$p_value, perc), "EnsemblId"])

for (i in 1:5){
    results_df = enrichment(pascal_genes, drug_targets[[i]], universe_genes, perc, "Pascal_UKBB", "Pascal_UKBB", drug_targets_db[i], results_df, j)
    results_df = AUC(pascal_UKBB_df, drug_targets[[i]], results_df, j)
    target_df = drug_target_info(pascal_genes, pascal_UKBB_df, drug_targets[[i]], drug_info[[i]], perc, "Pascal_UKBB", "Pascal_UKBB", drug_targets_db[i], target_df)
    j = j+1
}

## 3) Pascal - MR universe (GTEx & eQTLGen)
perc = 0.01
universe_genes = intersect(pascal_df$EnsemblId, mr_df$EnsemblId)
pascal_mr_df = pascal_df[pascal_df$EnsemblId %in% universe_genes,]
pascal_genes = unique(pascal_mr_df[pascal_mr_df$p_value <= quantile(pascal_mr_df$p_value, perc), "EnsemblId"])

for (i in 1:5){
    results_df = enrichment(pascal_genes, drug_targets[[i]], universe_genes, perc, "Pascal", "Pascal_eqtlgen_gtex", drug_targets_db[i], results_df, j)
    results_df = AUC(pascal_mr_df, drug_targets[[i]], results_df, j)
    target_df = drug_target_info(pascal_genes, pascal_mr_df, drug_targets[[i]], drug_info[[i]], perc, "Pascal", "Pascal_eqtlgen_gtex", drug_targets_db[i], target_df)
    j = j+1
}

## 4) Pascal - MR universe (deCODE proteins)
perc = 0.05
universe_genes = intersect(pascal_df$EnsemblId, protein_df$EnsemblId)
pascal_mr_df = pascal_df[pascal_df$EnsemblId %in% universe_genes,]
pascal_genes = unique(pascal_mr_df[pascal_mr_df$p_value <= quantile(pascal_mr_df$p_value, perc), "EnsemblId"])

for (i in 1:5){
    results_df = enrichment(pascal_genes, drug_targets[[i]], universe_genes, perc, "Pascal", "Pascal_decode", drug_targets_db[i], results_df, j)
    results_df = AUC(pascal_mr_df, drug_targets[[i]], results_df, j)
    target_df = drug_target_info(pascal_genes, pascal_mr_df, drug_targets[[i]], drug_info[[i]], perc, "Pascal", "Pascal_decode", drug_targets_db[i], target_df)
    j = j+1
}

## 5) MR - eQTLGen & GTEx
perc = 0.01
universe_genes = unique(mr_df$EnsemblId)
mr_genes = unique(mr_df[mr_df$p_value <= quantile(mr_df$p_value, perc), "EnsemblId"])

for (i in 1:5){
    results_df = enrichment(mr_genes, drug_targets[[i]], universe_genes, perc, "MR_eqtlgen_gtex", "eqtlgen_gtex", drug_targets_db[i], results_df, j)
    results_df = AUC(mr_df, drug_targets[[i]], results_df, j)
    target_df = drug_target_info(mr_genes, mr_df, drug_targets[[i]], drug_info[[i]], perc, "MR_eqtlgen_gtex", "eqtlgen_gtex", drug_targets_db[i], target_df)
    j = j+1
}

## 6) MR - eQTLGen & GTEx - UKBB
perc = 0.01
universe_genes = unique(mr_UKBB_df$EnsemblId)
mr_genes = unique(mr_UKBB_df[mr_UKBB_df$p_value <= quantile(mr_UKBB_df$p_value, perc), "EnsemblId"])

for (i in 1:5){
    results_df = enrichment(mr_genes, drug_targets[[i]], universe_genes, perc, "MR_eqtlgen_gtex_UKBB", "eqtlgen_gtex_UKBB", drug_targets_db[i], results_df, j)
    results_df = AUC(mr_UKBB_df, drug_targets[[i]], results_df, j)
    target_df = drug_target_info(mr_genes, mr_UKBB_df, drug_targets[[i]], drug_info[[i]], perc, "MR_eqtlgen_gtex_UKBB", "eqtlgen_gtex_UKBB", drug_targets_db[i], target_df)
    j = j+1
}

## 7) MR - eQTLGen
perc = 0.01
universe_genes = unique(mr_eqtlgen_df$EnsemblId)
mr_genes = unique(mr_eqtlgen_df[mr_eqtlgen_df$p_value <= quantile(mr_eqtlgen_df$p_value, perc), "EnsemblId"])

for (i in 1:5){
    results_df = enrichment(mr_genes, drug_targets[[i]], universe_genes, perc, "MR_eqtlgen", "eqtlgen", drug_targets_db[i], results_df, j)
    results_df = AUC(mr_eqtlgen_df, drug_targets[[i]], results_df, j)
    target_df = drug_target_info(mr_genes, mr_eqtlgen_df, drug_targets[[i]], drug_info[[i]], perc, "MR_eqtlgen", "eqtlgen", drug_targets_db[i], target_df)
    j = j+1
}

## 8) MR - GTEx
perc = 0.01
universe_genes = unique(mr_gtex_df$EnsemblId)
mr_genes = unique(mr_gtex_df[mr_gtex_df$p_value <= quantile(mr_gtex_df$p_value, perc), "EnsemblId"])

for (i in 1:5){
    results_df = enrichment(mr_genes, drug_targets[[i]], universe_genes, perc, "MR_gtex", "gtex", drug_targets_db[i], results_df, j)
    results_df = AUC(mr_gtex_df, drug_targets[[i]], results_df, j)
    target_df = drug_target_info(mr_genes, mr_gtex_df, drug_targets[[i]], drug_info[[i]], perc, "MR_gtex", "gtex", drug_targets_db[i], target_df)
    j = j+1
}

## 9) MR - eQTLGen (decode universe)
perc = 0.05
universe_genes = intersect(mr_eqtlgen_df$EnsemblId, protein_df$EnsemblId)
mr_eqtlgen_decode_df = mr_eqtlgen_df[mr_eqtlgen_df$EnsemblId %in% universe_genes,]
mr_genes = unique(mr_eqtlgen_decode_df[mr_eqtlgen_decode_df$p_value <= quantile(mr_eqtlgen_decode_df$p_value, perc), "EnsemblId"])

for (i in 1:5){
    results_df = enrichment(mr_genes, drug_targets[[i]], universe_genes, perc, "MR_eqtlgen", "eqtlgen_decode", drug_targets_db[i], results_df, j)
    results_df = AUC(mr_eqtlgen_decode_df, drug_targets[[i]], results_df, j)
    target_df = drug_target_info(mr_genes, mr_eqtlgen_decode_df, drug_targets[[i]], drug_info[[i]], perc, "MR_eqtlgen", "eqtlgen_decode", drug_targets_db[i], target_df)
    j = j+1
}

## 10) MR - eQTLGen & GTEx (decode universe)
perc = 0.05
universe_genes = intersect(mr_df$EnsemblId, protein_df$EnsemblId)
mr_decode_df = mr_df[mr_df$EnsemblId %in% universe_genes,]
mr_genes = unique(mr_decode_df[mr_decode_df$p_value <= quantile(mr_decode_df$p_value, perc), "EnsemblId"])

for (i in 1:5){
    results_df = enrichment(mr_genes, drug_targets[[i]], universe_genes, perc, "MR_eqtlgen_gtex", "eqtlgen_gtex_decode", drug_targets_db[i], results_df, j)
    results_df = AUC(mr_decode_df, drug_targets[[i]], results_df, j)
    target_df = drug_target_info(mr_genes, mr_decode_df, drug_targets[[i]], drug_info[[i]], perc, "MR_eqtlgen_gtex", "eqtlgen_gtex_decode", drug_targets_db[i], target_df)
    j = j+1
}

## 11) MR - decode proteins
perc = 0.05
universe_genes = unique(protein_df$EnsemblId)
mr_genes = unique(protein_df[protein_df$p_value <= quantile(protein_df$p_value, perc), "EnsemblId"])

for (i in 1:5){
    results_df = enrichment(mr_genes, drug_targets[[i]], universe_genes, perc, "MR_decode", "decode", drug_targets_db[i], results_df, j)
    results_df = AUC(protein_df, drug_targets[[i]], results_df, j)
    target_df = drug_target_info(mr_genes, protein_df, drug_targets[[i]], drug_info[[i]], perc, "MR_decode", "decode", drug_targets_db[i], target_df)
    j = j+1
}

## 12) MR - decode proteins (eQTLGen & GTEx universe)
perc = 0.05
universe_genes = intersect(mr_df$EnsemblId, protein_df$EnsemblId)
protein_eqtlgen_df = protein_df[protein_df$EnsemblId %in% universe_genes,]
mr_genes = unique(protein_eqtlgen_df[protein_eqtlgen_df$p_value <= quantile(protein_eqtlgen_df$p_value, perc), "EnsemblId"])

for (i in 1:5){
    results_df = enrichment(mr_genes, drug_targets[[i]], universe_genes, perc, "MR_decode", "eqtlgen_gtex_decode", drug_targets_db[i], results_df, j)
    results_df = AUC(protein_eqtlgen_df, drug_targets[[i]], results_df, j)
    target_df = drug_target_info(mr_genes, protein_eqtlgen_df, drug_targets[[i]], drug_info[[i]], perc, "MR_decode", "eqtlgen_gtex_decode", drug_targets_db[i], target_df)
    j = j+1
}

## 13) MR - decode proteins - UKBB
perc = 0.05
universe_genes = unique(protein_ukbb_df$EnsemblId)
mr_genes = unique(protein_ukbb_df[protein_ukbb_df$p_value <= quantile(protein_ukbb_df$p_value, perc), "EnsemblId"])

for (i in 1:5){
    results_df = enrichment(mr_genes, drug_targets[[i]], universe_genes, perc, "MR_decode_UKBB", "decode_UKBB", drug_targets_db[i], results_df, j)
    results_df = AUC(protein_ukbb_df, drug_targets[[i]], results_df, j)
    target_df = drug_target_info(mr_genes, protein_ukbb_df, drug_targets[[i]], drug_info[[i]], perc, "MR_decode_UKBB", "decode_UKBB", drug_targets_db[i], target_df)
    j = j+1
}

## 14) Exome
perc = 0.01
universe_genes = unique(exome_df$EnsemblId)
exome_genes = unique(exome_df[exome_df$p_value <= quantile(exome_df$p_value, perc), "EnsemblId"])

for (i in 1:5){
    results_df = enrichment(exome_genes, drug_targets[[i]], universe_genes, perc, "Exome", "Exome", drug_targets_db[i], results_df, j)
    results_df = AUC(exome_df, drug_targets[[i]], results_df, j)
    target_df = drug_target_info(exome_genes, exome_df, drug_targets[[i]], drug_info[[i]], perc, "Exome", "Exome", drug_targets_db[i], target_df)
    j = j+1
}

## 15) Random
perc = 0.01
universe_genes = unique(random_df$EnsemblId)
random_genes = unique(random_df[random_df$p_value <= quantile(random_df$p_value, perc), "EnsemblId"])

for (i in 1:5){
    results_df = enrichment(random_genes, drug_targets[[i]], universe_genes, perc, "Random", "All", drug_targets_db[i], results_df, j)
    results_df = AUC(random_df, drug_targets[[i]], results_df, j)
    j = j+1
}

write.table(results_df, snakemake@output[[1]], row.names = F, quote = F, sep = '\t')

# delete first line of target_df (line of NAs)
target_df = target_df[-1, ]
write.table(target_df, snakemake@output[[2]], row.names = F, quote = F, sep = '\t')