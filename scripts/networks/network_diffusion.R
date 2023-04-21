library(igraph)
library(Rcpp)
library(stringr)
library(pROC)

#### Network diffusion utilities and parameters ####
sourceCpp(paste0(getwd(), "/scripts/networks/mat_util.cpp"))

print("Get network parameters")
r = as.numeric(snakemake@wildcards[["r"]])/10
seed = as.integer(snakemake@params[["seed"]])
print(paste0("Seed: ", seed))
weighted = snakemake@wildcards[["weighted"]]
t = as.numeric(snakemake@wildcards[["thresh"]])

#### Network transition matrix ####
print("Calculate network")
edge_df = read.table(snakemake@input[["nw"]], header = T)
edge_df$weight = abs(edge_df$weight)
edge_df = edge_df[edge_df$weight >= t, ]

if (weighted == "no"){ 
    edge_df$weight = 1
}

graph_full = graph_from_data_frame(edge_df, directed = FALSE, vertices = NULL)

# extract adjacency matrix of biggest component
components <- clusters(graph_full, mode = "weak")
biggest_cluster_id <- which.max(components$csize)
vert_ids <- V(graph_full)[components$membership == biggest_cluster_id]

if (weighted == "no"){
    graph <- as_adjacency_matrix(induced_subgraph(graph_full, vert_ids), sparse = F)
    gene_order <- rownames(graph)
} else {
    vert_ids_names = as.vector(names(vert_ids))
    edge_df = edge_df[(edge_df$gene1 %in% vert_ids_names) & (edge_df$gene2 %in% vert_ids_names), ]

    graph <- as_adjacency_matrix(induced_subgraph(graph_full, vert_ids), sparse = F)
    gene_order <- rownames(graph)
    gene_order_df = data.frame(gene = gene_order, idx_R = seq(1:(length(gene_order))))
    gene_order_df$idx = gene_order_df$idx_R - 1

    edge_df = merge(edge_df, gene_order_df[, c("gene", "idx")], by.x = "gene1", by.y = "gene")
    edge_df = edge_df[, c("idx", "gene2", "weight")]
    names(edge_df) = c("idx1", "gene2", "weight")

    edge_df = merge(edge_df, gene_order_df[, c("gene", "idx")], by.x = "gene2", by.y = "gene")
    edge_df = edge_df[, c("idx1", "idx", "weight")]
    names(edge_df) = c("idx1", "idx2", "weight")

    graph = add_weights_adjacency(graph, edge_df$idx1, edge_df$idx2, edge_df$weight)

}

print("Biggest component is extracted")
print(graph[1:5,1:5])

# RWR
W = stoch_col_norm_(graph)
T = solve(diag(length(gene_order))-(1-r)*W) # transition matrix
print("Transition matrix is calculated")

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

#### MR - Protein ####
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

#### Random ####
set.seed(seed)
random_scores <- rnorm(length(gene_order))
random_df <- data.frame(EnsemblId = gene_order,
                         z = random_scores)
print(head(random_df))

#### Result dataframes & Enrichment function ####
diffusion_df <- data.frame(EnsemblId = gene_order) # store diffusion profiles

colnames = c("Percentile", "Method", "Universe", "Database", "TP", "OR", "Pval", "diagonal_padding", "ci_low", "ci_up", "AUC", "ci_low_AUC", "ci_up_AUC")
results_df = data.frame(matrix(NA, n_databases*8, length(colnames)))
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

#### Calculate diffusion profiles and enrichments ####

diffusion <- function(method_df, gene_order, T, z = FALSE){
    
    # set missing method genes which are present in network to p-value of 0
    if (z == FALSE){
        method_df = method_df[method_df$EnsemblId %in% gene_order, ]
        missing_genes <- setdiff(gene_order, method_df$EnsemblId)
        missing_df <- data.frame(EnsemblId = missing_genes,
                                p_value = rep(1, length(missing_genes)))
        method_df <- rbind(method_df, missing_df)
        method_df <- method_df[match(gene_order, method_df$EnsemblId), ] # order has to match graph order!
        method_df[method_df$p_value == 0, "p_value"] = 1e-300
        method_df$z2 <- qnorm(method_df$p_value / 2, lower.tail = F)**2
    } else {
        method_df = method_df[method_df$EnsemblId %in% gene_order, ]
        missing_genes <- setdiff(gene_order, method_df$EnsemblId)
        missing_df <- data.frame(EnsemblId = missing_genes,
                                 z = rep(0, length(missing_genes)))
        method_df <- rbind(method_df, missing_df)
        method_df <- method_df[match(gene_order, method_df$EnsemblId), ] # order has to match graph order!
        method_df$z2 <- method_df$z**2
    }

    # diffuse through network
    p0 <- method_df$z2 / sum(method_df$z2)
    d <- T %*% p0

    return (d)
}

j = 1

## 1) Pascal
print("Calculate Pascal diffusion profile")
d = diffusion(pascal_df, gene_order, T)
diffusion_df["Pascal"] = d

perc = 0.99
pascal_genes = unique(diffusion_df[diffusion_df$Pascal >= quantile(diffusion_df$Pascal, perc), "EnsemblId"])
diffusion_df$rank = rank(-diffusion_df$Pascal, ties.method = "min")

for (i in 1:5){
    results_df = enrichment(pascal_genes, drug_targets[[i]], gene_order, perc, "Pascal", "Pascal", drug_targets_db[i], results_df, j)
    results_df = AUC(diffusion_df, drug_targets[[i]], results_df, j)
    target_df = drug_target_info(pascal_genes, diffusion_df, drug_targets[[i]], drug_info[[i]], perc, "Pascal", "Pascal", drug_targets_db[i], target_df)
    j = j+1
}

## 2) Pascal - UKBB
print("Calculate Pascal diffusion profile - UKBB")
d = diffusion(pascal_UKBB_df, gene_order, T)
diffusion_df["Pascal_UKBB"] = d

perc = 0.99
pascal_genes = unique(diffusion_df[diffusion_df$Pascal_UKBB >= quantile(diffusion_df$Pascal_UKBB, perc), "EnsemblId"])
diffusion_df$rank = rank(-diffusion_df$Pascal_UKBB, ties.method = "min")

for (i in 1:5){
    results_df = enrichment(pascal_genes, drug_targets[[i]], gene_order, perc, "Pascal_UKBB", "Pascal_UKBB", drug_targets_db[i], results_df, j)
    results_df = AUC(diffusion_df, drug_targets[[i]], results_df, j)
    target_df = drug_target_info(pascal_genes, diffusion_df, drug_targets[[i]], drug_info[[i]], perc, "Pascal_UKBB", "Pascal_UKBB", drug_targets_db[i], target_df)
    j = j+1
}

## 3) MR - eQTLGen & GTEx
print("Calculate MR - eQTLGen & GTEx diffusion profile")
d = diffusion(mr_df, gene_order, T)
diffusion_df["MR_eqtlgen_gtex"] = d

perc = 0.99
mr_genes = unique(diffusion_df[diffusion_df$MR_eqtlgen_gtex >= quantile(diffusion_df$MR_eqtlgen_gtex, perc), "EnsemblId"])
diffusion_df$rank = rank(-diffusion_df$MR_eqtlgen_gtex, ties.method = "min")

for (i in 1:5){
    results_df = enrichment(mr_genes, drug_targets[[i]], gene_order, perc, "MR_eqtlgen_gtex", "eqtlgen_gtex", drug_targets_db[i], results_df, j)
    results_df = AUC(diffusion_df, drug_targets[[i]], results_df, j)
    target_df = drug_target_info(mr_genes, diffusion_df, drug_targets[[i]], drug_info[[i]], perc, "MR_eqtlgen_gtex", "eqtlgen_gtex", drug_targets_db[i], target_df)
    j = j+1
}

## 4) MR - eQTLGen & GTEx - UKBB
print("Calculate MR - eQTLGen & GTEx diffusion profile - UKBB")
d = diffusion(mr_UKBB_df, gene_order, T)
diffusion_df["MR_eqtlgen_gtex_UKBB"] = d

perc = 0.99
mr_genes = unique(diffusion_df[diffusion_df$MR_eqtlgen_gtex_UKBB >= quantile(diffusion_df$MR_eqtlgen_gtex_UKBB, perc), "EnsemblId"])
diffusion_df$rank = rank(-diffusion_df$MR_eqtlgen_gtex_UKBB, ties.method = "min")

for (i in 1:5){
    results_df = enrichment(mr_genes, drug_targets[[i]], gene_order, perc, "MR_eqtlgen_gtex_UKBB", "eqtlgen_gtex_UKBB", drug_targets_db[i], results_df, j)
    results_df = AUC(diffusion_df, drug_targets[[i]], results_df, j)
    target_df = drug_target_info(mr_genes, diffusion_df, drug_targets[[i]], drug_info[[i]], perc, "MR_eqtlgen_gtex_UKBB", "eqtlgen_gtex_UKBB", drug_targets_db[i], target_df)
    j = j+1
}

## 5) MR - eQTLGen
print("Calculate MR - eQTLGen diffusion profile")
d = diffusion(mr_eqtlgen_df, gene_order, T)
diffusion_df["MR_eqtlgen"] = d

perc = 0.99
mr_genes = unique(diffusion_df[diffusion_df$MR_eqtlgen >= quantile(diffusion_df$MR_eqtlgen, perc), "EnsemblId"])
diffusion_df$rank = rank(-diffusion_df$MR_eqtlgen, ties.method = "min")

for (i in 1:5){
    results_df = enrichment(mr_genes, drug_targets[[i]], gene_order, perc, "MR_eqtlgen", "eqtlgen", drug_targets_db[i], results_df, j)
    results_df = AUC(diffusion_df, drug_targets[[i]], results_df, j)
    target_df = drug_target_info(mr_genes, diffusion_df, drug_targets[[i]], drug_info[[i]], perc, "MR_eqtlgen", "eqtlgen", drug_targets_db[i], target_df)
    j = j+1
}

## 6) MR - decode proteins
print("Calculate MR - decode diffusion profile")
d = diffusion(protein_df, gene_order, T)
diffusion_df["MR_decode"] = d

perc = 0.99
mr_genes = unique(diffusion_df[diffusion_df$MR_decode >= quantile(diffusion_df$MR_decode, perc), "EnsemblId"])
diffusion_df$rank = rank(-diffusion_df$MR_decode, ties.method = "min")

for (i in 1:5){
    results_df = enrichment(mr_genes, drug_targets[[i]], gene_order, perc, "MR_decode", "decode", drug_targets_db[i], results_df, j)
    results_df = AUC(diffusion_df, drug_targets[[i]], results_df, j)
    target_df = drug_target_info(mr_genes, diffusion_df, drug_targets[[i]], drug_info[[i]], perc, "MR_decode", "decode", drug_targets_db[i], target_df)
    j = j+1
}

## 7) Exome
print("Calculate Exome diffusion profile")
d = diffusion(exome_df, gene_order, T)
diffusion_df["Exome"] = d

perc = 0.99
exome_genes = unique(diffusion_df[diffusion_df$Exome >= quantile(diffusion_df$Exome, perc), "EnsemblId"])
diffusion_df$rank = rank(-diffusion_df$Exome, ties.method = "min")

for (i in 1:5){
    results_df = enrichment(exome_genes, drug_targets[[i]], gene_order, perc, "Exome", "Exome", drug_targets_db[i], results_df, j)
    results_df = AUC(diffusion_df, drug_targets[[i]], results_df, j)
    target_df = drug_target_info(exome_genes, diffusion_df, drug_targets[[i]], drug_info[[i]], perc, "Exome", "Exome", drug_targets_db[i], target_df)
    j = j+1
}

## 8) Random
print("Calculate Random diffusion profile")
d = diffusion(random_df, gene_order, T, z = TRUE)
diffusion_df["Random"] = d

perc = 0.99
random_genes = unique(diffusion_df[diffusion_df$Random >= quantile(diffusion_df$Random, perc), "EnsemblId"])
diffusion_df$rank = rank(-diffusion_df$Random, ties.method = "min")

for (i in 1:5){
    results_df = enrichment(random_genes, drug_targets[[i]], gene_order, perc, "Random", "All", drug_targets_db[i], results_df, j)
    results_df = AUC(diffusion_df, drug_targets[[i]], results_df, j)
    j = j+1
}

write.table(diffusion_df, snakemake@output[["diffusion"]], row.names = F, quote = F, sep = '\t')
write.table(results_df, snakemake@output[["overlap"]], row.names = F, quote = F, sep = '\t')

target_df = target_df[-1, ]
write.table(target_df, snakemake@output[["target_info"]], row.names = F, quote = F, sep = '\t')