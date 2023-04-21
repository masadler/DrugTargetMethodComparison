library(igraph)
library(Rcpp)
library(stringr)
library(pROC)

sourceCpp(paste0(getwd(), "/scripts/networks/mat_util.cpp"))

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

graph <- graph_from_data_frame(edge_df, directed = FALSE, vertices = NULL)

degrees <- strength(graph) # weighted vertex degree

universe_genes <- names(V(graph))
diffusion_df <- data.frame(EnsemblId = universe_genes, degree = degrees)


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

#### Enrichment with node degree ####

colnames = c("Percentile", "Method", "Universe", "Database", "TP", "OR", "Pval", "diagonal_padding", "ci_low", "ci_up", "AUC", "ci_low_AUC", "ci_up_AUC")
results_df = data.frame(matrix(NA, n_databases, length(colnames)))
names(results_df) = colnames


enrichment <- function(method_genes, drug_genes, universe_genes, perc, method, universe, database, results_df, i){
        
        # intersect with universe
        method_genes = intersect(method_genes, universe_genes)
        drug_genes = intersect(drug_genes, universe_genes)

        # calculate enrichment  
        diagonal_padding = FALSE

        # calculate enrichment        
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

    method_df$rank = rank(-method_df$degree, ties.method = "min")
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

## Degree
j = 1
perc = 0.99
degree_genes = unique(diffusion_df[diffusion_df$degree >= quantile(diffusion_df$degree, perc), "EnsemblId"])

for (i in 1:5){
    results_df = enrichment(degree_genes, drug_targets[[i]], universe_genes, perc, "Degree", "All", drug_targets_db[i], results_df, j)
    results_df = AUC(diffusion_df, drug_targets[[i]], results_df, j)
    j = j+1
}

write.table(diffusion_df, snakemake@output[["diffusion"]], row.names = F, quote = F, sep = '\t')
write.table(results_df, snakemake@output[["overlap"]], row.names = F, quote = F, sep = '\t')