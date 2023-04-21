library(igraph)
library(stringr)

colnames = c("Network", "Nodes", "Edges", "Median_degree", "Average_log_degree", "sd_log_degree")
result_df = data.frame(matrix(NA, 3, length(colnames)))
names(result_df) = colnames

prepare_network <- function(file){

    edge_df = read.table(file, header = T)
    edge_df$weight = abs(edge_df$weight)

    graph_full = graph_from_data_frame(edge_df, directed = FALSE, vertices = NULL)

    # extract adjacency matrix of biggest component
    components <- clusters(graph_full, mode = "weak")
    biggest_cluster_id <- which.max(components$csize)
    vert_ids <- V(graph_full)[components$membership == biggest_cluster_id]

    vert_ids_names = as.vector(names(vert_ids))
    edge_df = edge_df[(edge_df$gene1 %in% vert_ids_names) & (edge_df$gene2 %in% vert_ids_names), ]

    graph <- induced_subgraph(graph_full, vert_ids)

    return (graph)
}

calculate_properties <- function(graph, result_df, i){
    # number of nodes

    result_df[i, "Nodes"] = length(V(graph))
    result_df[i, "Edges"] = gsize(graph)

    # degree distribution
    degrees = strength(graph)
    result_df[i, "Median_degree"] = median(degrees)
    result_df[i, "Average_log_degree"] = mean(log(degrees))
    result_df[i, "sd_log_degree"] = sd(log(degrees))

    return (result_df)
}

# STRING
print("STRING")
i = 1
file = snakemake@input[["nw_files"]][i]
result_df[i, "Network"] = "STRING"
graph = prepare_network(file)
result_df = calculate_properties(graph, result_df, i)


# FAVA
print("FAVA")
i = 2
file = snakemake@input[["nw_files"]][i]
result_df[i, "Network"] = "FAVA"
graph = prepare_network(file)
result_df = calculate_properties(graph, result_df, i)

# CoXRNAseq
print("CoXRNAseq")
i = 3
file = snakemake@input[["nw_files"]][i]
result_df[i, "Network"] = "CoXRNAseq"
graph = prepare_network(file)
result_df = calculate_properties(graph, result_df, i)

result_df = t(result_df)
write.table(result_df, snakemake@output[[1]], quote = F, sep = '\t')