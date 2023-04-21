library(biomaRt)

ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
genes <- getBM(attributes=c('ensembl_gene_id',
                            'ensembl_peptide_id','hgnc_symbol'), mart = ensembl)

dim(genes) # 159680
head(genes)
genes = genes[genes$ensembl_peptide_id != "",]
dim(genes) # 104770

length(unique(genes$ensembl_gene_id)) # 23393
length(unique(genes$ensembl_peptide_id)) # 104763
length(unique(genes$hgnc_symbol)) # 19684

write.table(genes, "Ensembl_transcript_peptide_HGNC.tsv", sep = '\t', row.names = F, quote = F)