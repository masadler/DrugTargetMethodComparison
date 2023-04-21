# Gene prioritization methods

Welcome to the Github accompanying the pre-print [**Multi-layered genetic approaches to identify approved drug targets**](https://doi.org/10.1101/2023.03.21.23285637) which compares gene prioritization methods in their ability to identify historical drug targets.

# Citations

If you use scripts from this Github please consider citing the pre-print:

Sadler MC, Auwerx C, Deelen P, Kutalik Z. Multi-layered genetic approaches to identify approved drug targets. [medRxiv. 2023:2023-03.](https://doi.org/10.1101/2023.03.21.23285637)

# Usage

This Github contains the workflow pipeline that was used to define drug targets, calculate method-specific gene scores, network diffusion scores and drug target enrichment scores (odds ratios and AUC scores). Code to calculate other statistics presented in the study such as agreement between different gene scoring methods, drug database and network statistics as well drug target heritability analyses is also part of the pipeline.

The workflow is organized as a [snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline. This allows for a reproducible analysis and parallel computations - this is especially useful since the analysis has been conducted on 30 different traits/drug indications and 3 different networks. However, individual scripts can also be run without the snakemake workflow manager by replacing snakemake variables by hard-coded input and output paths.

# Input data

In this study, input data comes from various public databases which can be large datasets and sometimes needs the signing of terms and conditions (i.e. to download the deCODE data). Here, we describe the different steps that need to be taken in order to assemble all required input data. When using this data, please consider citing the sources appropriately.

## Drug data

To define drug targets two different sources of input data is needed:

1) Indication - drug links
2) Drug - gene target links

Indication - drug links have been collected from [Ruiz et al., Nature Communications, 2021](https://doi.org/10.1038/s41467-021-21770-8), [DrugBank](https://go.drugbank.com) and [ChEMBL](https://www.ebi.ac.uk/chembl/). For Ruiz et al., Supplementary Table 6 needs to be downloaded, for DrugBank, the file (`scripts/drug_targets/GWAS_x_drug_drugBank.py`) should be run and for ChEMBL the indication - drug link file can be downloaded [here](https://www.ebi.ac.uk/chembl/g/#browse/drug_indications) (click on "Select All" and then on the right "TSV" button to download the full file).

Drug - gene target links have been collected from [DGIdb](https://www.dgidb.org), [STITCH](http://stitch.embl.de) and [ChEMBL](https://www.ebi.ac.uk/chembl/). On the DGIdb website, drug-target links can be downloaded from the [Downloads page](https://www.dgidb.org/downloads) (interactions, genes and drugs files), STITCH data can be downloaded from the [Download page](http://stitch.embl.de/cgi/download.pl?UserId=Ov0ZH2xj6cD6&sessionId=lre85GJq7ib4) (select "Homo sapiens" organism and then download the `9606.actions.v5.0.tsv.gz` file as well as the `chemical.sources.v5.0.tsv.gz` file which contains the mapping of chemical IDs to DrugBank - since this file is very large, this pipeline assumes that the file is already subsetted to entries relevant to DrugBank and named `stitch_chemical_drugbank_vocabulary.tsv`)

## Gene prioritization data

1) GWAS scores

GWAS scores should be computed using the [PascalX software](https://bergmannlab.github.io/PascalX/index.html). The documentation on installing the software is very detailed. With default settings, gene scores are calculated based on gene names and not EnsemblIds. The rule `run_pascal_ensembl` in the Snakefile calculates gene scores based on EnsemblIds. Downloading and formatting the reference is explained in the [documentation](https://bergmannlab.github.io/PascalX/usage.html).

In this pipeline, GWAS summary statistics are assumed to be in the [GCTA-COJO format](https://cnsgenomics.com/software/gcta/#COJO) (.ma), however, any other format works as well with the PascalX software as long as the field number containing p-values is specified.

2) QTL-GWAS scores

QTL-GWAS scores have been calculated using the [smr-ivw software](https://github.com/masadler/smrivw). The [Github wiki](https://github.com/masadler/smrivw/wiki/Univariable-MR) on univariable MR explains usage and theory behind MR gene causal effect calculations, as well as the formatting of input QTL data. For this software, GWAS summary statistics have to be in [GCTA-COJO format](https://cnsgenomics.com/software/gcta/#COJO) (.ma). The rules `transcript_eQTLGen_MR` and `Decode_MR` in the Snakefile calculate transcript and protein QTL-GWAS scores, respectively. A reference panel in [PLINK 1.9 Format](https://www.cog-genomics.org/plink/2.0/input#pheno) (.bed, .bim, .fam files) is required.

3) Exome scores

In this study, pre-calculated Exome scores from [Backman et al., Nature, 2021](https://doi.org/10.1038/s41586-021-04103-z) were used. All scores can be calculated from the [GWAS catalog](https://www.ebi.ac.uk/gwas/) and the full list of mapping between GWAS Catalog accession IDs and UK Biobank traits is available in the Supplementary Data 4 of the [Backman et al., Nature, 2021](https://doi.org/10.1038/s41586-021-04103-z) paper. Files downloaded for this study are specified in `data/Exome_data.csv`.

## Networks

Three networks have been used in this study:

1) STRING: The STRING protein-protein interaction network v11 can be downloaded [here](https://string-db.org/cgi/download?sessionId=brnMKItddrXp): select "Homo sapiens" organism and then download the `9606.protein.links.v11.5.txt.gz` file.

2) CoXRNAseq: The co-expression network calculated from RNA-seq samples can be downloaded [here](https://downloads.molgeniscloud.org/downloads/depict2/) (documentation can be found [here](https://github.com/molgenis/systemsgenetics/wiki/Downstreamer)). Only the network `gene_coregulation_165_eigenvectors_protein_coding.txt` from the bundle is needed.

3) FAVA: This network which is based on single cell RNA-seq read-count data from the Human Protein Atlas and proteomics data from the PRoteomics IDEntifications (PRIDE) database can be downloaded [here](https://doi.org/10.5281/zenodo.6803472).

## Heritability data

Cis-heritability estimates for transcripts and proteins are available in the Supplementary tables of this study. Alternatively, they can be calculated with the [LDAK software](https://dougspeed.com) using the following command:

```bash
ldak5.2.linux --reml ldak.reml --summary sumstat --bfile refpanel --power -.25 --ignore-weights YES --region-number 1 --region-prefix region_prefix
```

where ldak.reml is the prefix of your output file, sumstat the summary statistics in LDAK format, refpanel the reference panel in Plink format, region_prefix the prefix of the file containing your SNPs from the region. Further details are in the documentation of the [LDAK REML](https://dougspeed.com/reml-analysis/) method.

## License

MIT License
