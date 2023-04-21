import pandas as pd

df = pd.read_csv("data/GWAS.csv")
traits = df.Acronym
df_UKBB = pd.read_csv("data/GWAS_UKBB.csv")
traits_UKBB = df_UKBB.Acronym_UKBB
traits_UKBB_dic = df_UKBB.set_index("Acronym")["Acronym_UKBB"].to_dict()
df_Exome = pd.read_csv("data/Exome_data.csv")
traits_exome_dic = df_Exome.set_index("Acronym")["Acronym_Exome"].to_dict()

networks = ["STRING", "FAVA", "CoXRNAseqf"]
rs = ["2", "4", "6", "8"]
seed_dic = {trait: str(i) for i, trait in enumerate(traits)}

eQTLGen_data_path = "" # Path to eQTLGen QTL data in BESD format
deCODE_data_path = "" # Path to deCODE QTL data in BESD format
GTEx_MR_result_path = "" # Path do computed GTEx MR results in .msmr format
Pascal_result_path = "" # Path to computed Pascal GWAS scores
GWAS_path = "" # Path to GWAS summary statistics in .ma format
exome_path = "" # Path to Exome burden test summary statistics from Backman et al. 2021
ref_panel_path = "" # Path to a reference panel in Plink format, e.g. UK10K, 1000G
ref_panel_path_Pascal = "" # Path to reference panel for PascalX input: see documentation here: https://bergmannlab.github.io/PascalX/usage.html
genome_annotation_file = "" # Ensembl protein coding transcripts -> run scripts/method/get_ensembl_annotation_file.py
biomart_mapping_file = "" # Biomart mapping of transcript Ensembl, protein Ensembl and HGNC gene symbols -> run scripts/method/biomart_mapping.R
ruiz_path = "" # Path to supplementary tables from Ruiz et al. 2021, Nature Communications
chembl_path = "" # Path to ChEMBL files
dgidb_path = "" # Path to DGIdb files
stitch_path = "" # Path to STITCH files
nw_path = "" # Path to original networks
reml_path = "" # Path to REML-LDAK cis-heritability files

rule all:
    input:
        "output/method/pascal_mr_enrichment.tsv",
        "output/method/pascal_exome_enrichment.tsv",
        "output/method/mr_exome_enrichment.tsv",
        "output/method/mr_tissue_statistics.tsv",
        "output/method/mr_eqtl_pqtl_enrichment.tsv",
        "output/drug_targets/drug_database_statistics.tsv",
        "output/drug_targets/drug_database_correlations.txt",
        "output/direct_overlap/direct_overlap_aggregated_enrichment.tsv",
        "output/direct_overlap/direct_overlap_aggregated_target_info.tsv",
        "output/direct_overlap/direct_overlap_aggregated_running_OR.tsv",
        "output/direct_overlap/direct_overlap_pairwise_drug_overlap_Jaccard.tsv",
        "output/direct_overlap/direct_overlap_pairwise_drug_overlap_Jaccard_ukbb.tsv",
        "output/direct_overlap/transcript_remlheritability_aggregated_drug_targets.tsv",
        "output/direct_overlap/protein_remlheritability_aggregated_drug_targets.tsv",
        "output/networks/network_statistics.tsv",
        expand("output/networks/diffusion_overlap_aggregated_enrichment.{nw}.r{r}.wyes.t0.tsv", r = rs, nw = networks),
        expand("output/networks/diffusion_overlap_aggregated_target_info.{nw}.r{r}.wyes.t0.tsv", r = rs, nw = networks),
        expand("output/networks/diffusion_overlap_aggregated_degree_enrichment.{nw}.wyes.t0.tsv", nw = networks)

rule run_pascal_ensembl:
    input:
        gwas = f'{GWAS_path}' + "{trait}_gwas_summary_uk10kck.ma",
        annotation = genome_annotation_file
    params:
        ref = ref_panel_path_Pascal
    output:
        f'{Pascal_result_path}' + "Pascal_Ensembl_scoring_{trait}.tsv"
    shell:
        """
        python3 scripts/method/run_pascal_ensembl.py -r {params.ref} -g {input.gwas} -a {input.annotation} -o {output}
        """

rule transcript_eQTLGen_MR:
    input:
        gwas = f'{GWAS_path}' + "{trait}_gwas_summary_uk10kck.ma",
        eqtl= f'{eQTLGen_data_path}' + "cis-eQTL.besd",
        uk10k= f'{ref_panel_path}' + "uk10k.autosomal.bed"
    output:
        "output/omics_MR/transcript_eQTLGen_{trait}_MR.msmr"
    params:
        out = "output/omics_MR/transcript_eQTLGen_{trait}_MR",
        eqtl= f'{eQTLGen_data_path}' + "cis-eQTL",
        uk10k = f'{ref_panel_path}' + "uk10k.autosomal"       
    threads: 4
    resources:
        mem_mb = 10000
    shell:
        "smrivw --bfile {params.uk10k} --gwas-summary {input.gwas} --beqtl-summary {params.eqtl} "
        "--out {params.out} --smr-ivw --peqtl-smr 1e-6 --ld-multi-snp 0.01 "
        "--thread-num {threads} --diff-freq-prop 0.3 --diff-freq 0.05 --get-snp-effects"

rule Decode_MR:
    input:
        gwas = f'{GWAS_path}' + "{trait}_gwas_summary_uk10kck.ma",
        pqtl= f'{deCODE_data_path}' + "cis-pQTL.besd",
        uk10k= f'{ref_panel_path}' + "uk10k.autosomal.bed"
    output:
        "output/omics_MR/protein_deCODE_{trait}_MR.msmr"
    params:
        out = "output/omics_MR/protein_deCODE_{trait}_MR",
        pqtl= f'{deCODE_data_path}' + "cis-pQTL",
        uk10k= f'{ref_panel_path}' + "uk10k.autosomal.bed"
    threads: 4
    resources:
        mem_mb = 10000
    shell:
        "smrivw --bfile {params.uk10k} --gwas-summary {input.gwas} --beqtl-summary {params.pqtl} "
        "--out {params.out} --smr-ivw --peqtl-smr 0.05 --ld-multi-snp 0.01 "
        "--thread-num {threads} --diff-freq-prop 0.3 --diff-freq 0.05 --get-snp-effects"

rule prepare_exome_data:
    input:
        trait_info = "data/Exome_data.csv"
    params:
        data_path = exome_path
    output:
        ["output/exome_genes/"+traits_exome_dic[trait]+"_gene_burden_M3_1_associations.tsv" for trait in traits]
    script:
        "scripts/method/prepare_exome_data.py"

rule aggregate_GTEx_eQTLGen_MR:
    input:
        "output/omics_MR/transcript_eQTLGen_{trait}_MR.msmr"
    params:
        path = GTEx_MR_result_path
    output:
        "output/omics_MR/{trait}_aggregated_gtex_eQTLGen_MR.tsv",
        "output/omics_MR/{trait}_aggregated_gtex_eQTLGen_MR_mintissue.tsv",
        "output/omics_MR/{trait}_aggregated_gtex_MR.tsv",
    script:
        "scripts/method/aggregate_GTEx_tissue_MR.py"

rule pascal_mr_enrichment:
    input:
        pascal = expand(f'{Pascal_result_path}' + "Pascal_Ensembl_scoring_{trait}.tsv", trait = traits),
        mr = expand("output/omics_MR/{trait}_aggregated_gtex_eQTLGen_MR.tsv", trait = traits),
        mr_gtex = expand("output/omics_MR/{trait}_aggregated_gtex_MR.tsv", trait = traits),
        gene_info = genome_annotation_file # contains only protein-coding genes
    output:
        "output/method/pascal_mr_enrichment.tsv"
    script:
        "scripts/method/pascal_mr_enrichment.R"

rule pascal_exome_enrichment:
    input:
        pascal = expand(f'{Pascal_result_path}' + "Pascal_Ensembl_scoring_{trait}.tsv", trait = traits),
        pascal_UKBB = expand(f'{Pascal_result_path}' + "Pascal_Ensembl_scoring_{trait}.tsv", trait = traits_UKBB),
        exome = ["output/exome_genes/" + traits_exome_dic[trait] + "_gene_burden_M3_1_associations.tsv" for trait in traits]
    output:
        "output/method/pascal_exome_enrichment.tsv"
    script:
        "scripts/method/pascal_exome_enrichment.R"

rule mr_exome_enrichment:
    input:
        mr = expand("output/omics_MR/{trait}_aggregated_gtex_eQTLGen_MR.tsv", trait = traits),
        mr_UKBB = expand("output/omics_MR/{trait}_aggregated_gtex_eQTLGen_MR.tsv", trait = traits_UKBB),
        exome = ["output/exome_genes/" + traits_exome_dic[trait] + "_gene_burden_M3_1_associations.tsv" for trait in traits]
    output:
        "output/method/mr_exome_enrichment.tsv"
    script:
        "scripts/method/mr_exome_enrichment.R"

rule mr_tissue_statistics:
    input:
        mr = expand("output/omics_MR/{trait}_aggregated_gtex_eQTLGen_MR.tsv", trait = traits),
        gene_info = genome_annotation_file,
        tissue = "data/Tissue_categories.csv"
    output:
        "output/method/mr_tissue_statistics.tsv"
    script:
        "scripts/method/mr_tissue_statistics.py"

rule mr_eqtl_pqtl_enrichment: 
    input:
        mr_e = expand("output/omics_MR/transcript_eQTLGen_{trait}_MR.msmr", trait = traits),
        mr_p = expand("output/omics_MR/protein_deCODE_{trait}_MR.msmr", trait = traits),
        protein_info = f'{deCODE_data_path}' + "Decode_Somamer_info_TSS_GRCh38.csv" # mapping of protein IDs to EnsemblIds
    output:
        "output/method/mr_eqtl_pqtl_enrichment.tsv"
    script:
        "scripts/method/mr_eqtl_pqtl_enrichment.R"

rule GWAS_x_drug_Ruiz:
    input:
        gwas_indication = "data/GWAS_indication.csv",
        drug_indication_data = f'{ruiz_path}' + "6_drug_indication_df_Ruiz_et_al_2021.csv"
    output:
        "output/drug_targets/{trait}_drugs_Ruiz.csv"
    script:
        "scripts/drug_targets/GWAS_x_drug_Ruiz.py"

rule GWAS_x_drug_Chembl:
    input:
        gwas_indication = "data/GWAS_indication.csv",
        drug_indication_data = f'{chembl_path} + "Chembl_drug_indications_20220516.tsv"
    output:
        "output/drug_targets/{trait}_drugs_Chembl.csv"
    script:
        "scripts/drug_targets/GWAS_x_drug_Chembl.py"

rule drug_Ruiz_x_gene_dgidb:
    input:
        drugs = "output/drug_targets/{trait}_drugs_Ruiz.csv",
        dgidb_data = f'{dgidb_path}' + "dgidb_interactions_202101.tsv",
        gene_voc = f'{dgidb_path}' + "dgidb_genes_202101.tsv",
        drug_voc = f'{dgidb_path}' + "dgidb_drugs_202101.tsv",
        drug_voc_chembl_db = "data/Chembl_to_drugbank_dgidb_202202.txt"
    output:
        "output/drug_targets/{trait}_drug_Ruiz_gene_dgidb.tsv"
    script:
        "scripts/drug_targets/drug_x_gene_dgidb.py"

rule drug_Ruiz_x_gene_stitch:
    input:
        drugs = "output/drug_targets/{trait}_drugs_Ruiz.csv",
        stitch_voc = f'{stitch_path}' + "stitch_chemical_drugbank_vocabulary.tsv",
        stitch_data = f'{stitch_path}' + "9606.actions.v5.0.tsv.gz",
        prot_gene_ensembl_voc = biomart_mapping_file
    output:
        "output/drug_targets/{trait}_drug_Ruiz_gene_stitch.tsv"
    script:
        "scripts/drug_targets/drug_x_gene_stitch.py"

rule drug_drugbank_x_gene_dgidb:
    input:
        drugs = "output/drug_targets/{trait}_drugs_drugbank.csv",
        dgidb_data = f'{dgidb_path}' + "dgidb_interactions_202101.tsv",
        gene_voc = f'{dgidb_path}' + "dgidb_genes_202101.tsv",
        drug_voc = f'{dgidb_path}' + "dgidb_drugs_202101.tsv",
        drug_voc_chembl_db = "data/Chembl_to_drugbank_dgidb_202202.txt"
    output:
        "output/drug_targets/{trait}_drug_drugbank_gene_dgidb.tsv"
    script:
        "scripts/drug_targets/drug_x_gene_dgidb.py"

rule drug_drugbank_x_gene_stitch:
    input:
        drugs = "output/drug_targets/{trait}_drugs_drugbank.csv", # generated by scripts/drug_targets/GWAS_x_drug_drugBank.py
        stitch_voc = f'{stitch_path}' + "stitch_chemical_drugbank_vocabulary.tsv",
        stitch_data = f'{stitch_path}' + "9606.actions.v5.0.tsv.gz",
        prot_gene_ensembl_voc = biomart_mapping_file
    output:
        "output/drug_targets/{trait}_drug_drugbank_gene_stitch.tsv"
    script:
        "scripts/drug_targets/drug_x_gene_stitch.py"

rule drug_Chembl_x_gene_Chembl:
    input:
        drugs = "output/drug_targets/{trait}_drugs_Chembl.csv",
        target_map = "data/Chembl_drug_targets_ChEMBLID_UniProtIDs_EnsemblIDs.csv", # generated by scripts/drug_targets/Chembl_uniprot_to_ensembl.py
        drug_target = f'{chembl_path}' "Chembl_drug_mechanisms_20220517.tsv"
    output:
        "output/drug_targets/{trait}_drug_Chembl_gene_Chembl.tsv"
    script: 
        "scripts/drug_targets/drug_Chembl_x_gene_Chembl.py"

rule drug_database_statistics:
    input:
        drug_Ruiz = expand("output/drug_targets/{trait}_drugs_Ruiz.csv", trait = traits),
        drug_drugbank = expand("output/drug_targets/{trait}_drugs_drugbank.csv", trait = traits),
        drug_chembl = expand("output/drug_targets/{trait}_drugs_Chembl.csv", trait = traits),
        drug_Ruiz_dgidb = expand("output/drug_targets/{trait}_drug_Ruiz_gene_dgidb.tsv", trait = traits),
        drug_Ruiz_stitch = expand("output/drug_targets/{trait}_drug_Ruiz_gene_stitch.tsv", trait = traits),
        drug_drugbank_dgidb = expand("output/drug_targets/{trait}_drug_drugbank_gene_dgidb.tsv", trait = traits),
        drug_drugbank_stitch = expand("output/drug_targets/{trait}_drug_drugbank_gene_stitch.tsv", trait = traits),
        drug_chembl_chembl = expand("output/drug_targets/{trait}_drug_Chembl_gene_Chembl.tsv", trait = traits)
    output:
        "output/drug_targets/drug_database_statistics.tsv"
    script:
        "scripts/drug_targets/drug_database_statistics.py"

rule drug_database_correlations:
    input:
        drug_Ruiz_dgidb = expand("output/drug_targets/{trait}_drug_Ruiz_gene_dgidb.tsv", trait = traits),
        drug_Ruiz_stitch = expand("output/drug_targets/{trait}_drug_Ruiz_gene_stitch.tsv", trait = traits),
        drug_drugbank_dgidb = expand("output/drug_targets/{trait}_drug_drugbank_gene_dgidb.tsv", trait = traits),
        drug_drugbank_stitch = expand("output/drug_targets/{trait}_drug_drugbank_gene_stitch.tsv", trait = traits),
        drug_chembl = expand("output/drug_targets/{trait}_drug_Chembl_gene_Chembl.tsv", trait = traits),
        gene_info = genome_annotation_file
    output:
        "output/drug_targets/drug_database_correlations.txt"
    script:
        "scripts/drug_targets/drug_database_correlations.R"

rule generate_random_scores:
    input:
        gene_info = genome_annotation_file
    params:
        seed = lambda wildcards: int("{}".format(str(seed_dic[wildcards.trait])))
    output:
        "output/method/random_gene_{trait}_scores.tsv"
    script:
        "scripts/direct_overlap/generate_random_scores.R"

rule direct_overlap_enrichment:
    input:
        drug_Ruiz_dgidb = "output/drug_targets/{trait}_drug_Ruiz_gene_dgidb.tsv",
        drug_Ruiz_stitch = "output/drug_targets/{trait}_drug_Ruiz_gene_stitch.tsv",
        drug_drugbank_dgidb = "output/drug_targets/{trait}_drug_drugbank_gene_dgidb.tsv",
        drug_drugbank_stitch = "output/drug_targets/{trait}_drug_drugbank_gene_stitch.tsv",
        drug_chembl = "output/drug_targets/{trait}_drug_Chembl_gene_Chembl.tsv",
        prot_gene_ensembl_voc = biomart_mapping_file,
        pascal = f'{Pascal_result_path}' + "Pascal_Ensembl_scoring_{trait}.tsv",
        pascal_ukbb = lambda wildcards: f'{Pascal_result_path}' + "Pascal_Ensembl_scoring_{}.tsv".format(str(traits_UKBB_dic[wildcards.trait])),
        mr = "output/omics_MR/{trait}_aggregated_gtex_eQTLGen_MR.tsv",
        mr_gtex = "output/omics_MR/{trait}_aggregated_gtex_MR.tsv",
        mr_UKBB = lambda wildcards: "output/omics_MR/{}_aggregated_gtex_eQTLGen_MR.tsv".format(str(traits_UKBB_dic[wildcards.trait])),
        gene_info = genome_annotation_file,
        mr_protein = "output/omics_MR/protein_deCODE_{trait}_MR.msmr",
        mr_protein_UKBB = lambda wildcards: "output/omics_MR/protein_deCODE_{}_MR.msmr".format(str(traits_UKBB_dic[wildcards.trait])),
        protein_info = f'{deCODE_data_path}' + "Decode_Somamer_info_TSS_GRCh38.csv", # mapping of protein IDs to EnsemblIds
        exome = lambda wildcards: "output/exome_genes/{}_gene_burden_M3_1_associations.tsv".format(str(traits_exome_dic[wildcards.trait])),
        random = "output/method/random_gene_{trait}_scores.tsv"
    output:
        temp("output/direct_overlap/direct_overlap_{trait}_enrichment.tsv"),
        temp("output/direct_overlap/direct_overlap_{trait}_target_info.tsv")
    script:
        "scripts/direct_overlap/direct_overlap_enrichment_withtargetinfo.R"

rule aggregate_direct_overlap_enrichment:
    input:
        files = expand("output/direct_overlap/direct_overlap_{trait}_enrichment.tsv", trait = traits)
    params:
        3
    output:
        "output/direct_overlap/direct_overlap_aggregated_enrichment.tsv"
    script:
        "scripts/utilities/aggregate_tables.py"

rule aggregate_direct_overlap_target_info:
    input:
        files = expand("output/direct_overlap/direct_overlap_{trait}_target_info.tsv", trait = traits)
    params:
        3
    output:
        "output/direct_overlap/direct_overlap_aggregated_target_info.tsv"
    script:
        "scripts/utilities/aggregate_tables.py"

rule direct_overlap_running_OR:
    input:
        drug_Ruiz_dgidb = "output/drug_targets/{trait}_drug_Ruiz_gene_dgidb.tsv",
        drug_Ruiz_stitch = "output/drug_targets/{trait}_drug_Ruiz_gene_stitch.tsv",
        drug_drugbank_dgidb = "output/drug_targets/{trait}_drug_drugbank_gene_dgidb.tsv",
        drug_drugbank_stitch = "output/drug_targets/{trait}_drug_drugbank_gene_stitch.tsv",
        drug_chembl = "output/drug_targets/{trait}_drug_Chembl_gene_Chembl.tsv",
        pascal = f'{Pascal_result_path}' + "Pascal_Ensembl_scoring_{trait}.tsv",
        pascal_ukbb = lambda wildcards: f'{Pascal_result_path}' + "Pascal_Ensembl_scoring_{}.tsv".format(str(traits_UKBB_dic[wildcards.trait])),
        mr = "output/omics_MR/{trait}_aggregated_gtex_eQTLGen_MR.tsv",
        mr_gtex = "output/omics_MR/{trait}_aggregated_gtex_MR.tsv",
        mr_UKBB = lambda wildcards: "output/omics_MR/{}_aggregated_gtex_eQTLGen_MR.tsv".format(str(traits_UKBB_dic[wildcards.trait])),
        gene_info = genome_annotation_file,
        mr_protein = "output/omics_MR/protein_deCODE_{trait}_MR.msmr",
        protein_info = f'{deCODE_data_path}' + "Decode_Somamer_info_TSS_GRCh38.csv", # mapping of protein IDs to EnsemblIds
        exome = lambda wildcards: "output/exome_genes/{}_gene_burden_M3_1_associations.tsv".format(str(traits_exome_dic[wildcards.trait]))
    output:
        temp("output/direct_overlap/direct_overlap_{trait}_running_OR.tsv")
    script:
        "scripts/direct_overlap/direct_overlap_running_OR.R"

rule aggregate_direct_overlap_running_OR:
    input:
        files = expand("output/direct_overlap/direct_overlap_{trait}_running_OR.tsv", trait = traits)
    params:
        3
    output:
        "output/direct_overlap/direct_overlap_aggregated_running_OR.tsv"
    script:
        "scripts/utilities/aggregate_tables.py"

rule direct_overlap_pairwise_drug_overlap:
    input:
        drug_drugbank_dgidb = expand("output/drug_targets/{trait}_drug_drugbank_gene_dgidb.tsv", trait = traits),
        pascal = expand(pascal = f'{Pascal_result_path}' + "Pascal_Ensembl_scoring_{trait}.tsv", trait = traits),
        mr = expand("output/omics_MR/{trait}_aggregated_gtex_eQTLGen_MR.tsv", trait = traits),
        gene_info = genome_annotation_file,
        mr_protein = expand("output/omics_MR/protein_deCODE_{trait}_MR.msmr", trait = traits),
        protein_info = f'{deCODE_data_path}' + "Decode_Somamer_info_TSS_GRCh38.csv", # mapping of protein IDs to EnsemblIds
        exome = ["output/exome_genes/" + traits_exome_dic[trait] + "_gene_burden_M3_1_associations.tsv" for trait in traits]
    output:
        "output/direct_overlap/direct_overlap_pairwise_drug_overlap_Jaccard.tsv"
    script:
        "scripts/direct_overlap/direct_overlap_pairwise_drug_overlap.py"

rule direct_overlap_pairwise_drug_overlap_ukbb:
    input:
        drug_drugbank_dgidb = expand("output/drug_targets/{trait}_drug_drugbank_gene_dgidb.tsv", trait = traits),
        pascal = [f'{Pascal_result_path}' + "Pascal_Ensembl_scoring_" + traits_UKBB_dic[trait] + ".tsv" for trait in traits],
        mr = ["output/omics_MR/" + traits_UKBB_dic[trait] + "_aggregated_gtex_eQTLGen_MR.tsv" for trait in traits],
        gene_info = genome_annotation_file,
        mr_protein = ["output/omics_MR/protein_deCODE_" + traits_UKBB_dic[trait] + "_MR.msmr" for trait in traits],
        protein_info = f'{deCODE_data_path}' + "Decode_Somamer_info_TSS_GRCh38.csv", # mapping of protein IDs to EnsemblIds
        exome = ["output/exome_genes/" + traits_exome_dic[trait] + "_gene_burden_M3_1_associations.tsv" for trait in traits]
    output:
        "output/direct_overlap/direct_overlap_pairwise_drug_overlap_Jaccard_ukbb.tsv"
    script:
        "scripts/direct_overlap/direct_overlap_pairwise_drug_overlap.py"

rule transcript_heritability_drug_targets:
    input:
        eqtlgen = f'{eQTLGen_data_path}' + "cis-eQTL.epi",
        mr_eqtlgen = "output/omics_MR/transcript_eQTLGen_{trait}_MR.msmr",
        her = f'{reml_path}' + "reml_result_protein_coding_eQTLGen.tsv",
        gene_info = genome_annotation_file,
        drug_Ruiz_dgidb = "output/drug_targets/{trait}_drug_Ruiz_gene_dgidb.tsv",
        drug_Ruiz_stitch = "output/drug_targets/{trait}_drug_Ruiz_gene_stitch.tsv",
        drug_drugbank_dgidb = "output/drug_targets/{trait}_drug_drugbank_gene_dgidb.tsv",
        drug_drugbank_stitch = "output/drug_targets/{trait}_drug_drugbank_gene_stitch.tsv",
        drug_chembl = "output/drug_targets/{trait}_drug_Chembl_gene_Chembl.tsv"
    output:
        temp("output/direct_overlap/transcript_remlheritability_{trait}_drug_targets.tsv")
    threads: 1
    resources:
        mem_mb = 10000
    script:
        "scripts/direct_overlap/transcript_remlheritability_drug_targets.R"

rule aggregate_transcript_heritability_drug_targets:
    input:
        files = expand("output/direct_overlap/transcript_remlheritability_{trait}_drug_targets.tsv", trait = traits)
    params:
        3
    output:
        "output/direct_overlap/transcript_remlheritability_aggregated_drug_targets.tsv"
    script:
        "scripts/utilities/aggregate_tables.py"

rule decode_enrichment_drug_target:
    input:
        mr = "output/omics_MR/heritability_decode_MR_ld01.msmr",
        protein_info = f'{deCODE_data_path}' + "Decode_Somamer_info_TSS_GRCh38.csv", # mapping of protein IDs to EnsemblIds
        decode = f'{deCODE_data_path}' + "cis-pQTL.epi",
        her = f'{reml_path}' + "converged_reml_results_data_decode.tsv",
        gene_info = genome_annotation_file,
        drug_Ruiz_dgidb = "output/drug_targets/{trait}_drug_Ruiz_gene_dgidb.tsv",
        drug_Ruiz_stitch = "output/drug_targets/{trait}_drug_Ruiz_gene_stitch.tsv",
        drug_drugbank_dgidb = "output/drug_targets/{trait}_drug_drugbank_gene_dgidb.tsv",
        drug_drugbank_stitch = "output/drug_targets/{trait}_drug_drugbank_gene_stitch.tsv",
        drug_chembl = "output/drug_targets/{trait}_drug_Chembl_gene_Chembl.tsv"
    output:
        temp("output/direct_overlap/protein_remlheritability_{trait}_drug_targets.tsv")
    threads: 1
    resources:
        mem_mb = 10000
    script:
        "scripts/direct_overlap/protein_remlheritability_drug_targets.R"

rule aggregate_protein_heritability_drug_targets:
    input:
        files = expand("output/direct_overlap/protein_remlheritability_{trait}_drug_targets.tsv", trait = traits)
    params:
        3
    output:
        "output/direct_overlap/protein_remlheritability_aggregated_drug_targets.tsv"
    script:
        "scripts/utilities/aggregate_tables.py"

rule prepare_network_STRING:
    input:
        network = f'{nw_path}' + "9606.protein.links.v11.5.txt",
        prot_gene_ensembl_voc = biomart_mapping_file
    output:
        "output/networks/network_STRING.tsv"
    script:
        "scripts/networks/prepare_network_STRING.py"

rule prepare_network_CoXRNAseqf:
    input:
        network = f'{nw_path}' + "gene_coregulation_165_eigenvectors_protein_coding.txt",
        gene_info = genome_annotation_file # all protein coding
    threads: 1
    resources:
        mem_mb = 30000
    output:
        "output/networks/network_CoXRNAseqf.tsv"
    script:
        "scripts/networks/prepare_network_CoXRNAseqf.py"

rule prepare_network_FAVA:
    input:
        network = f'{nw_path}' + "Fava_Network_SingleCells_Proteomics_CombinedScores_Over15percent.tsv.gz",
        prot_gene_ensembl_voc = biomart_mapping_file
    output:
         "output/networks/network_FAVA.tsv"
    script:
        "scripts/networks/prepare_network_FAVA.py"

rule network_statistics:
    input:
        nw_files = expand("output/networks/network_{nw}.tsv", nw = networks)
    output:
        "output/networks/network_statistics.tsv"
    script:
        "scripts/networks/network_statistics.R"

rule network_diffusion_overlap_enrichment:
    input:
        drug_Ruiz_dgidb = "output/drug_targets/{trait}_drug_Ruiz_gene_dgidb.tsv",
        drug_Ruiz_stitch = "output/drug_targets/{trait}_drug_Ruiz_gene_stitch.tsv",
        drug_drugbank_dgidb = "output/drug_targets/{trait}_drug_drugbank_gene_dgidb.tsv",
        drug_drugbank_stitch = "output/drug_targets/{trait}_drug_drugbank_gene_stitch.tsv",
        drug_chembl = "output/drug_targets/{trait}_drug_Chembl_gene_Chembl.tsv",
        pascal = f'{Pascal_result_path}' + "Pascal_Ensembl_scoring_{trait}.tsv",
        pascal_ukbb = lambda wildcards: f'{Pascal_result_path}' + "Pascal_Ensembl_scoring_{}.tsv".format(str(traits_UKBB_dic[wildcards.trait])),
        mr = "output/omics_MR/{trait}_aggregated_gtex_eQTLGen_MR.tsv",
        mr_UKBB = lambda wildcards: "output/omics_MR/{}_aggregated_gtex_eQTLGen_MR.tsv".format(str(traits_UKBB_dic[wildcards.trait])),
        gene_info = genome_annotation_file,
        mr_protein = "output/omics_MR/protein_deCODE_{trait}_MR.msmr",
        protein_info = f'{deCODE_data_path}' + "Decode_Somamer_info_TSS_GRCh38.csv", # mapping of protein IDs to EnsemblIds
        exome = lambda wildcards: "output/exome_genes/{}_gene_burden_M3_1_associations.tsv".format(str(traits_exome_dic[wildcards.trait])),
        prot_gene_ensembl_voc = biomart_mapping_file,
        nw = "output/networks/network_{nw}.tsv"
    params:
        seed = lambda wildcards: int("{}".format(str(seed_dic[wildcards.trait])))
    output:
        diffusion = "output/networks/diffusion_profiles/diffusion_{trait}_mainmethods.{nw}.r{r}.w{weighted}.t{thresh}.tsv",
        overlap = temp("output/networks/diffusion_overlap_{trait}_enrichment.{nw}.r{r}.w{weighted}.t{thresh}.tsv"),
        target_info = temp("output/networks/diffusion_overlap_{trait}_target_info.{nw}.r{r}.w{weighted}.t{thresh}.tsv")
    script:
        "scripts/networks/network_diffusion.R"

rule aggregate_network_diffusion_overlap_enrichment:
    input:
        files = expand("output/networks/diffusion_overlap_{trait}_enrichment.{{nw}}.r{{r}}.w{{weighted}}.t{{thresh}}.tsv", trait = traits)
    params:
        2
    output:
        "output/networks/diffusion_overlap_aggregated_enrichment.{nw}.r{r}.w{weighted}.t{thresh}.tsv"
    script:
        "scripts/utilities/aggregate_tables.py"

rule aggregate_network_diffusion_overlap_target_info:
    input:
        files = expand("output/networks/diffusion_overlap_{trait}_target_info.{{nw}}.r{{r}}.w{{weighted}}.t{{thresh}}.tsv", trait = traits)
    params:
        2
    output:
        "output/networks/diffusion_overlap_aggregated_target_info.{nw}.r{r}.w{weighted}.t{thresh}.tsv"
    script:
        "scripts/utilities/aggregate_tables.py"

rule network_degree_overlap_enrichment:
    input:
        drug_Ruiz_dgidb = "output/drug_targets/{trait}_drug_Ruiz_gene_dgidb.tsv",
        drug_Ruiz_stitch = "output/drug_targets/{trait}_drug_Ruiz_gene_stitch.tsv",
        drug_drugbank_dgidb = "output/drug_targets/{trait}_drug_drugbank_gene_dgidb.tsv",
        drug_drugbank_stitch = "output/drug_targets/{trait}_drug_drugbank_gene_stitch.tsv",
        drug_chembl = "output/drug_targets/{trait}_drug_Chembl_gene_Chembl.tsv",
        nw = "output/networks/network_{nw}.tsv"
    output:
        diffusion = "output/networks/diffusion_profiles/diffusion_{trait}_degree.{nw}.w{weighted}.t{thresh}.tsv",
        overlap = temp("output/networks/diffusion_overlap_{trait}_degree_enrichment.{nw}.w{weighted}.t{thresh}.tsv")
    script:
        "scripts/networks/network_degree.R"

rule aggregate_degree_overlap_enrichment:
    input:
        files = expand("output/networks/diffusion_overlap_{trait}_degree_enrichment.{{nw}}.w{{weighted}}.t{{thresh}}.tsv", trait = traits)
    params:
        2
    output:
        "output/networks/diffusion_overlap_aggregated_degree_enrichment.{nw}.w{weighted}.t{thresh}.tsv"
    script:
        "scripts/utilities/aggregate_tables.py"

