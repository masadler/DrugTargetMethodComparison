import pandas as pd

"""

STITCH:

Target selection:
1) Score >= 700
2) Only inhibition and activation (i.e. "action" is not null and a_is_acting == t)

Gene ID:
- Protein ensembl identifiers were mapped to gene ensembl identifiers using biomart (GRCh37)

Drug ID:
- Chemical ID (from Pubmed) to DrugBank ID mapping is provided by STITCH (chemical.sources.v5.0.tsv.gz).

"""

drug_df = pd.read_csv(snakemake.input["drugs"])

# Map drugs to STITCH identifiers
stitch_df = pd.read_csv(snakemake.input["stitch_voc"], sep = '\t', names = ["CID1", "CID2", "source", "drugbank_id"])
stitch_df["CID"] = stitch_df["CID1"].apply(lambda x: "CID" + x[4:12])
stitch_df = stitch_df[["drugbank_id", "CID"]]
stitch_df = stitch_df.drop_duplicates()
drug_stitch_df = drug_df.merge(stitch_df, on = "drugbank_id")

# get drug targets from STITCH
stitch_df = pd.read_csv(snakemake.input['stitch_data'], sep = '\t', compression = "gzip")
stitch_df = stitch_df[stitch_df.score >= 700]
stitch_df["item_id_a"] = stitch_df["item_id_a"].apply(lambda x: "CID" + x[4:12]) # does not work for proteins, but we are only interested in drugs for item a
stitch_df = stitch_df.merge(drug_stitch_df, left_on = "item_id_a", right_on = "CID")

# only activating & inhibiting
stitch_df = stitch_df[(~stitch_df.action.isnull()) & (stitch_df.a_is_acting == 't')]

# merge Ensembl Identifier & name
gene_df = pd.read_csv(snakemake.input["prot_gene_ensembl_voc"], sep = '\t')
stitch_df["ensembl_peptide_id"] = stitch_df["item_id_b"].apply(lambda x: x[5:])
stitch_df = stitch_df.merge(gene_df, on = "ensembl_peptide_id")
stitch_df = stitch_df.rename(columns = {'ensembl_peptide_id': 'Ensembl_peptide', 'ensembl_gene_id': 'EnsemblId', 'hgnc_symbol': 'gene_name'})

stitch_df = stitch_df[['drugbank_id', 'db_name', 'EnsemblId', 'gene_name', 'Ensembl_peptide']]
stitch_df = stitch_df.drop_duplicates()
stitch_df["weight"] = 1

stitch_df.to_csv(snakemake.output[0], index = False)