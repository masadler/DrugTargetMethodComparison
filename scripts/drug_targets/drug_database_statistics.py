import pandas as pd

statistics = ["trait", "drugs_chembl", "drugs_drugbank", "drugs_ruiz", "targets_chembl", "targets_drugbank_dgidb", 
              "targets_drugbank_stitch", "targets_ruiz_dgidb", "targets_ruiz_stitch"]
result = {key : [] for key in statistics}

drug_chembl_files = snakemake.input["drug_chembl"]
drug_drugbank_files = snakemake.input["drug_drugbank"]
drug_ruiz_files = snakemake.input["drug_Ruiz"]

target_chembl = snakemake.input["drug_chembl_chembl"]
target_drugbank_dgidb = snakemake.input["drug_drugbank_dgidb"]
target_drugbank_stitch = snakemake.input["drug_drugbank_stitch"]
target_ruiz_dgidb = snakemake.input["drug_Ruiz_dgidb"]
target_ruiz_stitch = snakemake.input["drug_Ruiz_stitch"]

traits = [f.split("output/drug_targets/")[1].split("_")[0] for f in drug_chembl_files]

for i, trait in enumerate(traits):

    result["trait"].append(trait)

    chembl_df = pd.read_csv(drug_chembl_files[i])
    drugbank_df = pd.read_csv(drug_drugbank_files[i])
    ruiz_df = pd.read_csv(drug_ruiz_files[i])

    result["drugs_chembl"].append(chembl_df.Chembl_id.nunique())
    result["drugs_drugbank"].append(drugbank_df.drugbank_id.nunique())
    result["drugs_ruiz"].append(ruiz_df.drugbank_id.nunique())

    chembl_df = pd.read_csv(target_chembl[i])
    drugbank_dgidb_df = pd.read_csv(target_drugbank_dgidb[i])
    drugbank_stitch_df = pd.read_csv(target_drugbank_stitch[i])
    ruiz_dgidb_df = pd.read_csv(target_ruiz_dgidb[i])
    ruiz_stitch_df = pd.read_csv(target_ruiz_stitch[i])

    result["targets_chembl"].append(chembl_df.EnsemblId.nunique())
    result["targets_drugbank_dgidb"].append(drugbank_dgidb_df.EnsemblId.nunique())
    result["targets_drugbank_stitch"].append(drugbank_stitch_df.EnsemblId.nunique())
    result["targets_ruiz_dgidb"].append(ruiz_dgidb_df.EnsemblId.nunique())
    result["targets_ruiz_stitch"].append(ruiz_stitch_df.EnsemblId.nunique())

result_df = pd.DataFrame(result)
result_df.to_csv(snakemake.output[0], index = False, sep = '\t')

