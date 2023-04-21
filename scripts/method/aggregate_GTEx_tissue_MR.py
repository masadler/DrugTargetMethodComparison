import numpy as np
import pandas as pd

from glob import glob
from pathlib import Path

trait = snakemake.wildcards["trait"]
dir_path = Path(snakemake.params["path"]).resolve()
mr_files = list(dir_path.glob('transcript_*_tissue_' + trait + '_MR.msmr'))

df = pd.read_csv(snakemake.input[0], sep = '\t')
df = df[["ProbeID", "ProbeName", "b_ivw", "p_ivw", "h2cis", "n_iv", "topSNP"]]
df = df.rename(columns = {"p_ivw": "p_eqtlgen", "b_ivw": "b_eqtlgen", "h2cis": "h2cis_eqtlgen", "n_iv": "n_iv_eqtlgen", "topSNP": "topSNP_eqtlgen"})
df['pmin'] = df.p_eqtlgen

for f in mr_files:
    tdf = pd.read_csv(f, sep = '\t')
    tissue = str(f).split('transcript_')[1].split('_tissue')[0]
    tdf = tdf[["ProbeID", "b_ivw", "p_ivw", "h2cis", "n_iv", "topSNP"]]
    tdf['ProbeID'] = tdf['ProbeID'].apply(lambda x: x.split('.')[0])
    df = df.merge(tdf, on = "ProbeID", how = "outer")
    df["pmin"] = df.apply(lambda x: np.nanmin(np.array([x["pmin"], x["p_ivw"]])), axis = 1)
    df = df.rename(columns = {"b_ivw": "b_" + tissue, "p_ivw": "p_" + tissue, "h2cis": "h2cis_" + tissue, "n_iv": "n_iv_" + tissue, "topSNP": "topSNP_" + tissue})

df.to_csv(snakemake.output[0], index = False, sep = '\t')

## get dataset with only the tissue with minimum p-value

# get tissue for each gene
def get_col_name(row):    
    b = (df.ix[row.name] == row['pmin'])
    return b.idxmax()

df_t = df[["ProbeID", "pmin"]]
df = df.drop('pmin', axis = 1)
df_t['p_tissue'] = df_t.apply(get_col_name, axis=1)
df_t['Tissue'] = df_t.apply(lambda x: x['p_tissue'].split("p_")[1], axis = 1)

df = df.merge(df_t[["ProbeID", "pmin", "Tissue"]], on = "ProbeID")
df = df.rename(columns = {"pmin": "p_tissue"})

df['b_tissue'] = df.apply(lambda x: x["b_" + x["Tissue"]], axis = 1)
df['h2cis_tissue'] = df.apply(lambda x: x["h2cis_" + x["Tissue"]], axis = 1)
df['n_iv_tissue'] = df.apply(lambda x: x["n_iv_" + x["Tissue"]], axis = 1)
df['topSNP_tissue'] = df.apply(lambda x: x["topSNP_" + x["Tissue"]], axis = 1)

df[["ProbeID", "Tissue", "b_tissue", "p_tissue", "h2cis_tissue", "n_iv_tissue", "topSNP_tissue"]].to_csv(snakemake.output[1], index = False, sep = '\t')

## use only GTEx mr files

f = mr_files[0]
df = pd.read_csv(f, sep = '\t')
tissue = str(f).split('transcript_')[1].split('_tissue')[0]
df = df[["ProbeID", "ProbeName", "b_ivw", "p_ivw", "h2cis"]]
df = df.rename(columns = {"b_ivw": "b_" + tissue, "p_ivw": "p_" + tissue, "h2cis": "h2cis_" + tissue})
df['pmin'] = df["p_" + tissue]

for f in mr_files[1:]:
    tdf = pd.read_csv(f, sep = '\t')
    tissue = str(f).split('transcript_')[1].split('_tissue')[0]
    tdf = tdf[["ProbeID", "b_ivw", "p_ivw", "h2cis"]]
    tdf['ProbeID'] = tdf['ProbeID'].apply(lambda x: x.split('.')[0])
    df = df.merge(tdf, on = "ProbeID", how = "outer")
    df["pmin"] = df.apply(lambda x: np.nanmin(np.array([x["pmin"], x["p_ivw"]])), axis = 1)
    df = df.rename(columns = {"b_ivw": "b_" + tissue, "p_ivw": "p_" + tissue, "h2cis": "h2cis_" + tissue})

df.to_csv(snakemake.output[2], index = False, sep = '\t')
