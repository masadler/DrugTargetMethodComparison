import pandas as pd

result_df = pd.DataFrame()
num_split = snakemake.params[0]

for file in snakemake.input["files"]:
    df = pd.read_csv(file, sep = '\t')
    trait = file.split("_")[num_split]
    df['trait'] = trait
    result_df = pd.concat([result_df, df])

result_df.to_csv(snakemake.output[0], sep = '\t', index = False)