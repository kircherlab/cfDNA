import pandas as pd

input_file = snakemake.input["DHS_summit"]
input_header = snakemake.input["header"]
output_file_continuous = snakemake.output["DHS_continuous"]
output_file_binary = snakemake.output["DHS_binary"]


header = pd.read_csv(input_header).columns.to_list()
df = pd.read_csv(input_file, sep="\t", names=header)
# filter for top quartile 
df_filtered = df.loc[df["mean_signal"] > float(7.617222e-01) ]

# create continuous matrix with ID (name) as index and mean signal as value
# contains 16 columns, one for each component
DHS_matrix = pd.DataFrame()
for key in df_filtered["component"].unique():
    name = key.replace(" ","_").replace("/_","").replace(".","").lower()
    component = df_filtered.loc[df_filtered["component"] == key].sort_values(by="mean_signal", ascending=False).iloc[:1000] #.set_index("identifier")
    component[name] = component["mean_signal"]
    component = component.set_index("name")
    DHS_matrix = pd.concat([DHS_matrix, component[name]], axis=1).fillna(0)

DHS_matrix.to_csv(output_file_continuous, sep="\t")
# binarize continuos matrix
DHS_matrix.mask(DHS_matrix > 0, 1).to_csv(output_file_binary, sep="\t")