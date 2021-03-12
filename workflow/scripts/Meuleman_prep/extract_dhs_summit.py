import pandas as pd
import csv

input_file = snakemake.input["DHS_Index"]
output_file = snakemake.output["DHS_summit"]

df = pd.read_csv(input_file, sep="\t")

df["identifier"] = "chr" + df["identifier"].astype(str).str.replace(".","-")

df_summit = pd.DataFrame()
df_summit["chrom"] = df["seqname"]
df_summit["start"] = df["summit"]
df_summit["end"] = df["summit"] +1
df_summit["name"] = df["identifier"]
df_summit["mean_signal"] = df["mean_signal"]
df_summit["strand"] = "."
df_summit["component"] = df["component"]

df_summit.to_csv(output_file, sep="\t", header=None, index=None)

with open(output_file + ".header", "w") as file:
    writer=csv.writer(file)
    writer.writerow(df_summit.columns.to_list())