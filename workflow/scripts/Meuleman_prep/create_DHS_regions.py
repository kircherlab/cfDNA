import pandas as pd

input_file = snakemake.input["DHS_summit"]
input_header = snakemake.input["header"]
output_summitAnno_all = snakemake.output["summitAnno_all"]
output_summitAnno_pattern = snakemake.params["summitAnno_pattern"]
output_DHS_all = snakemake.output["DHS_all"]
output_DHS_pattern = snakemake.params["DHS_pattern"]

header = pd.read_csv(input_header).columns.to_list()
df = pd.read_csv(input_file, sep="\t", names=header)

# filter for top quartile 
df_filtered = df.loc[df["mean_signal"] > float(7.617222e-01) ]

# filter df for Top1000
Top1000_DHS_summitAnno_all = pd.DataFrame()
for key in df_filtered["component"].unique():
    name = key.replace(" ","_").replace("/_","").replace(".","").lower()
    outfile = output_summitAnno_pattern.format(name)
    component = df_filtered.loc[df_filtered["component"] == key].sort_values(by="mean_signal", ascending=False).iloc[:1000]
    component["chrom"] = component["chrom"].str.replace("chr","")
    region = pd.DataFrame()
    region["ID"] = component["name"]
    region["chrom"] = component["chrom"]
    region["start"] = component["start"]
    region["end"] = component["end"]
    region["strand"] = component["strand"]
    Top1000_DHS_summitAnno_all = pd.concat([Top1000_DHS_summitAnno_all,region]) # append the component dfs to the Top1000_DHS_summitAnno_all df
    region.to_csv(outfile, sep="\t", header=None, index=None)

# save 
Top1000_DHS_summitAnno_all.to_csv(output_summitAnno_all, sep="\t", header=None, index=None)

# filter df for Top1000
Top1000_DHS_all = pd.DataFrame()
for key in df_filtered["component"].unique():
    name = key.replace(" ","_").replace("/_","").replace(".","").lower()
    outfile = output_DHS_pattern.format(name)
    component = df_filtered.loc[df_filtered["component"] == key].sort_values(by="mean_signal", ascending=False).iloc[:1000]
    #component["chrom"] = component["chrom"].str.replace("chr","")
    region = pd.DataFrame()
    region["chrom"] = component["chrom"]
    region["start"] = component["start"]
    region["end"] = component["end"]
    region["ID"] = component["name"]
    region["score"] = component["mean_signal"]
    region["strand"] = component["strand"]
    Top1000_DHS_all = pd.concat([Top1000_DHS_all,region]) # append the component dfs to the Top1000_DHS_summitAnno_all df
    region.to_csv(outfile, sep="\t", header=None, index=None)

# save 
Top1000_DHS_summitAnno_all.to_csv(output_DHS_all, sep="\t", header=None, index=None)