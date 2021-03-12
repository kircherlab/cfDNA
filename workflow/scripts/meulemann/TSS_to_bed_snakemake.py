import pandas as pd

input = snakemake.input[0]
output = snakemake.output[0]

name = input.split("_Meulemann")[0].split("/")[-1]


TSSdf = pd.read_csv(input, sep=":", header=None)
regions = pd.DataFrame()
regions = regions.assign(chrom =  "chr" + TSSdf[0].astype(str))
regions["start"] = TSSdf[1]
regions["end"] = TSSdf[1]+1
regions["name"] = name
regions["score"] = 0
regions["strand"] = TSSdf[2]
regions.to_csv(output, sep="\t", header=None, index=None)