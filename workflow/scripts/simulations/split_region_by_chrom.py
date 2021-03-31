#!/usr/bin/env python

"""
:Author: Sebastian Röner
:Contact: sebastian.roener@charite.de
:Date: 04.01.2021
"""

import pandas as pd

bedfile = snakemake.input["Region_file"]
outfile = snakemake.output[0]
ref_chrom = snakemake.params["ref_chrom"]

print(ref_chrom)

#chroms_raw = bed_df["chrom"].unique().tolist()
#    chroms_ref = [f"chr{i}" for i in config["simulations"]["chroms"]]
#    chroms = [chrom for chrom in chroms_raw if chrom in chroms_ref]


def split_by_chrom(bedfile, ref_chrom, outfile):
    BEDCOLS = [
        "chrom",
        "chromStart",
        "chromEnd",
        "name",
        "score",
        "strand",
        "thickStart",
        "thickEnd",
        "itemRGB",
        "blockCount",
        "blockSizes",
        "blockStarts",
    ]
    bed_df = pd.read_csv(bedfile, sep="\t", header=None)
    colnames = []
    for i in range(len(bed_df.columns)):
        colnames.append(BEDCOLS[i])
    bed_df.columns = colnames
    chroms_raw = bed_df["chrom"].unique().tolist()
    chroms = [chrom for chrom in chroms_raw if chrom in ref_chrom]
    for chrom in chroms:
        c_outfile = outfile.replace("__snakemake_dynamic__", str(chrom))
        bed_df.loc[bed_df["chrom"] == chrom].to_csv(
        c_outfile, sep="\t", header=None, index=None
        )


if __name__ == "__main__":
    split_by_chrom(bedfile, ref_chrom, outfile)

