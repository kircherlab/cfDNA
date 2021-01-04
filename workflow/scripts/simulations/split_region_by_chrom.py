#!/usr/bin/env python

"""
:Author: Sebastian Röner
:Contact: sebastian.roener@charite.de
:Date: 04.01.2021
"""

import pandas as pd

bedfile = snakemake.input["Region_file"]
chrom = snakemake.params["chrom"]
outfile = snakemake.output[0]

print(bedfile)
print(outfile)
print(chrom)


def split_by_chrom(bedfile, chrom, outfile):
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
    bed_df.loc[bed_df["chrom"] == chrom].to_csv(
        outfile, sep="\t", header=None, index=None
    )


if __name__ == "__main__":
    split_by_chrom(bedfile, chrom, outfile)
