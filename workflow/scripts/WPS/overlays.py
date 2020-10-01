#!/usr/bin/env python

"""
:Author: Sebastian Röner
:Contact: sebastian.roener@charite.de
:Date: 08.09.2020
"""

import matplotlib
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats

matplotlib.use("pdf")


# get variables from snakefile

WPS = snakemake.input["WPS"]
WPS_refs = snakemake.input["WPS_ref"]
COV = snakemake.input["COV"]
COV_refs = snakemake.input["COV_ref"]
sample_ID = snakemake.params["sample"]
ref_IDs = snakemake.params["ref_IDs"]
target = snakemake.params["target"]
outfile = snakemake.output[0]

# def functions


def calculate_flanking_regions(val: int):
    """Calculates flanking regions for point of interest.

    Args:
        val (int): should be length of value vector

    Raises:
        TypeError: Only integers are allowed

    Returns:
        [iterator]: range of values around center point (e.g. range(-1000,1000))
    """

    if not isinstance(val, int):
        raise TypeError("Only integers are allowed")

    if val % 2 == 0:
        flank = int(val / 2)
        region = range(-flank, flank)
    elif val % 2 == 1:
        flank_l = int(val / 2 - 0.5)
        flank_r = int(val / 2 + 0.5)
        region = range(-flank_l, flank_r)
    return region


def add_sample(path: str):
    """Reads .csv file, calculates mean over all rows and divides by trimmed mean.

    Args:
        path (str): Path to a .csv file

    Returns:
        [type]: [description]
    """
    sample = pd.read_csv(path, header=None).mean()
    sample = sample / stats.trim_mean(sample, 0.25)
    return sample


# load tables containing position specific scores for all defined target regions
# average over all regions per sample and substract the trimmed mean to normalise


av_WPS = pd.DataFrame()
av_WPS[sample_ID] = add_sample(WPS)
for (ref_ID, WPS_ref) in zip(ref_IDs, WPS_refs):
    av_WPS[ref_ID] = add_sample(WPS_ref)

av_WPS["position"] = calculate_flanking_regions(len(av_WPS))
av_WPS = av_WPS.set_index("position")

av_COV = pd.DataFrame()
av_COV[sample_ID] = add_sample(COV)
for (ref_ID, COV_ref) in zip(ref_IDs, COV_refs):
    av_COV[ref_ID] = add_sample(COV_ref)

av_COV["position"] = calculate_flanking_regions(len(av_WPS))
av_COV = av_COV.set_index("position")

# create line plots and save to a single pdf

with PdfPages(outfile) as pdf:
    Fig_WPS = av_WPS.plot(
        title=f"adjusted WPS: {target} target regions",
        xlabel="Position relative to target site",
        ylabel="normalized WPS",
    )
    Fig_Cov = av_COV.plot(
        title=f"adjusted read coverage: {target} target regions",
        xlabel="Position relative to target site",
        ylabel="normalized read coverage",
    )
    pdf.savefig(Fig_WPS.get_figure())
    pdf.savefig(Fig_Cov.get_figure())
