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
WPS_back = snakemake.input["WPS_back"]
WPS_back_refs = snakemake.input["WPS_back_ref"]
COV_back = snakemake.input["COV_back"]
COV_back_refs = snakemake.input["COV_back_ref"]

sample_ID = snakemake.params["sample"]
ref_IDs = snakemake.params["ref_IDs"]
target = snakemake.params["target"]
outfile = snakemake.output[0]
overlay_mode = snakemake.params["overlay_mode"]

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


def add_sample(path_a: str, path_b: str, overlay_mode:str = "mean",):
    """Reads .csv file, calculates mean over all rows and divides by trimmed mean.

    Args:
        path_a (str): [path to targets]
        path_b (str): [path to background]
        overlay_mode (str): [Mode of operation, e.g. mean, median]

    Returns:
        [Pandas Object]: [Pandas object containing normalized/processed data]
    """

    if overlay_mode == "mean":
        sample_a = pd.read_csv(path_a, header=None).mean()
    elif overlay_mode == "median":
        sample_a = pd.read_csv(path_a, header=None).median()
    else:
        raise ValueError(f"{overlay_mode} is not a valid keyword.")
    sample_b = pd.read_csv(path_b, header=None).mean(axis=1)
    sample = sample_a / stats.trim_mean(sample_b, 0.1)
    return sample


# load tables containing position specific scores for all defined target regions
# average over all regions per sample and substract the trimmed mean to normalise


av_WPS = pd.DataFrame()
av_WPS[sample_ID] = add_sample(WPS, WPS_back)
for (ref_ID, WPS_ref, WPS_back_ref) in zip(ref_IDs, WPS_refs, WPS_back_refs):
    av_WPS[ref_ID] = add_sample(WPS_ref, WPS_back_ref)

av_WPS["position"] = calculate_flanking_regions(len(av_WPS))
av_WPS = av_WPS.set_index("position")

av_COV = pd.DataFrame()
av_COV[sample_ID] = add_sample(COV, COV_back)
for (ref_ID, COV_ref, COV_back_ref) in zip(ref_IDs, COV_refs, COV_back_refs):
    av_COV[ref_ID] = add_sample(COV_ref, COV_ref)

av_COV["position"] = calculate_flanking_regions(len(av_COV))
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
