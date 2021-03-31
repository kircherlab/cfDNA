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
import numpy as np

matplotlib.use("pdf")


# get variables from snakefile

WPS = snakemake.input["WPS"]
WPS_refs = snakemake.input["WPS_ref"]
WPS_back = snakemake.input["WPS_back"]
WPS_back_refs = snakemake.input["WPS_back_ref"]
WPS_sim = snakemake.input["WPS_sim"]
WPS_sim_refs = snakemake.input["WPS_sim_ref"]
WPS_sim_back = snakemake.input["WPS_sim_back"]
WPS_sim_back_refs = snakemake.input["WPS_sim_back_ref"]
COV = snakemake.input["COV"]
COV_refs = snakemake.input["COV_ref"]
COV_back = snakemake.input["COV_back"]
COV_back_refs = snakemake.input["COV_back_ref"]
COV_sim = snakemake.input["COV_sim"]
COV_sim_refs = snakemake.input["COV_sim_ref"]
COV_sim_back = snakemake.input["COV_sim_back"]
COV_sim_back_refs = snakemake.input["COV_sim_back_ref"]

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


def add_sample(target_path_a: str, target_path_b: str, sim_path_a: str, sim_path_b: str):
    """Reads .csv file, calculates mean over all rows and divides by trimmed mean.

    Args:
        path (str): Path to a .csv file

    Returns:
        [type]: [description]
    """
    sam_a = pd.read_csv(target_path_a, header=None)
    sam_b = pd.read_csv(target_path_b, header=None)
    sim_a = pd.read_csv(sim_path_a, header=None)
    sim_b = pd.read_csv(sim_path_b, header=None)
    sam_a_noo = sam_a[(np.abs(stats.zscore(sam_a)) < 3).all(axis=1)]
    sam_b_noo = sam_b[(np.abs(stats.zscore(sam_b)) < 3).all(axis=1)]
    sim_a_noo = sim_a[(np.abs(stats.zscore(sim_a)) < 3).all(axis=1)]
    sim_b_noo = sim_b[(np.abs(stats.zscore(sim_b)) < 3).all(axis=1)]
    sam = sam_a_noo.mean() / stats.trim_mean(sam_b_noo, 0.1)
    sim = sim_a_noo.mean() / stats.trim_mean(sim_b_noo, 0.1)
    sample = sam - sim
    return sample


# load tables containing position specific scores for all defined target regions
# average over all regions per sample and substract the trimmed mean to normalise


av_WPS = pd.DataFrame()
av_WPS[sample_ID] = add_sample(WPS, WPS_back, WPS_sim, WPS_sim_back)
for (ref_ID, WPS_ref, WPS_back_ref, WPS_sim_ref, WPS_sim_back_ref) in zip(ref_IDs, WPS_refs, WPS_back_refs, WPS_sim_refs, WPS_sim_back_refs):
    av_WPS[ref_ID] = add_sample(WPS_ref, WPS_back_ref, WPS_sim_ref, WPS_sim_back_ref)

av_WPS["position"] = calculate_flanking_regions(len(av_WPS))
av_WPS = av_WPS.set_index("position")

av_COV = pd.DataFrame()
av_COV[sample_ID] = add_sample(COV, COV_back, COV_sim, COV_sim_back)
for (ref_ID, COV_ref, COV_back_ref, COV_sim_ref, COV_sim_back_ref) in zip(ref_IDs, COV_refs, COV_back_refs, COV_sim_refs, COV_sim_back_refs):
    av_COV[ref_ID] = add_sample(COV_ref, COV_ref, COV_sim_ref, COV_sim_back_ref)

av_COV["position"] = calculate_flanking_regions(len(av_COV))
av_COV = av_COV.set_index("position")

# create line plots and save to a single pdf

with PdfPages(outfile) as pdf:
    Fig_WPS = av_WPS.iloc[300:-300].plot(
        title=f"adjusted WPS: {target} target regions",
        xlabel="Position relative to target site",
        ylabel="normalized WPS",
    )
    Fig_Cov = av_COV.iloc[300:-300].plot(
        title=f"adjusted read coverage: {target} target regions",
        xlabel="Position relative to target site",
        ylabel="normalized read coverage",
    )
    pdf.savefig(Fig_WPS.get_figure())
    pdf.savefig(Fig_Cov.get_figure())
