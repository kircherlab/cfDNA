#!/usr/bin/env python

"""
:Author: Sebastian Röner
:Contact: sebastian.roener@charite.de
:Date: 08.09.2020
"""

# imports

import pandas as pd

from scipy import stats
import matplotlib
matplotlib.use('pdf')
from matplotlib.backends.backend_pdf import PdfPages

# get variables from snakefile

WPS = snakemake.input["WPS"]
WPS_refs = snakemake.input["WPS_ref"]
COV = snakemake.input["COV"]
COV_refs = snakemake.input["COV_ref"]
sample_ID = snakemake.params["sample"]
ref_IDs = snakemake.params["ref_IDs"]
target = snakemake.params["target"]
outfile = snakemake.output[0]

# load tables containing position specific scores for all defined target regions +-1000 bp
# average over all regions per sample and substract the trimmed mean to normalise

av_WPS = pd.DataFrame()
av_WPS[sample_ID] = pd.read_csv(WPS, header=None).mean()
av_WPS[sample_ID] = av_WPS[sample_ID]-stats.trim_mean(av_WPS[sample_ID], 0.25)
for (ref_ID, WPS_ref) in zip(ref_IDs, WPS_refs):
    av_WPS[ref_ID] = pd.read_csv(WPS_ref, header=None).mean()
    av_WPS[ref_ID] = av_WPS[ref_ID]-stats.trim_mean(av_WPS[ref_ID], 0.25)

av_WPS["position"] = range(-1001,1001)
av_WPS = av_WPS.set_index("position")


av_COV = pd.DataFrame()
av_COV[sample_ID] = pd.read_csv(COV, header=None).mean()
av_COV[sample_ID] = av_COV[sample_ID]-stats.trim_mean(av_COV[sample_ID], 0.25)
for (ref_ID, COV_ref) in zip(ref_IDs, COV_refs):
    print(ref_ID, COV_ref)
    av_COV[ref_ID] = pd.read_csv(COV_ref, header=None).mean()
    av_COV[ref_ID] = av_COV[ref_ID]-stats.trim_mean(av_COV[ref_ID], 0.25)

av_COV["position"] = range(-1001,1001)
av_COV = av_COV.set_index("position")

# create line plots and save to a single pdf

with PdfPages(outfile) as pdf:
    Fig_WPS = av_WPS.plot(title=f"adjusted WPS: {target} target regions", xlabel="Position relative to target site", ylabel="adjusted WPS")
    Fig_Cov = av_COV.plot(title=f"adjusted read coverage: {target} target regions", xlabel="Position relative to target site", ylabel="adjusted read coverage")
    pdf.savefig(Fig_WPS.get_figure())
    pdf.savefig(Fig_Cov.get_figure())
