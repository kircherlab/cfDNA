#!/usr/bin/env python

"""
:Author: Sebastian Röner
:Contact: sebastian.roener@charite.de
:Date: 08.09.2020
"""

import matplotlib
import seaborn as sns
from matplotlib import pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats
from scipy.signal import savgol_filter

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
smoothing = snakemake.params["smoothing"]
rolling = snakemake.params["rolling"]
background_norm = snakemake.params["background_norm"]
edge_norm = snakemake.params["edge_norm"]
flank_edge = 500

win_len = snakemake.params["win_len"]
poly = snakemake.params["poly"]
compute_midpoint_coverage = snakemake.params["compute_midpoint_coverage"]
if compute_midpoint_coverage == True:
    MP = snakemake.input["MP"]
    MP_refs = snakemake.input["MP_ref"]
    MP_back = snakemake.input["MP_back"]
    MP_back_refs = snakemake.input["MP_back_ref"]



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


def add_sample(path_a: str, path_b: str, overlay_mode:str = "mean",smoothing:bool = False, rolling:bool = False, background_norm:bool = False, window:int = 1000, flank=500):
    """Reads .csv file, calculates mean over all rows and divides by trimmed mean.

    Args:
        path_a (str): [path to targets]
        path_b (str): [path to background]
        overlay_mode (str): [Mode of operation, e.g. mean, median]

    Returns:
        [Pandas Object]: [Pandas object containing normalized/processed data]
    """

    if window % 2 == 0:
        fstart=int(window/2+1)
        fstop=int(-window/2)
    elif window % 2 == 1:
        fstart=int(window/2-0.5+1)
        fstop=int(-window/2+0.5)
    
    if overlay_mode.lower() == "mean":
        sample_a = pd.read_csv(path_a, header=None).iloc[:,1:].mean()
        if background_norm:
            sample_b = pd.read_csv(path_b, header=None).iloc[:,1:].mean(axis=1)
            sample = pd.DataFrame(sample_a / stats.trim_mean(sample_b, 0.1), columns=["value"])
        else:
            sample = pd.DataFrame(sample_a, columns=["value"])
        sample["position"] = calculate_flanking_regions(len(sample))
        sample=sample.set_index("position")
    elif overlay_mode.lower() == "median":
        sample_a = pd.read_csv(path_a, header=None).iloc[:,1:].median()
        if background_norm:
            sample_b = pd.read_csv(path_b, header=None).iloc[:,1:].median(axis=1)
            sample = pd.DataFrame(sample_a / stats.trim_mean(sample_b, 0.1), columns=["value"])
        else:
            sample = pd.DataFrame(sample_a, columns=["value"])
        sample["position"] = calculate_flanking_regions(len(sample))
        sample=sample.set_index("position")
    elif overlay_mode.lower() == "confidence":
        sample_a = pd.read_csv(path_a, header=None).iloc[:,1:].T
        if background_norm:
            sample_b = pd.read_csv(path_b, header=None).iloc[:,1:].mean(axis=1)
            sample = sample_a / stats.trim_mean(sample_b, 0.1)
        else:
            sample = sample_a
        sample["position"] = calculate_flanking_regions(len(sample))
        sample = sample.set_index("position")
        
    else:
        raise ValueError(f"{overlay_mode} is not a valid keyword.")

    if smoothing:
        if rolling:
            sample = sample.apply(lambda x:savgol_filter(x,window_length=win_len, polyorder=poly)) - sample.apply(lambda x:savgol_filter(x,window_length=21, polyorder=2)).rolling(1000, center=True).median() 
        else:
            sample = sample.apply(lambda x:savgol_filter(x,window_length=win_len, polyorder=poly))
    else:
        if rolling:
            sample = sample - sample.rolling(1000, center=True).median()
    
    sample = sample.iloc[fstart:fstop,:]
        
    if overlay_mode.lower() == "confidence":
        sample = sample.melt(ignore_index=False, var_name="sample_nr")
        
    return sample


# load tables containing position specific scores for all defined target regions
# average over all regions per sample and substract the trimmed mean to normalise


data_sample = pd.DataFrame(add_sample(WPS, WPS_back,overlay_mode,smoothing,rolling,background_norm))
sys.stderr.write("WPS [%s]: %s %s\n"%(sample_ID,COV, COV_back))
if(edge_norm == True):
    edge_mean = pd.concat([data_sample.head(100), data_sample.tail(100)]).mean()
    av_WPS = pd.DataFrame(data_sample/edge_mean)
else:
    av_WPS = pd.DataFrame(data_sample)
av_WPS.columns = av_WPS.columns.astype(str)
av_WPS.columns.values[-1] = sample_ID
for (ref_ID, WPS_ref, WPS_back_ref) in zip(ref_IDs, WPS_refs, WPS_back_refs):
    sys.stderr.write("WPS [%s]: %s %s\n"%(ref_ID, WPS_ref, WPS_back_ref))
    data_reference = add_sample(WPS_ref, WPS_back_ref,overlay_mode,smoothing,rolling,background_norm)["value"]
    if(edge_norm == True):
        edge_mean = pd.concat([data_reference.head(100), data_reference.tail(100)]).mean()
        av_WPS[ref_ID] = data_reference/edge_mean
    else:
        av_WPS[ref_ID] = data_reference

data_sample = pd.DataFrame(add_sample(COV, COV_back,overlay_mode,smoothing,rolling,background_norm))
sys.stderr.write("COV [%s]: %s %s\n"%(sample_ID, COV, COV_back))
if(edge_norm == True):
    edge_mean = pd.concat([data_sample.head(100), data_sample.tail(100)]).mean()
    av_COV = pd.DataFrame(data_sample/edge_mean)
else:
    av_COV = pd.DataFrame(data_sample)
av_COV.columns = av_COV.columns.astype(str)
av_COV.columns.values[-1] = sample_ID
for (ref_ID, COV_ref, COV_back_ref) in zip(ref_IDs, COV_refs, COV_back_refs):
    sys.stderr.write("COV [%s]: %s %s\n"%(ref_ID, COV_ref, COV_back_ref))
    data_reference = add_sample(COV_ref, COV_back_ref,overlay_mode,smoothing,rolling,background_norm)["value"]
    if(edge_norm == True):
        edge_mean = pd.concat([data_reference.head(100), data_reference.tail(100)]).mean()
        av_COV[ref_ID] = data_reference/edge_mean
    else:
        av_COV[ref_ID] = data_reference

if compute_midpoint_coverage == True:
    data_sample = pd.DataFrame(add_sample(MP, MP_back,overlay_mode,smoothing,rolling,background_norm))
    sys.stderr.write("MP [%s]: %s %s\n"%(sample_ID, MP, MP_back))
    if(edge_norm == True):
        edge_mean = pd.concat([data_sample.head(100), data_sample.tail(100)]).mean()
        av_MP = pd.DataFrame(data_sample/edge_mean)
    else:
        av_MP = pd.DataFrame(data_sample)
    av_MP.columns = av_MP.columns.astype(str)
    av_MP.columns.values[-1] = sample_ID
    for (ref_ID, MP_ref, MP_back_ref) in zip(ref_IDs, MP_refs, MP_back_refs):
        sys.stderr.write("MP [%s]: %s %s\n"%(ref_ID, MP_ref, MP_back_ref))
        data_reference = add_sample(MP_ref, MP_back_ref,overlay_mode,smoothing,rolling,background_norm)["value"]
        if(edge_norm == True):
            edge_mean = pd.concat([data_reference.head(100), data_reference.tail(100)]).mean()
            av_MP[ref_ID] = data_reference/edge_mean
        else:
            av_MP[ref_ID] = data_reference

# create line plots and save to a single pdf

if compute_midpoint_coverage == True:
    if overlay_mode.lower() == "confidence":
        av_WPS_long=av_WPS.reset_index().melt(id_vars=["position", "sample_nr"],  value_name="score", var_name="phenotype").sort_values(by="position")
        av_COV_long=av_COV.reset_index().melt(id_vars=["position", "sample_nr"],  value_name="score", var_name="phenotype").sort_values(by="position")
        av_MP_long=av_MP.reset_index().melt(id_vars=["position", "sample_nr"],  value_name="score", var_name="phenotype").sort_values(by="position")
        Fig_WPS = sns.lineplot(data=av_WPS_long, x="position", y="score", hue="phenotype",)
        plt.suptitle("adjusted WPS: {target} target regions")
        plt.xlabel("Position relative to target site")
        plt.ylabel("normalized WPS")
        plt.close()
        Fig_Cov = sns.lineplot(data=av_COV_long, x="position", y="score", hue="phenotype",)
        plt.suptitle(f"adjusted read coverage: {target} target regions")
        plt.xlabel("Position relative to target site")
        plt.ylabel("normalized read coverage")
        plt.close()
        Fig_MP = sns.lineplot(data=av_MP_long, x="position", y="score", hue="phenotype",)
        plt.suptitle(f"adjusted midpoint coverage: {target} target regions")
        plt.xlabel("Position relative to target site")
        plt.ylabel("normalized midpoint coverage")
        plt.close()
    else:
        Fig_WPS = sns.lineplot(data=av_WPS)
        plt.suptitle(f"adjusted WPS: {target} target regions")
        plt.xlabel("Position relative to target site")
        plt.ylabel("normalized WPS")
        plt.close()
        Fig_Cov = sns.lineplot(data=av_COV)
        plt.suptitle(f"adjusted read coverage: {target} target regions")
        plt.xlabel("Position relative to target site")
        plt.ylabel("normalized read coverage")
        plt.close()
        Fig_MP = sns.lineplot(data=av_MP)
        plt.suptitle(f"adjusted midpoint coverage: {target} target regions")
        plt.xlabel("Position relative to target site")
        plt.ylabel("normalized midpoint coverage")
        plt.close()

    with PdfPages(outfile) as pdf:
        pdf.savefig(Fig_WPS.get_figure())
        pdf.savefig(Fig_Cov.get_figure())
        pdf.savefig(Fig_MP.get_figure())

else:
    if overlay_mode.lower() == "confidence":
        av_WPS_long=av_WPS.reset_index().melt(id_vars=["position", "sample_nr"],  value_name="score", var_name="phenotype").sort_values(by="position")
        av_COV_long=av_COV.reset_index().melt(id_vars=["position", "sample_nr"],  value_name="score", var_name="phenotype").sort_values(by="position")
        Fig_WPS = sns.lineplot(data=av_WPS_long, x="position", y="score", hue="phenotype",)
        plt.suptitle("adjusted WPS: {target} target regions")
        plt.xlabel("Position relative to target site")
        plt.ylabel("normalized WPS")
        plt.close()
        Fig_Cov = sns.lineplot(data=av_COV_long, x="position", y="score", hue="phenotype",)
        plt.suptitle(f"adjusted read coverage: {target} target regions")
        plt.xlabel("Position relative to target site")
        plt.ylabel("normalized read coverage")
        plt.close()
    else:
        Fig_WPS = sns.lineplot(data=av_WPS)
        plt.suptitle(f"adjusted WPS: {target} target regions")
        plt.xlabel("Position relative to target site")
        plt.ylabel("normalized WPS")
        plt.close()
        Fig_Cov = sns.lineplot(data=av_COV)
        plt.suptitle(f"adjusted read coverage: {target} target regions")
        plt.xlabel("Position relative to target site")
        plt.ylabel("normalized read coverage")
        plt.close()

    with PdfPages(outfile) as pdf:
        pdf.savefig(Fig_WPS.get_figure())
        pdf.savefig(Fig_Cov.get_figure())
