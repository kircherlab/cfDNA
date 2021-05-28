import pandas as pd
import numpy as np
from scipy import stats


target_WPS = snakemake.input["target_WPS"]
background_WPS = snakemake.input["background_WPS"]

target_COV = snakemake.input["target_COV"]
background_COV = snakemake.input["background_COV"]

output_WPS = snakemake.output["output_WPS"]
output_COV = snakemake.output["output_COV"]

index_choice=None

def normalize_sample(path_a: str, path_b: str):
    """Reads .csv file, calculates mean over all rows and divides by trimmed mean.

    Args:
        path (str): Path to a .csv file

    Returns:
        [type]: [description]
    """
    sample_a = pd.read_csv(path_a, header=None)
    if pd.api.types.is_string_dtype(sample_a[0]):
        sample_a = sample_a.set_index(0)
        sample_a.index.name="ID"
    sample_b = pd.read_csv(path_b, header=None).mean(axis=1)
    sample = sample_a / stats.trim_mean(sample_b, 0.1)
    return sample.round(4)

normalized_WPS = normalize_sample(path_a=target_WPS, path_b=background_WPS)

normalized_COV = normalize_sample(path_a=target_COV, path_b=background_COV)

if normalized_WPS.index.is_object():
    normalized_WPS.to_csv(output_WPS, sep="\t",header=None)
else:
    normalized_WPS.to_csv(output_WPS, sep="\t",header=None, index=None)

if normalized_COV.index.is_object():
    normalized_COV.to_csv(output_COV, sep="\t",header=None)
else:
    normalized_COV.to_csv(output_COV, sep="\t",header=None, index=None)





