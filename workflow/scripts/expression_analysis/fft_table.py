import pandas as pd
import numpy as np
from scipy import stats, signal, fft
from statsmodels.tsa.filters.filtertools import recursive_filter

from spec_pgram_mod import spec_pgram

input_file = snakemake.input["normalized_WPS"]
output_file = snakemake.output["FFT_table"]


taps = 1/np.arange(start=5,stop=100,step=4)
def rec_fil(x):
    xtmp = x.iloc[1:301].append(x.iloc[1:])
    ftmp = recursive_filter(xtmp,taps)[300:]
    return ftmp

def calculate_FFT(path):
    indf = pd.read_csv(path, sep="\t", header=None)
    if pd.api.types.is_string_dtype(indf[0]):
        indf = indf.set_index(0)
        indf.index.name="#Region"
    
    #indf = indf.apply(lambda x: signal.savgol_filter(x, window_length=21, polyorder=2), axis=1)
    #indf = pd.DataFrame(indf.to_list(), index=indf.index)
    indf = indf.apply(lambda x: rec_fil(x), axis=1)
    indf = indf.apply(lambda x: x-stats.trim_mean(x,0.1), axis=1)
    specDF = indf.apply(lambda x: spec_pgram(x,pad=0.3,taper=0.3,spans=2,plot=False,detrend=True,demean=True), axis=1)
    freq = specDF.apply(lambda x: 1/x["freq"])
    header = pd.DataFrame(freq.values.tolist(), index=freq.index).apply(lambda x: np.unique(x)).iloc[0,:].astype(int)
    datadf = specDF.apply(lambda x: x["spec"])
    datadf = pd.DataFrame(datadf.values.tolist(), columns=header)
    datadf.index = indf.index
    tmpdf = datadf.loc[:,(datadf.columns >= 120) & (datadf.columns <= 280) ]
    return tmpdf

FFT_df = calculate_FFT(input_file)

#print("FFT_df")
#print(FFT_df.head())

if FFT_df.index.is_object():
    FFT_df.to_csv(output_file, sep="\t")
else:
    FFT_df.to_csv(output_file, sep="\t", index=None)