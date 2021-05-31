import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.use("pdf")

input_files = snakemake.input["FFT_tables"]
output_plot = snakemake.output["plot"]
input_IDs = snakemake.params["input_IDs"]
min_freq = snakemake.params["min"]
max_freq = snakemake.params["max"]


FFT_df = pd.DataFrame()
for ID,path in zip(input_IDs,input_files):
    print(ID, path)
    indf=pd.read_csv(path,sep="\t",index_col="#Region")
    indf = indf.loc[:,(indf.columns.astype(int) >= int(min_freq)) & (indf.columns.astype(int) <= int(max_freq))]
    FFT_df[ID]=indf.values.flatten()

sns.set_theme(context='notebook', style='darkgrid', palette='deep', font='sans-serif', font_scale=1.5, color_codes=True, rc={'figure.figsize':(18.7,16.27)})#(14.7,12.27)})

FFT_spearman = FFT_df.sort_index(axis=1).corr(method="spearman")

fig_spearman = sns.clustermap(data=FFT_spearman,vmin=-1, vmax=1, center=0, xticklabels=True, yticklabels=True)
plt.title(f'Spearman correlation between FFT values ({min_freq}-{max_freq})', fontsize=20)
plt.setp(fig_spearman.ax_heatmap.get_yticklabels(), rotation=0)
plt.setp(fig_spearman.ax_heatmap.get_xticklabels(), rotation=90)
#plt.yticks(rotation=90)
#plt.xticks(rotation=90) 
plt.tight_layout()


fig_spearman.savefig(output_plot)