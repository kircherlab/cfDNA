import pandas as pd
import numpy as np
from sklearn import cluster
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score

from umap import UMAP
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import seaborn as sns


input_files = snakemake.input["FFT_tables"]
input_IDs = snakemake.params["input_IDs"]
outfile = snakemake.output["plots"]
out_metrics = snakemake.output["metrics"]

metric_ID = snakemake.params["metric_ID"]

min_freq = snakemake.params["min"]
max_freq = snakemake.params["max"]
umap_epochs = snakemake.params["umap_epochs"]
umap_comps = snakemake.params["umap_comps"]
kn_clust = snakemake.params["n_clust"]
kn_init = snakemake.params["kn_init"]
clust_type = snakemake.params["clust_type"]

FFT_df = pd.DataFrame()
for ID,path in zip(input_IDs,input_files):
    indf=pd.read_csv(path,sep="\t",index_col="#Region")
    indf = indf.loc[:,(indf.columns.astype(int) >= int(min_freq)) & (indf.columns.astype(int) <= int(max_freq))]
    FFT_df[ID]=indf.values.flatten()

#conds = [FFT_df.T.index.str.contains("NPH")==True,
#        (FFT_df.T.index.str.contains("P") == True) & (FFT_df.T.index.str.contains("NPH") == False),
#        (FFT_df.T.index.str.contains("C") == True),
#        (FFT_df.T.index.str.contains("B") == True)]
#
#choices = ["healthy","P","C","B"]

conds = [ (FFT_df.T.index.str.contains("NPH")==True) ,
             ((FFT_df.T.index.str.contains("BH0")==True) | (FFT_df.T.index.str.contains("IH0")==True)),
        (FFT_df.T.index.str.contains("P") == True) & (FFT_df.T.index.str.contains("NPH") == False),
        (FFT_df.T.index.str.contains("C2") == True),
        (FFT_df.T.index.str.contains("B") == True) & (FFT_df.T.index.str.contains("BH") == False),
        (FFT_df.T.index.str.contains("IC") == True)   ]
choices = ["healthy","healthy_Cell","prostate","colon","breast","disease_Cell"]

conds_source = [((FFT_df.T.index.str.contains("BH0")==True) | (FFT_df.T.index.str.contains("IH0")==True)),
            (FFT_df.T.index.str.contains("IC")==True),
             (FFT_df.T.index.str.contains("NPH")==True),
            ((FFT_df.T.index.str.contains("C2_")==True) | (FFT_df.T.index.str.contains("P")==True) | (FFT_df.T.index.str.contains("NPH")==False)),
            ]
choices_source = ["healthy_Cell","disease_Cell","healthy_Graz","disease_Graz"]

conds_bin = conds_bin = [(FFT_df.T.index.str.contains("NPH")==True),
            ((FFT_df.T.index.str.contains("C2_")==True) | (FFT_df.T.index.str.contains("P")==True) | (FFT_df.T.index.str.contains("NPH")==False)),
            ((FFT_df.T.index.str.contains("BH0")==True) | (FFT_df.T.index.str.contains("IH0")==True)),
            (FFT_df.T.index.str.contains("IC")==True)
            ]
choices_bin = ["healthy_Graz","disease_Graz","healthy_Cell","disease_Cell"]



FFT_df_T = FFT_df.T
FFT_df_T["status"] = np.select(conds, choices)
FFT_df_T["status_bin"] = np.select(conds_bin, choices_bin)
FFT_df_T["status_source"] = np.select(conds_source, choices_source)


umap_2d = UMAP(n_components=2, n_epochs=300)

umap_projections = umap_2d.fit_transform(FFT_df.T)

projection = plt.figure(figsize=(15,8))
sns.scatterplot(x=umap_projections[:,0],
                y=umap_projections[:,1],
                hue=FFT_df_T["status"])
plt.title(f"UMAP 2D projection: FFT values ({min_freq}-{max_freq}bp) ground truth", fontsize=20)
plt.xlabel("UMAP_0")
plt.ylabel("UMAP_1")
plt.close()

projection_bin = plt.figure(figsize=(15,8))
sns.scatterplot(x=umap_projections[:,0],
                y=umap_projections[:,1],
                hue=FFT_df_T["status_bin"])
plt.title(f"UMAP 2D projection: FFT values ({min_freq}-{max_freq}bp) binary labels", fontsize=20)
plt.xlabel("UMAP_0")
plt.ylabel("UMAP_1")
plt.close()

projection_source = plt.figure(figsize=(15,8))
sns.scatterplot(x=umap_projections[:,0],
                y=umap_projections[:,1],
                hue=FFT_df_T["status_source"])
plt.title(f"UMAP 2D projection: FFT values ({min_freq}-{max_freq}bp) labels by dataset", fontsize=20)
plt.xlabel("UMAP_0")
plt.ylabel("UMAP_1")
plt.close()


if clust_type == "UMAP":
    clust_input=FFT_df.T
elif clust_type == "all":
    clust_input=FFT_df.T


plot_list=[]
metric_list = list()

for nclust  in kn_clust:
    Kclust_init = cluster.KMeans(n_clusters=int(nclust), init="random", n_init=int(kn_init))
    cluster_all_labels = Kclust_init.fit_predict(clust_input)

    Kclust_plot = plt.figure(figsize=(15,8))
    sns.scatterplot(x=umap_projections[:,0],
                    y=umap_projections[:,1],
                    hue=cluster_all_labels.astype(str))
    plt.title(f"UMAP 2D projection labeled by k-means(n={nclust}): FFT values ({min_freq}-{max_freq}bp)", fontsize=20)
    plt.xlabel("UMAP_0")
    plt.ylabel("UMAP_1")
    plt.close()

    plot_list.append(Kclust_plot)

    clustered = (cluster_all_labels >= 0)

    metric_dict = {"ID":[metric_ID],
                "raw_ARS":[adjusted_rand_score(FFT_df_T["status"], cluster_all_labels)],
                "noOutliers_ARS":[adjusted_rand_score(FFT_df_T["status"][clustered], cluster_all_labels[clustered])],
                "raw_AMI":[adjusted_mutual_info_score(FFT_df_T["status"], cluster_all_labels)],
                "noOutliers_AMI":[adjusted_mutual_info_score(FFT_df_T["status"][clustered], cluster_all_labels[clustered])],
                "raw_ARS_bin":[adjusted_rand_score(FFT_df_T["status_bin"], cluster_all_labels)],
                "noOutliers_ARS_bin":[adjusted_rand_score(FFT_df_T["status_bin"][clustered], cluster_all_labels[clustered])],
                "raw_AMI_bin":[adjusted_mutual_info_score(FFT_df_T["status_bin"], cluster_all_labels)],
                "noOutliers_AMI_bin":[adjusted_mutual_info_score(FFT_df_T["status_bin"][clustered], cluster_all_labels[clustered])],
                "raw_ARS_source":[adjusted_rand_score(FFT_df_T["status_source"], cluster_all_labels)],
                "noOutliers_ARS_source":[adjusted_rand_score(FFT_df_T["status_source"][clustered], cluster_all_labels[clustered])],
                "raw_AMI_source":[adjusted_mutual_info_score(FFT_df_T["status_source"], cluster_all_labels)],
                "noOutliers_AMI_source":[adjusted_mutual_info_score(FFT_df_T["status_source"][clustered], cluster_all_labels[clustered])]
               }
    metric_list.append(metric_dict)

metric_df = pd.DataFrame.from_records(metric_list)

metric_df.to_csv(out_metrics, sep="\t",index=None)

with PdfPages(outfile) as pdf:
    pdf.savefig(projection.get_figure())
    pdf.savefig(projection_bin.get_figure())
    pdf.savefig(projection_source.get_figure())
    for plot in plot_list:
        pdf.savefig(plot)
