from snakemake.utils import validate
import pandas as pd


configfile: "config/config.yml" # "config/config_testing.yml" # "config/config_components_bugfix3.yml"


include: "snakefile_GE_analysis.smk"

validate(config, schema="workflow/schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
samples.index.names = ["sample_id"]
validate(samples, schema="workflow/schemas/samples.schema.yaml")

rule target:
    input:
        expand(expand("results/plots/unsupervised/{ID}/heatmaps/FFT_{{min}}-{{max}}_heatmap.png",
            ID=samples["ID"].unique() ),
            zip,
            min=[x[0] for x in config["unsupervised"]["frequencies"]],
            max=[x[1] for x in config["unsupervised"]["frequencies"]]),
        expand(expand("results/plots/unsupervised/{ID}/clustermaps/FFT_{{min}}-{{max}}_clustermap.png",
            ID=samples["ID"].unique() ),
            zip,
            min=[x[0] for x in config["unsupervised"]["frequencies"]],
            max=[x[1] for x in config["unsupervised"]["frequencies"]]),
        expand(expand("results/plots/unsupervised/{ID}/HDBSCAN/FFT_{{min}}-{{max}}_UMAP_projection_HDBSCAN.pdf",
            ID=samples["ID"].unique(), ),
            zip,
            min=[x[0] for x in config["unsupervised"]["frequencies"]],
            max=[x[1] for x in config["unsupervised"]["frequencies"]]),
        expand(expand("results/plots/unsupervised/{ID}/kmeans/FFT_{clust_TYPE}_{{min}}-{{max}}_UMAP_projection_kmeans.pdf",
            ID=samples["ID"].unique(),
            clust_TYPE=["all"]),
            zip,
            min=[x[0] for x in config["unsupervised"]["frequencies"]],
            max=[x[1] for x in config["unsupervised"]["frequencies"]],
            ),
        expand(expand("results/tables/unsupervised/{ID}/kmeans/FFT_{clust_TYPE}_{{min}}-{{max}}_kmeans_metrics.tsv",
            ID=samples["ID"].unique(),
            clust_TYPE=["all"]),
            zip,
            min=[x[0] for x in config["unsupervised"]["frequencies"]],
            max=[x[1] for x in config["unsupervised"]["frequencies"]],
            ),



rule heatmaps:
    input:
        FFT_tables=lambda wc: expand(expand("results/intermediate/{{ID}}/FFT_table/transcriptanno-{SAMPLE}-FFT_table.{GENOME}.tsv",
            zip,
            SAMPLE=samples["sample"].loc[samples["ID"] == wc.ID],
            GENOME=samples["genome_build"].loc[samples["ID"] == wc.ID],
        ),ID=wc.ID,)
    output:
         plot="results/plots/unsupervised/{ID}/heatmaps/FFT_{min}-{max}_heatmap.png"# Plots with heatmaps for freqs
    params:
        min="{min}",
        max="{max}",
        input_IDs=lambda wc: expand("{SAMPLE}", SAMPLE=samples["sample"].loc[samples["ID"] == wc.ID])
    conda: "workflow/envs/unsupervised.yml"
    script:
        "workflow/scripts/unsupervised/heatmaps.py"

rule clustermaps:
    input:
        FFT_tables=lambda wc: expand(expand("results/intermediate/{{ID}}/FFT_table/transcriptanno-{SAMPLE}-FFT_table.{GENOME}.tsv",
            zip,
            SAMPLE=samples["sample"].loc[samples["ID"] == wc.ID],
            GENOME=samples["genome_build"].loc[samples["ID"] == wc.ID],
        ),ID=wc.ID,)
    output:
         plot="results/plots/unsupervised/{ID}/clustermaps/FFT_{min}-{max}_clustermap.png"# Plots with heatmaps for freqs
    params:
        min="{min}",
        max="{max}",
        input_IDs=lambda wc: expand("{SAMPLE}", SAMPLE=samples["sample"].loc[samples["ID"] == wc.ID])
    conda: "workflow/envs/unsupervised.yml"
    script:
        "workflow/scripts/unsupervised/clustermaps.py"

rule kmeans_clustering:
    input:
        FFT_tables=lambda wc: expand(expand("results/intermediate/{{ID}}/FFT_table/transcriptanno-{SAMPLE}-FFT_table.{GENOME}.tsv",
            zip,
            SAMPLE=samples["sample"].loc[samples["ID"] == wc.ID],
            GENOME=samples["genome_build"].loc[samples["ID"] == wc.ID],
        ),ID=wc.ID,)
    output:
        plots="results/plots/unsupervised/{ID}/kmeans/FFT_{clust_TYPE}_{min}-{max}_UMAP_projection_kmeans.pdf",
        metrics = "results/tables/unsupervised/{ID}/kmeans/FFT_{clust_TYPE}_{min}-{max}_kmeans_metrics.tsv"
    params:
        min="{min}",
        max="{max}",
        input_IDs=lambda wc: expand("{SAMPLE}", SAMPLE=samples["sample"].loc[samples["ID"] == wc.ID]),
        umap_epochs=300,
        umap_comps=10,
        n_clust = config["unsupervised"]["kmeans"]["n_clusters"],
        kn_init = 20,
        clust_type = "{clust_TYPE}",
        metric_ID = "{ID}_{min}-{max}"
    conda: "workflow/envs/unsupervised.yml"
    script:
        "workflow/scripts/unsupervised/kmeans_clustering.py"

rule HDBSCAN_clustering:
    input:
        FFT_tables=lambda wc: expand(expand("results/intermediate/{{ID}}/FFT_table/transcriptanno-{SAMPLE}-FFT_table.{GENOME}.tsv",
            zip,
            SAMPLE=samples["sample"].loc[samples["ID"] == wc.ID],
            GENOME=samples["genome_build"].loc[samples["ID"] == wc.ID],
        ),ID=wc.ID,)
    output:
        plots="results/plots/unsupervised/{ID}/HDBSCAN/FFT_{min}-{max}_UMAP_projection_HDBSCAN.pdf",
        metrics = "results/tables/unsupervised/{ID}/HDBSCAN/FFT_{min}-{max}_UMAP_projection_HDBSCAN_metrics.tsv"
    params:
        min="{min}",
        max="{max}",
        input_IDs=lambda wc:expand("{SAMPLE}", SAMPLE=samples["sample"].loc[samples["ID"] == wc.ID]),
        umap_epochs=300,
        umap_comps=10,
        HDBSCAN_min_clustsize=2,
        metric_ID = "{ID}_{min}-{max}"
    conda: "workflow/envs/unsupervised.yml"
    script:
        "workflow/scripts/unsupervised/HDBSCAN_clustering.py"