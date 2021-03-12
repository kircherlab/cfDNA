from snakemake.utils import validate
import pandas as pd

configfile: "config/config_GE_bugfix.yml"
validate(config, schema="workflow/schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
samples.index.names = ["sample_id"]
validate(samples, schema="workflow/schemas/samples.schema.yaml")

rule all:
    input:
        expand("results/bugfix/GE/plots/{ID}/{tissue}_allFreq_correlation_plot_overlap_meuleman.pdf",
                tissue=config["tissue"],
                ID=samples["ID"]),
        expand("results/bugfix/GE/tables/{ID}/Ave193-199bp_correlation_overlap_meuleman.pdf",
                ID=samples["ID"]),
        expand("results/bugfix/GE/tables/{ID}/{refSample}_Ave193-199bp_correlation_rank_overlap_meuleman.pdf",
                refSample = config["refSample"],
                ID=samples["ID"])


rule correlation_plots:
    input:
        samples = expand("results/bugfix/fft_{SAMPLE}_summary_overlap_meuleman.tsv.gz",
                          SAMPLE=samples["sample"]),
        proteinAtlas = "resources/protein_atlas/RNAtable.tsv",
        labels = "resources/labels.txt",
    output:
        allFreq = "results/bugfix/GE/plots/{ID}/{tissue}_allFreq_correlation_plot_overlap_meuleman.pdf",
    conda: "workflow/envs/cfDNA.yml"
    params: 
        IDs = expand("{SAMPLE}",SAMPLE=samples["sample"]),
        WPSprefix = "results/bugfix/fft_%s_summary_overlap_meuleman.tsv.gz",
        tissue = "{tissue}"
    script:
        "workflow/scripts/expression_analysis/snakemake_correlation_plots.R"


rule correlation_table:
    input:
        samples = expand("results/bugfix/fft_{SAMPLE}_summary_overlap_meuleman.tsv.gz",
                          SAMPLE=samples["sample"]),
        proteinAtlas = "resources/protein_atlas/RNAtable.tsv",
        labels = "resources/labels.txt",
    output:
        aveCor = "results/bugfix/GE/tables/{ID}/Ave193-199bp_correlation_overlap_meuleman.pdf",
    conda: "workflow/envs/cfDNA.yml"
    params: 
        IDs = expand("{SAMPLE}",SAMPLE=samples["sample"]),
        WPSprefix = "results/bugfix/fft_%s_summary_overlap_meuleman.tsv.gz",
    script:
        "workflow/scripts/expression_analysis/snakemake_correlation_table.R"


rule rank_correlation_table:
    input:
        samples = expand("results/bugfix/fft_{SAMPLE}_summary_overlap_meuleman.tsv.gz",
                          SAMPLE=samples["sample"]),
        proteinAtlas = "resources/protein_atlas/RNAtable.tsv",
        labels = "resources/labels.txt",
    output:
        aveCorRank = "results/bugfix/GE/tables/{ID}/{refSample}_Ave193-199bp_correlation_rank_overlap_meuleman.pdf",
    conda: "workflow/envs/cfDNA.yml"
    params: 
        IDs = expand("{SAMPLE}",SAMPLE=samples["sample"]),
        WPSprefix = "results/bugfix/fft_%s_summary_overlap_meuleman.tsv.gz",
        refSample = config["refSample"]
    script:
        "workflow/scripts/expression_analysis/snakemake_correlation_rank_table.R"
