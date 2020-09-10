from snakemake.utils import validate
import pandas as pd

configfile: "config/config.yml"
validate(config, schema="workflow/schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
samples.index.names = ["sample_id"]
validate(samples, schema="workflow/schemas/samples.schema.yaml")

rule all:
    input:
        "resources/annotations/transcriptAnno-GRCh37.75.upstream.tsv",
        "resources/annotations/transcriptAnno-GRCh37.75.downstream.tsv",
        "resources/annotations/transcriptAnno-GRCh37.75.body.tsv",
        expand("results/intermediate/body/fft_summaries/fft_{SAMPLE}_WPS.tsv.gz",
                SAMPLE=samples["sample"]),
        expand("results/plots/{ID}/{tissue}_allFreq_correlation_plot.pdf",
                tissue=config["tissue"],
                ID=samples["ID"]),
        expand("results/tables/{ID}/Ave193-199bp_correlation.pdf",
                ID=samples["ID"]),
        expand("results/tables/{ID}/{refSample}_Ave193-199bp_correlation_rank.pdf",
                refSample = config["refSample"],
                ID=samples["ID"])


rule prep:
    input:
        transcriptAnno="resources/transcriptAnno-GRCh37.75.tsv.gz"
    output:
        upstream="resources/annotations/transcriptAnno-GRCh37.75.upstream.tsv",
        downstream="resources/annotations/transcriptAnno-GRCh37.75.downstream.tsv",
        body="resources/annotations/transcriptAnno-GRCh37.75.body.tsv"
    conda: "workflow/envs/cfDNA.yml"
    shell:
        """
        zcat {input.transcriptAnno} | awk 'BEGIN{{ FS="\\t"; OFS="\\t" }}{{ if ($5 == "+") {{ print $1,$2,$3-10000,$3,$5 }} else {{ print $1,$2,$4,$4+10000,$5 }} }}' > {output.upstream}
        zcat {input.transcriptAnno} | awk 'BEGIN{{ FS="\\t"; OFS="\\t" }}{{ if ($5 == "+") {{ print $1,$2,$4,$4+10000,$5 }} else {{ print $1,$2,$3-10000,$3,$5 }} }}' > {output.downstream}
        zcat {input.transcriptAnno} | awk 'BEGIN{{ FS="\\t"; OFS="\\t" }}{{ if ($5 == "+") {{ print $1,$2,$3-1,$3-1+10000,$5 }} else {{ print $1,$2,$4-1-10000,$4-1,$5 }} }}' >{output.body}
        """

rule extract_counts:
    input:
        body="resources/annotations/transcriptAnno-GRCh37.75.body.tsv",
        BAMFILE= lambda wildcards: samples["path"][wildcards.SAMPLE]
    output:
        "results/intermediate/body/fft_summaries/fft_{SAMPLE}_WPS.tsv.gz"
    params:
        minRL=config["minRL"],
        maxRL=config["maxRL"],
        out_pre="results/tmp/body/{SAMPLE}/block_%s.tsv.gz"
    conda: "workflow/envs/cfDNA.yml"
    shell:
        """
        mkdir -p results/tmp/body/{wildcards.SAMPLE}

        workflow/scripts/expression_analysis/extractReadStartsFromBAM_Region_WPS.py \
        --minInsert={params.minRL} \
        --maxInsert={params.maxRL} \
        -i {input.body} \
        -o {params.out_pre} {input.BAMFILE}

        mkdir -p results/tmp/body/{wildcards.SAMPLE}/fft

        ( cd results/tmp/body/{wildcards.SAMPLE}; ls block_*.tsv.gz ) | \
        xargs -n 500 Rscript workflow/scripts/expression_analysis/fft_path.R \
        results/tmp/body/{wildcards.SAMPLE}/ \
        results/tmp/body/{wildcards.SAMPLE}/fft

        mkdir -p results/tmp/body/fft_summaries/

        workflow/scripts/expression_analysis/convert_files.py \
        -a {input.body} \
        -t results/tmp/ \
        -r results/intermediate/ \
        -p body \
        -i {wildcards.SAMPLE}
        
        rm -fR results/tmp/body/{wildcards.SAMPLE}/fft
        """


rule correlation_plots:
    input:
        samples = expand("results/intermediate/body/fft_summaries/fft_{SAMPLE}_WPS.tsv.gz",
                          SAMPLE=samples["sample"]),
        proteinAtlas = "resources/protein_atlas/RNAtable.tsv",
        labels = "resources/labels.txt",
    output:
        allFreq = "results/plots/{ID}/{tissue}_allFreq_correlation_plot.pdf",
    conda: "workflow/envs/cfDNA.yml"
    params: 
        IDs = expand("{SAMPLE}",SAMPLE=samples["sample"]),
        WPSprefix = "results/intermediate/body/fft_summaries/fft_%s_WPS.tsv.gz",
        tissue = "{tissue}"
    script:
        "workflow/scripts/expression_analysis/snakemake_correlation_plots.R"


rule correlation_table:
    input:
        samples = expand("results/intermediate/body/fft_summaries/fft_{SAMPLE}_WPS.tsv.gz",
                          SAMPLE=samples["sample"]),
        proteinAtlas = "resources/protein_atlas/RNAtable.tsv",
        labels = "resources/labels.txt",
    output:
        aveCor = "results/tables/{ID}/Ave193-199bp_correlation.pdf",
    conda: "workflow/envs/cfDNA.yml"
    params: 
        IDs = expand("{SAMPLE}",SAMPLE=samples["sample"]),
        WPSprefix = "results/intermediate/body/fft_summaries/fft_%s_WPS.tsv.gz",
    script:
        "workflow/scripts/expression_analysis/snakemake_correlation_table.R"


rule rank_correlation_table:
    input:
        samples = expand("results/intermediate/body/fft_summaries/fft_{SAMPLE}_WPS.tsv.gz",
                          SAMPLE=samples["sample"]),
        proteinAtlas = "resources/protein_atlas/RNAtable.tsv",
        labels = "resources/labels.txt",
    output:
        aveCorRank = "results/tables/{ID}/{refSample}_Ave193-199bp_correlation_rank.pdf",
    conda: "workflow/envs/cfDNA.yml"
    params: 
        IDs = expand("{SAMPLE}",SAMPLE=samples["sample"]),
        WPSprefix = "results/intermediate/body/fft_summaries/fft_%s_WPS.tsv.gz",
        refSample = config["refSample"]
    script:
        "workflow/scripts/expression_analysis/snakemake_correlation_rank_table.R"
