from snakemake.utils import validate
import pandas as pd

configfile: "config/config_DHS_analysis.yml"
validate(config, schema="workflow/schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
samples.index.names = ["sample_id"]
validate(samples, schema="workflow/schemas/samples.schema.yaml")
regions = pd.read_csv(config["regions"], sep="\t").set_index("target", drop=False)
regions.index.names = ["region_id"]
validate(regions, schema="workflow/schemas/regions.schema.yaml")

rule all:
    input:
        #"resources/annotations/transcriptAnno-GRCh37.75.upstream.tsv",
        #"resources/annotations/transcriptAnno-GRCh37.75.downstream.tsv",
        #"resources/annotations/transcriptAnno-GRCh37.75.body.tsv",
        expand("results/intermediate/{ID}/DHS_analysis/{region}-transcriptAnno-GRCh38.upstream.tsv",
                ID=samples["ID"],
                region=regions["target"]),
        expand("results/intermediate/{ID}/DHS_analysis/{region}-transcriptAnno-GRCh38.downstream.tsv",
                ID=samples["ID"],
                region=regions["target"]),
        expand("results/intermediate/{ID}/DHS_analysis/{region}-transcriptAnno-GRCh38.body.tsv",
                ID=samples["ID"],
                region=regions["target"]),
        expand("results/plots/{ID}/{region}/{tissue}_allFreq_correlation_plot.pdf",
                tissue=config["tissue"],
                ID=samples["ID"],
                region=regions["target"]),
        expand("results/tables/{ID}/{region}/Ave193-199bp_correlation.pdf",
                ID=samples["ID"],
                region=regions["target"]),
        expand("results/tables/{ID}/{region}/{refSample}_Ave193-199bp_correlation_rank.pdf",
                refSample = config["refSample"],
                ID=samples["ID"],
                region=regions["target"])


rule prep:
    input:
        transcriptAnno=lambda wildcards: regions["path"][wildcards.region],
        #transcriptAnno="resources/transcriptAnno-GRCh37.75.tsv.gz"
    output:
        upstream="results/intermediate/{ID}/DHS_analysis/{region}-transcriptAnno-GRCh38.upstream.tsv",
        downstream="results/intermediate/{ID}/DHS_analysis/{region}-transcriptAnno-GRCh38.downstream.tsv",
        body="results/intermediate/{ID}/DHS_analysis/{region}-transcriptAnno-GRCh38.body.tsv"
    conda: "workflow/envs/cfDNA.yml"
    shell:
        """
        zcat {input.transcriptAnno} | awk 'BEGIN{{ FS="\\t"; OFS="\\t" }}{{ if ($5 == "+") {{ print $1,$2,$3-10000,$3,$5 }} else {{ print $1,$2,$4,$4+10000,$5 }} }}' > {output.upstream}
        zcat {input.transcriptAnno} | awk 'BEGIN{{ FS="\\t"; OFS="\\t" }}{{ if ($5 == "+") {{ print $1,$2,$4,$4+10000,$5 }} else {{ print $1,$2,$3-10000,$3,$5 }} }}' > {output.downstream}
        zcat {input.transcriptAnno} | awk 'BEGIN{{ FS="\\t"; OFS="\\t" }}{{ if ($5 == "+") {{ print $1,$2,$3-1,$3-1+10000,$5 }} else {{ print $1,$2,$4-1-10000,$4-1,$5 }} }}' >{output.body}
        """

rule extract_counts:
    input:
        body="results/intermediate/{ID}/DHS_analysis/{region}-transcriptAnno-GRCh38.body.tsv",
        BAMFILE= lambda wildcards: samples["path"][wildcards.SAMPLE]
    output:
        "results/intermediate/{ID}/DHS_analysis/body/{region}/fft_summaries/fft_{SAMPLE}_WPS.tsv.gz"
    params:
        minRL=config["minRL"],
        maxRL=config["maxRL"],
        out_pre="results/tmp/body/{region}/{SAMPLE}/block_%s.tsv.gz"
    conda: "workflow/envs/cfDNA.yml"
    shell:
        """
        mkdir -p results/tmp/body/{wildcards.region}/{wildcards.SAMPLE}

        workflow/scripts/expression_analysis/extractReadStartsFromBAM_Region_WPS.py \
        --minInsert={params.minRL} \
        --maxInsert={params.maxRL} \
        -i {input.body} \
        -o {params.out_pre} {input.BAMFILE}

        mkdir -p results/tmp/{wildcards.ID}/DHS_analysis/body/{wildcards.region}/{wildcards.SAMPLE}/fft

        ( cd results/tmp/body/{wildcards.region}/{wildcards.SAMPLE}; ls block_*.tsv.gz ) | \
        xargs -n 500 Rscript workflow/scripts/expression_analysis/fft_path.R \
        results/tmp/body/{wildcards.region}/{wildcards.SAMPLE} \
        results/tmp/{wildcards.ID}/DHS_analysis/body/{wildcards.region}/{wildcards.SAMPLE}/fft

        mkdir -p results/tmp/body/fft_summaries/

        workflow/scripts/expression_analysis/convert_files.py \
        -a {input.body} \
        -t results/tmp/ \
        -r results/intermediate/ \
        -p {wildcards.ID}/DHS_analysis/body/{wildcards.region} \
        -i {wildcards.SAMPLE}
        
        rm -fR results/tmp/body/{wildcards.SAMPLE}/fft
        """
# change rm -fR to:         results/tmp/{wildcards.ID}/DHS_analysis/body/{wildcards.region}/{wildcards.SAMPLE}/fft

rule correlation_plots:
    input:
        samples = expand("results/intermediate/{{ID}}/DHS_analysis/body/{{region}}/fft_summaries/fft_{SAMPLE}_WPS.tsv.gz",
                          SAMPLE=samples["sample"]),
        proteinAtlas = "resources/protein_atlas/RNAtable.tsv",
        labels = "resources/labels.txt",
    output:
        allFreq = "results/plots/{ID}/{region}/{tissue}_allFreq_correlation_plot.pdf",
    conda: "workflow/envs/cfDNA.yml"
    params: 
        IDs = expand("{SAMPLE}",SAMPLE=samples["sample"]),
        WPSprefix = "results/intermediate/{ID}/DHS_analysis/body/{region}/fft_summaries//fft_%s_WPS.tsv.gz", # wo/wie wird der prefix hier genutzt
        tissue = "{tissue}"
    script:
        "workflow/scripts/expression_analysis/snakemake_correlation_plots.R"


rule correlation_table:
    input:
        samples = expand("results/intermediate/{{ID}}/DHS_analysis/body/{{region}}/fft_summaries/fft_{SAMPLE}_WPS.tsv.gz",
                          SAMPLE=samples["sample"]),
        proteinAtlas = "resources/protein_atlas/RNAtable.tsv",
        labels = "resources/labels.txt",
    output:
        aveCor = "results/tables/{ID}/{region}/Ave193-199bp_correlation.pdf",
    conda: "workflow/envs/cfDNA.yml"
    params: 
        IDs = expand("{SAMPLE}",SAMPLE=samples["sample"]),
        WPSprefix = "results/intermediate/{ID}/DHS_analysis/body/{region}/fft_summaries//fft_%s_WPS.tsv.gz",
    script:
        "workflow/scripts/expression_analysis/snakemake_correlation_table.R"


rule rank_correlation_table:
    input:
        samples = expand("results/intermediate/{{ID}}/DHS_analysis/body/{{region}}/fft_summaries/fft_{SAMPLE}_WPS.tsv.gz",
                          SAMPLE=samples["sample"]),
        proteinAtlas = "resources/protein_atlas/RNAtable.tsv",
        labels = "resources/labels.txt",
    output:
        aveCorRank = "results/tables/{ID}/{region}/{refSample}_Ave193-199bp_correlation_rank.pdf",
    conda: "workflow/envs/cfDNA.yml"
    params: 
        IDs = expand("{SAMPLE}",SAMPLE=samples["sample"]),
        WPSprefix = "results/intermediate/{ID}/DHS_analysis/body/{region}/fft_summaries//fft_%s_WPS.tsv.gz",
        refSample = config["refSample"]
    script:
        "workflow/scripts/expression_analysis/snakemake_correlation_rank_table.R"
