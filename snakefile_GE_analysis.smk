from snakemake.utils import validate
import pandas as pd

configfile: "config/config.yml"
validate(config, schema="workflow/schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
samples.index.names = ["sample_id"]
validate(samples, schema="workflow/schemas/samples.schema.yaml")

def get_refsample(refsample):
    return "results/intermediate/body/fft_summaries/fft_{}-{}_WPS.tsv.gz".format(refsample,samples["genome_build"].loc[samples["sample"] == refsample].values[0])

rule all:
    input:
        expand("results/intermediate/transcriptAnno/transcriptAnno-{GENOME}.103.filtered.tsv.gz",
             GENOME=samples["genome_build"].unique()
        ),
        expand("results/intermediate/transcriptAnno/transcriptAnno-{GENOME}.103.body.bed.gz",
                GENOME=samples["genome_build"].unique()),
        expand("results/intermediate/transcriptAnno/transcriptAnno_background-{GENOME}.103.body.bed.gz",
                GENOME=samples["genome_build"].unique()),
        #expand("results/intermediate/body/fft_summaries/fft_{SAMPLE}-{GENOME}_WPS.tsv.gz",
        #        zip,
        #        GENOME=samples["genome_build"],
        #        SAMPLE=samples["sample"]),
        #expand("results/plots/{ID}/{tissue}_allFreq_correlation_plot.pdf",
        #        tissue=config["tissue"],
        #        ID=samples["ID"]),
        #expand("results/tables/{ID}/Ave193-199bp_correlation.pdf",
        #        ID=samples["ID"]),
        #expand("results/tables/{ID}/{refSample}_Ave193-199bp_correlation_rank.pdf",
        #        refSample = config["refSample"],
        #        ID=samples["ID"])



rule join: 
    input:
        transcriptAnno= lambda wc: config[wc.GENOME]["transcriptAnno"],
        proteinAtlas= expand("resources/protein_atlas/RNAtable{SOURCE}.tsv.gz",
                            SOURCE=config["proteinAtlas"]),
    output:
        filteredTranscriptAnno="results/intermediate/transcriptAnno/transcriptAnno-{GENOME}.103.filtered.tsv.gz"
    conda: "workflow/envs/cfDNA.yml"
    shell:
        """
        set +o pipefail;
        (zcat {input.transcriptAnno} | \
        head -n 1; join -t"$(echo -e "\\t")" \
        <(\zcat {input.transcriptAnno}| tail -n +2 | sort -k1,1 ) \
        <(zcat {input.proteinAtlas} | tail -n +2 |cut -f 1 | sort )) | \
        gzip -c > {output.filteredTranscriptAnno}
        """

rule prep:
    input:
        transcriptAnno="results/intermediate/transcriptAnno/transcriptAnno-{GENOME}.103.filtered.tsv.gz"
    output:
        body="results/intermediate/transcriptAnno/transcriptAnno-{GENOME}.103.body.bed.gz"
    conda: "workflow/envs/cfDNA.yml"
    shell:
        """
        zcat {input.transcriptAnno} | tail -n +2 | \
        awk 'BEGIN{{ FS="\\t"; OFS="\\t" }}{{ if ($5 == "+") {{ print $2,$3-1,$3-1+10000,$1,0,$5 }} else {{ print $2,$4-1-10000,$4-1,$1,0,$5 }} }}'| \
        gzip -c > {output.body}
        """

rule generate_random_background:
    input:
        region=(
            "results/intermediate/transcriptAnno/transcriptAnno-{GENOME}.103.body.bed.gz"
        ),
        genome=lambda wildcards: config[wildcards.GENOME]["genome_autosomes"], 
        gap=lambda wildcards: config[wildcards.GENOME]["UCSC_gap"],
    output:
        "results/intermediate/transcriptAnno/transcriptAnno_background-{GENOME}.103.body.bed.gz"
    params:
        length=10000, #lambda wildcards, input: get_length(input.region)
    conda:
        "workflow/envs/background.yml"
    shell:
        """
        bedtools random -n 1000 -l {params.length} -g {input.genome} | \
        bedtools shuffle -i stdin -g {input.genome} -excl {input.gap} -noOverlapping | \
        gzip -c > {output}
        """

rule extract_counts:
    input:
        body="results/intermediate/transcriptAnno/transcriptAnno-{GENOME}.103.body.tsv.gz",
        BAMFILE= lambda wildcards: samples["path"][wildcards.SAMPLE]
    output:
        "results/intermediate/body/fft_summaries/fft_{SAMPLE}-{GENOME}_WPS.tsv.gz"
    params:
        minRL=config["minRL"],
        maxRL=config["maxRL"],
        out_pre="results/tmp/body/{SAMPLE}-{GENOME}/block_%s.tsv.gz"
    conda: "workflow/envs/cfDNA.yml"
    shell:
        """
        mkdir -p results/tmp/body/{wildcards.SAMPLE}-{wildcards.GENOME}

        workflow/scripts/expression_analysis/extractReadStartsFromBAM_Region_WPS.py \
        --minInsert={params.minRL} \
        --maxInsert={params.maxRL} \
        -i {input.body} \
        -o {params.out_pre} {input.BAMFILE}

        mkdir -p results/tmp/body/{wildcards.SAMPLE}-{wildcards.GENOME}/fft

        ( cd results/tmp/body/{wildcards.SAMPLE}-{wildcards.GENOME}; ls block_*.tsv.gz ) | \
        xargs -n 500 Rscript workflow/scripts/expression_analysis/fft_path.R \
        results/tmp/body/{wildcards.SAMPLE}-{wildcards.GENOME}/ \
        results/tmp/body/{wildcards.SAMPLE}-{wildcards.GENOME}/fft

        mkdir -p results/tmp/body/fft_summaries/

        workflow/scripts/expression_analysis/convert_files.py \
        -a {input.body} \
        -t results/tmp/ \
        -r results/intermediate/ \
        -p body \
        -i {wildcards.SAMPLE}-{wildcards.GENOME}
        
        #rm -fR results/tmp/body/{wildcards.SAMPLE}-{wildcards.GENOME}/fft
        """


rule correlation_plots:
    input:
        samples = expand("results/intermediate/body/fft_summaries/fft_{SAMPLE}-{GENOME}_WPS.tsv.gz",
                        zip,
                        GENOME=samples["genome_build"],
                        SAMPLE=samples["sample"]),        
        proteinAtlas= expand("resources/protein_atlas/RNAtable{SOURCE}.tsv.gz",
                            SOURCE=config["proteinAtlas"]),
        labels = expand("resources/protein_atlas/labels_{SOURCE}.tsv",
                        SOURCE=config["proteinAtlas"]),
    output:
        allFreq = "results/plots/{ID}/{tissue}_allFreq_correlation_plot.pdf",
    conda: "workflow/envs/cfDNA.yml"
    params: 
        tissue = "{tissue}"
    script:
        "workflow/scripts/expression_analysis/snakemake_correlation_plots.R"


rule correlation_table:
    input:
        samples = expand("results/intermediate/body/fft_summaries/fft_{SAMPLE}-{GENOME}_WPS.tsv.gz",
                        zip,
                        GENOME=samples["genome_build"],
                        SAMPLE=samples["sample"]),
        proteinAtlas= expand("resources/protein_atlas/RNAtable{SOURCE}.tsv.gz",
                            SOURCE=config["proteinAtlas"]),
        labels = expand("resources/protein_atlas/labels_{SOURCE}.tsv",
                        SOURCE=config["proteinAtlas"]),
    output:
        aveCor = "results/tables/{ID}/Ave193-199bp_correlation.pdf",
    conda: "workflow/envs/cfDNA.yml"
    script:
        "workflow/scripts/expression_analysis/snakemake_correlation_table.R"


rule rank_correlation_table:
    input:
        samples = expand("results/intermediate/body/fft_summaries/fft_{SAMPLE}-{GENOME}_WPS.tsv.gz",
                        zip,
                        GENOME=samples["genome_build"],
                        SAMPLE=samples["sample"]),
        proteinAtlas= expand("resources/protein_atlas/RNAtable{SOURCE}.tsv.gz",
                            SOURCE=config["proteinAtlas"]),
        labels = expand("resources/protein_atlas/labels_{SOURCE}.tsv",
                        SOURCE=config["proteinAtlas"]),
    output:
        aveCorRank = "results/tables/{ID}/{refSample}_Ave193-199bp_correlation_rank.pdf",
    conda: "workflow/envs/cfDNA.yml"
    params: 
        refSample = get_refsample(config["refSample"])
    script:
        "workflow/scripts/expression_analysis/snakemake_correlation_rank_table.R"
