from snakemake.utils import validate
import pandas as pd


configfile: "config/config.yml"

validate(config, schema="workflow/schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
samples.index.names = ["sample_id"]
validate(samples, schema="workflow/schemas/samples.schema.yaml")


def get_length(input):
    df = pd.read_csv(input, sep="\t", header=None)
    print(df.head())
    length = df[3] - df[2]
    return length[0]


rule all:
    input:
        "resources/annotations/transcriptAnno-GRCh37.75.body.tsv",
        expand(
            "results/intermediate/{ID}/table/transcriptanno_{SAMPLE}_WPS_normalized.tsv",
            SAMPLE=samples["sample"],
            ID=samples["ID"],
        ),
        expand(
            "results/intermediate/{ID}/table/transcriptanno_{SAMPLE}_COV_normalized.tsv",
            SAMPLE=samples["sample"],
            ID=samples["ID"],
        ),
        expand(
            "results/intermediate/{ID}/FFT_table/transcriptanno_{SAMPLE}_FFT_table.tsv",
            SAMPLE=samples["sample"],
            ID=samples["ID"],
        ),
        expand("results/intermediate/{ID}/corr_plots/{tissue}_allFreq_correlation_plot.pdf",
                tissue=config["tissue"],
                ID=samples["ID"].unique()),
        expand(
            "results/intermediate/{ID}/corr_tables/Ave193-199bp_correlation.pdf",
            ID=samples["ID"].unique(),
        ), 
        expand("results/intermediate/{ID}/corr_tables/{refSample}_Ave193-199bp_correlation_rank.pdf",
                refSample = config["refSample"],
                ID=samples["ID"].unique())


rule prep:
    input:
        transcriptAnno="resources/transcriptAnno-GRCh37.75.tsv.gz",
    output:
        body="results/intermediate/{ID}/annotations/transcriptAnno-GRCh37.75.body.bed",
    conda:
        "workflow/envs/cfDNA.yml"
    shell:
        """
        zcat {input.transcriptAnno} | awk 'BEGIN{{ FS="\\t"; OFS="\\t" }}{{ if ($5 == "+") {{ print $2,$3-1,$3-1+10000,$1,0,$5 }} else {{ print $2,$4-1-10000,$4-1,$1,0,$5 }} }}' >{output.body}
        """


rule generate_random_background:
    input:
        region=(
            "results/intermediate/{ID}/annotations/transcriptAnno-GRCh37.75.body.bed"
        ),
        genome=config["GRCh37_genome"],
        gap=config["UCSC_gap_GRCH37"],
    output:
        "results/intermediate/{ID}/background_region/transcriptanno_background_regions.bed",
    params:
        length=10000, #lambda wildcards, input: get_length(input.region)
    conda:
        "workflow/envs/background.yml"
    shell:
        """
        bedtools random -n 1000 -l {params.length} -g {input.genome} | \
        bedtools shuffle -i stdin -g {input.genome} -excl {input.gap} -noOverlapping \
        > {output}
        """


rule extract_counts:
    input:
        target=(
            "results/intermediate/{ID}/annotations/transcriptAnno-GRCh37.75.body.bed"
        ),
        BAMFILE=lambda wildcards: samples["path"][wildcards.SAMPLE],
    output:
        WPS="results/intermediate/{ID}/table/transcriptanno_{SAMPLE}_WPS.csv",
        COV="results/intermediate/{ID}/table/transcriptanno_{SAMPLE}_COV.csv",
        STARTS="results/intermediate/{ID}/table/transcriptanno_{SAMPLE}_STARTS.csv",
    params:
        minRL=config["minRL"],
        maxRL=config["maxRL"],
        out_pre="results/intermediate/{ID}/table/transcriptanno_{SAMPLE}_%s.csv",
    conda:
        "workflow/envs/cfDNA.yml"
    shell:
        """
        workflow/scripts/WPS/extractFromBAM_RegionBed_WPS_Cov.py \
        --minInsert={params.minRL} \
        --maxInsert={params.maxRL} \
        -i {input.target} \
        -o {params.out_pre} {input.BAMFILE}
        """


rule extract_counts_background:
    input:
        background="results/intermediate/{ID}/background_region/transcriptanno_background_regions.bed",
        BAMFILE=lambda wildcards: samples["path"][wildcards.SAMPLE],
    output:
        WPS="results/intermediate/{ID}/background_region/table/transcriptanno_{SAMPLE}_WPS.background.csv",
        COV="results/intermediate/{ID}/background_region/table/transcriptanno_{SAMPLE}_COV.background.csv",
        STARTS="results/intermediate/{ID}/background_region/table/transcriptanno_{SAMPLE}_STARTS.background.csv",
    params:
        minRL=config["minRL"],
        maxRL=config["maxRL"],
        out_pre="results/intermediate/{ID}/background_region/table/transcriptanno_{SAMPLE}_%s.background.csv",
    conda:
        "workflow/envs/cfDNA.yml"
    shell:
        """
        workflow/scripts/WPS/extractFromBAM_RegionBed_WPS_Cov.py \
        --minInsert={params.minRL} \
        --maxInsert={params.maxRL} \
        -i {input.background} \
        -o {params.out_pre} {input.BAMFILE}
        """


rule normalize_WPS:
    input:
        target_WPS="results/intermediate/{ID}/table/transcriptanno_{SAMPLE}_WPS.csv",
        background_WPS="results/intermediate/{ID}/background_region/table/transcriptanno_{SAMPLE}_WPS.background.csv",
        target_COV="results/intermediate/{ID}/table/transcriptanno_{SAMPLE}_COV.csv",
        background_COV="results/intermediate/{ID}/background_region/table/transcriptanno_{SAMPLE}_COV.background.csv",
    output:
        output_WPS=(
            "results/intermediate/{ID}/table/transcriptanno_{SAMPLE}_WPS_normalized.tsv"
        ),
        output_COV=(
            "results/intermediate/{ID}/table/transcriptanno_{SAMPLE}_COV_normalized.tsv"
        ),
    conda:
        "workflow/envs/cfDNA.yml"
    script:
        """workflow/scripts/expression_analysis/normalize.py"""


rule FFT_table:
    input:
        normalized_WPS=(
            "results/intermediate/{ID}/table/transcriptanno_{SAMPLE}_WPS_normalized.tsv"
        ),
    output:
        FFT_table=(
            "results/intermediate/{ID}/FFT_table/transcriptanno_{SAMPLE}_FFT_table.tsv"
        ),
    conda:
        "workflow/envs/overlays.yml"
    script:
        """workflow/scripts/expression_analysis/fft_table.py"""


# rule extract_counts:
#    input:
#        body="resources/annotations/transcriptAnno-GRCh37.75.body.tsv",
#        BAMFILE= lambda wildcards: samples["path"][wildcards.SAMPLE]
#    output:
#        "results/intermediate/body/fft_summaries/fft_{SAMPLE}_WPS.tsv.gz"
#    params:
#        minRL=config["minRL"],
#        maxRL=config["maxRL"],
#        out_pre="results/tmp/body/{SAMPLE}/block_%s.tsv.gz"
#    conda: "workflow/envs/cfDNA.yml"
#    shell:
#        """
#        mkdir -p results/tmp/body/{wildcards.SAMPLE}
#
#        workflow/scripts/expression_analysis/extractReadStartsFromBAM_Region_WPS.py \
#        --minInsert={params.minRL} \
#        --maxInsert={params.maxRL} \
#        -i {input.body} \
#        -o {params.out_pre} {input.BAMFILE}
#
#        mkdir -p results/tmp/body/{wildcards.SAMPLE}/fft
#
#        ( cd results/tmp/body/{wildcards.SAMPLE}; ls block_*.tsv.gz ) | \
#        xargs -n 500 Rscript workflow/scripts/expression_analysis/fft_path.R \
#        results/tmp/body/{wildcards.SAMPLE}/ \
#        results/tmp/body/{wildcards.SAMPLE}/fft
#
#        mkdir -p results/tmp/body/fft_summaries/
#
#        workflow/scripts/expression_analysis/convert_files.py \
#        -a {input.body} \
#        -t results/tmp/ \
#        -r results/intermediate/ \
#        -p body \
#        -i {wildcards.SAMPLE}
#
#        rm -fR results/tmp/body/{wildcards.SAMPLE}/fft
#        """


rule correlation_plots:
    input:
        samples=lambda wildcards: expand(
            "results/intermediate/{ID}/FFT_table/transcriptanno_{SAMPLE}_FFT_table.tsv",
            SAMPLE=samples["sample"].loc[samples.ID == wildcards.ID],
            ID=wildcards.ID,
        ),
        proteinAtlas="resources/protein_atlas/RNAtable.tsv",
        labels="resources/labels.txt",
    output:
        allFreq="results/intermediate/{ID}/corr_plots/{tissue}_allFreq_correlation_plot.pdf",
    conda:
        "workflow/envs/cfDNA.yml"
    params:
        IDs=lambda wildcards: expand(
            "{SAMPLE}", SAMPLE=samples["sample"].loc[samples.ID == wildcards.ID]
        ),
        WPSprefix="results/intermediate/%s/FFT_table/transcriptanno_%s_FFT_table.tsv",
        tissue="{tissue}",
    script:
        "workflow/scripts/expression_analysis/snakemake_correlation_plots.R"


rule correlation_table:
    input:
        samples=lambda wildcards: expand(
            "results/intermediate/{ID}/FFT_table/transcriptanno_{SAMPLE}_FFT_table.tsv",
            SAMPLE=samples["sample"].loc[samples.ID == wildcards.ID],
            ID=wildcards.ID,
        ),
        proteinAtlas="resources/protein_atlas/RNAtable.tsv",
        labels="resources/labels.txt",
    output:
        aveCor="results/intermediate/{ID}/corr_tables/Ave193-199bp_correlation.pdf",
    conda:
        "workflow/envs/cfDNA.yml"
    params:
        samples=lambda wildcards: expand(
            "{SAMPLE}", SAMPLE=samples["sample"].loc[samples.ID == wildcards.ID]
        ),
        ID="{ID}",
        WPSprefix="results/intermediate/%s/FFT_table/transcriptanno_%s_FFT_table.tsv",
    script:
        "workflow/scripts/expression_analysis/snakemake_correlation_table.R"


rule rank_correlation_table:
    input:
        samples=lambda wildcards: expand(
            "results/intermediate/{ID}/FFT_table/transcriptanno_{SAMPLE}_FFT_table.tsv",
            SAMPLE=samples["sample"].loc[samples.ID == wildcards.ID],
            ID=wildcards.ID,
        ),
        proteinAtlas="resources/protein_atlas/RNAtable.tsv",
        labels="resources/labels.txt",
    output:
        aveCorRank="results/intermediate/{ID}/corr_tables/{refSample}_Ave193-199bp_correlation_rank.pdf",
    conda:
        "workflow/envs/cfDNA.yml"
    params:
        IDs=lambda wildcards: expand(
            "{SAMPLE}", SAMPLE=samples["sample"].loc[samples.ID == wildcards.ID]
        ),
        WPSprefix="results/intermediate/%s/FFT_table/transcriptanno_%s_FFT_table.tsv",
        refSample=config["refSample"],
    script:
        "workflow/scripts/expression_analysis/snakemake_correlation_rank_table.R"
