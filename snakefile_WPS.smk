from snakemake.utils import validate
import pandas as pd


configfile: "config/config.yml"


validate(config, schema="workflow/schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
samples.index.names = ["sample_id"]
validate(samples, schema="workflow/schemas/samples.schema.yaml")
regions = pd.read_csv(config["regions"], sep="\t").set_index("target", drop=False)
regions.index.names = ["region_id"]
validate(regions, schema="workflow/schemas/regions.schema.yaml")


def get_WPS_ref(sample):
    return expand(
        "results/intermediate/{{ID}}/table/{{target_region}}.{ref_SAMPLE}_WPS.csv",
        ref_SAMPLE=samples["ref_samples"][sample].split(","),
    )


def get_COV_ref(sample):
    return expand(
        "results/intermediate/{{ID}}/table/{{target_region}}.{ref_SAMPLE}_COV.csv",
        ref_SAMPLE=samples["ref_samples"][sample].split(","),
    )


def get_STARTS_ref(sample):
    return expand(
        "results/intermediate/{{ID}}/table/{{target_region}}.{ref_SAMPLE}_STARTS.csv",
        ref_SAMPLE=samples["ref_samples"][sample].split(","),
    )

def get_WPS_background_ref(sample):
    return expand(
        "results/intermediate/{{ID}}/background_region/table/{{target_region}}.{ref_SAMPLE}_WPS.background.csv",
        ref_SAMPLE=samples["ref_samples"][sample].split(","),
    )


def get_COV_background_ref(sample):
    return expand(
        "results/intermediate/{{ID}}/background_region/table/{{target_region}}.{ref_SAMPLE}_COV.background.csv",
        ref_SAMPLE=samples["ref_samples"][sample].split(","),
    )


def get_STARTS_background_ref(sample):
    return expand(
        "results/intermediate/{{ID}}/background_region/table/{{target_region}}.{ref_SAMPLE}_STARTS.background.csv",
        ref_SAMPLE=samples["ref_samples"][sample].split(","),
    )



def get_length(input):
    df = pd.read_csv(input, sep="\t", header=None)
    length = df[2] - df[1]
    return length[0]


rule all:
    input:
        expand(
            "results/intermediate/{ID}/table/{target_region}.{SAMPLE}_WPS.csv",
            SAMPLE=samples["sample"],
            ID=samples["ID"],
            target_region=regions["target"],
        ),
        expand(
            "results/intermediate/{ID}/table/{target_region}.{SAMPLE}_COV.csv",
            SAMPLE=samples["sample"],
            ID=samples["ID"],
            target_region=regions["target"],
        ),
        expand(
            "results/intermediate/{ID}/table/{target_region}.{SAMPLE}_STARTS.csv",
            SAMPLE=samples["sample"],
            ID=samples["ID"],
            target_region=regions["target"],
        ),
        expand(
            "results/intermediate/{ID}/background_region/{target_region}.{SAMPLE}_background_regions.bed",
            SAMPLE=samples["sample"],
            ID=samples["ID"],
            target_region=regions["target"],
        ),
        expand(
            "results/plots/overlays/{ID}/{target_region}.{SAMPLE}_overlays.pdf",
            SAMPLE=samples["sample"],
            ID=samples["ID"],
            target_region=regions["target"],
        ),
        


rule generate_random_background:
    input:
        region=lambda wildcards: regions["path"][wildcards.target_region],
        genome=config["GRCh37_genome"],
        gap=config["UCSC_gap"],
    output:
        "results/intermediate/{ID}/background_region/{target_region}.{SAMPLE}_background_regions.bed",
    params:
        length = lambda wildcards, input: get_length(input.region)
    conda:"workflow/envs/background.yml"
    shell:
        """
        bedtools random -n 1000 -l {params.length} -g {input.genome} | \
        bedtools shuffle -i stdin -g {input.genome} -excl {input.gap} -noOverlapping \
        > {output}
        """


rule extract_counts:
    input:
        target=lambda wildcards: regions["path"][wildcards.target_region],
        BAMFILE=lambda wildcards: samples["path"][wildcards.SAMPLE],
    output:
        WPS="results/intermediate/{ID}/table/{target_region}.{SAMPLE}_WPS.csv",
        COV="results/intermediate/{ID}/table/{target_region}.{SAMPLE}_COV.csv",
        STARTS="results/intermediate/{ID}/table/{target_region}.{SAMPLE}_STARTS.csv",
    params:
        minRL=config["minRL"],
        maxRL=config["maxRL"],
        out_pre="results/intermediate/{ID}/table/{target_region}.{SAMPLE}_%s.csv",
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
        background="results/intermediate/{ID}/background_region/{target_region}.{SAMPLE}_background_regions.bed",
        BAMFILE=lambda wildcards: samples["path"][wildcards.SAMPLE],
    output:
        WPS="results/intermediate/{ID}/background_region/table/{target_region}.{SAMPLE}_WPS.background.csv",
        COV="results/intermediate/{ID}/background_region/table/{target_region}.{SAMPLE}_COV.background.csv",
        STARTS="results/intermediate/{ID}/background_region/table/{target_region}.{SAMPLE}_STARTS.background.csv",
    params:
        minRL=config["minRL"],
        maxRL=config["maxRL"],
        out_pre="results/intermediate/{ID}/background_region/table/{target_region}.{SAMPLE}_%s.background.csv",
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


rule plot_overlays:
    input:
        WPS="results/intermediate/{ID}/table/{target_region}.{SAMPLE}_WPS.csv",
        WPS_ref=lambda wildcards: get_WPS_ref(wildcards.SAMPLE),
        COV="results/intermediate/{ID}/table/{target_region}.{SAMPLE}_COV.csv",
        COV_ref=lambda wildcards: get_COV_ref(wildcards.SAMPLE),
        WPS_back="results/intermediate/{ID}/background_region/table/{target_region}.{SAMPLE}_WPS.background.csv",
        WPS_back_ref=lambda wildcards: get_WPS_background_ref(wildcards.SAMPLE),
        COV_back="results/intermediate/{ID}/background_region/table/{target_region}.{SAMPLE}_COV.background.csv",
        COV_back_ref=lambda wildcards: get_COV_background_ref(wildcards.SAMPLE),
    output:
        "results/plots/overlays/{ID}/{target_region}.{SAMPLE}_overlays.pdf",
    params:
        target="{target_region}",
        sample="{SAMPLE}",
        ref_IDs=lambda wildcards: samples["ref_samples"][wildcards.SAMPLE].split(","),
    conda:
        "workflow/envs/overlays.yml"
    script:
        "workflow/scripts/WPS/overlays.py"