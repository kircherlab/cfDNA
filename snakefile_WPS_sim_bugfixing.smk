from snakemake.utils import validate
import pandas as pd

configfile: "config/config_components_bugfix2.yml"


validate(config, schema="workflow/schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
samples.index.names = ["sample_id"]
validate(samples, schema="workflow/schemas/samples.schema.yaml")
regions = pd.read_csv(config["regions"], sep="\t").set_index("target", drop=False)
regions.index.names = ["region_id"]
validate(regions, schema="workflow/schemas/regions.schema.yaml")


def get_WPS_ref(sample):
    return expand(
        "results/bugfix/{{ID}}/table/sample/target/{ref_SAMPLE}-{{target_region}}_target_WPS.csv",
        ref_SAMPLE=samples["ref_samples"][sample].split(","),
    )


def get_COV_ref(sample):
    return expand(
        "results/bugfix/{{ID}}/table/sample/target/{ref_SAMPLE}-{{target_region}}_target_COV.csv",
        ref_SAMPLE=samples["ref_samples"][sample].split(","),
    )

def get_WPS_sim_ref(sample):
    return expand(
        "results/bugfix/{{ID}}/table/simulations/target/{ref_SAMPLE}-{{target_region}}_target_sim_WPS.csv",
        ref_SAMPLE=samples["ref_samples"][sample].split(","),
    )


def get_COV_sim_ref(sample):
    return expand(
        "results/bugfix/{{ID}}/table/simulations/target/{ref_SAMPLE}-{{target_region}}_target_sim_COV.csv",
        ref_SAMPLE=samples["ref_samples"][sample].split(","),
    )


def get_STARTS_ref(sample):
    return expand(
        "results/bugfix/{{ID}}/table/sample/target/{ref_SAMPLE}-{{target_region}}_target_STARTS.csv",
        ref_SAMPLE=samples["ref_samples"][sample].split(","),
    )


def get_WPS_background_ref(sample):
    return expand(
        "results/bugfix/{{ID}}/table/sample/background/{ref_SAMPLE}-{{target_region}}_background_WPS.csv",
        ref_SAMPLE=samples["ref_samples"][sample].split(","),
    )


def get_COV_background_ref(sample):
    return expand(
        "results/bugfix/{{ID}}/table/sample/background/{ref_SAMPLE}-{{target_region}}_background_COV.csv",
        ref_SAMPLE=samples["ref_samples"][sample].split(","),
    )

##
def get_WPS_sim_background_ref(sample):
    return expand(
        "results/bugfix/{{ID}}/table/simulations/background/{ref_SAMPLE}-{{target_region}}_background_sim_WPS.csv",
        ref_SAMPLE=samples["ref_samples"][sample].split(","),
    )


def get_COV_sim_background_ref(sample):
    return expand(
        "results/bugfix/{{ID}}/table/simulations/background/{ref_SAMPLE}-{{target_region}}_background_sim_COV.csv",
        ref_SAMPLE=samples["ref_samples"][sample].split(","),
    )


def get_STARTS_background_ref(sample):
    return expand(
        "results/bugfix/{{ID}}/table/sample/background/{ref_SAMPLE}-{{target_region}}_background_STARTS.csv",
        ref_SAMPLE=samples["ref_samples"][sample].split(","),
    )



def get_length(wildcards):
    input=checkpoints.exclude_blacklist.get(ID=wildcards.ID, target_region=wildcards.target_region).output[0]
    #print(input)
    df = pd.read_csv(input, sep="\t", header=None)
    length = df[2] - df[1]
    #print(length)
    return length[0]

def get_num(wildcards):
    input=checkpoints.exclude_blacklist.get(ID=wildcards.ID,  target_region=wildcards.target_region).output[0]
    #print(input)
    df = pd.read_csv(input, sep="\t", header=None)
    num = len(df)
    #print(num)
    return num

rule all:
    input:
        expand(
            "results/bugfix/{ID}/regions/target_region/{target_region}_blacklist-excluded.bed",
            ID=samples["ID"],
            target_region=regions["target"],
        ),
        expand("results/bugfix/{ID}/regions/background_region/{target_region}_background_regions.bed",
                ID=samples["ID"],
                target_region=regions["target"],
        ),
        expand("results/bugfix/{ID}/simulations/sample/{SAMPLE}-{target_region}_sim.sorted.bam.bai",
            SAMPLE=samples["sample"],
            ID=samples["ID"],
            target_region=regions["target"],
            kmer=config["simulations"]["kmer"],
        ),
        expand(
            "results/bugfix/{ID}/table/sample/target/{SAMPLE}-{target_region}_target_WPS.csv",
            SAMPLE=samples["sample"],
            ID=samples["ID"],
            target_region=regions["target"],
        ),
        expand(
            "results/bugfix/{ID}/table/sample/target/{SAMPLE}-{target_region}_target_COV.csv",
            SAMPLE=samples["sample"],
            ID=samples["ID"],
            target_region=regions["target"],
        ),
        expand(
            "results/bugfix/{ID}/table/sample/target/{SAMPLE}-{target_region}_target_STARTS.csv",
            SAMPLE=samples["sample"],
            ID=samples["ID"],
            target_region=regions["target"],
        ),
        expand(
            "results/bugfix/{ID}/table/sample/background/{SAMPLE}-{target_region}_background_WPS.csv",
            SAMPLE=samples["sample"],
            ID=samples["ID"],
            target_region=regions["target"],
        ),
        expand(
            "results/bugfix/{ID}/table/sample/background/{SAMPLE}-{target_region}_background_COV.csv",
            SAMPLE=samples["sample"],
            ID=samples["ID"],
            target_region=regions["target"],
        ),
        expand(
            "results/bugfix/{ID}/table/sample/background/{SAMPLE}-{target_region}_background_STARTS.csv",
            SAMPLE=samples["sample"],
            ID=samples["ID"],
            target_region=regions["target"],
        ),
        expand(
            "results/bugfix/{ID}/table/simulations/target/{SAMPLE}-{target_region}_target_sim_WPS.csv",
            SAMPLE=samples["sample"],
            ID=samples["ID"],
            target_region=regions["target"],
        ),
        expand(
            "results/bugfix/{ID}/table/simulations/target/{SAMPLE}-{target_region}_target_sim_COV.csv",
            SAMPLE=samples["sample"],
            ID=samples["ID"],
            target_region=regions["target"],
        ),
        expand(
            "results/bugfix/{ID}/table/simulations/target/{SAMPLE}-{target_region}_target_sim_STARTS.csv",
            SAMPLE=samples["sample"],
            ID=samples["ID"],
            target_region=regions["target"],
        ),
        expand(
            "results/bugfix/{ID}/table/simulations/background/{SAMPLE}-{target_region}_background_sim_WPS.csv",
            SAMPLE=samples["sample"],
            ID=samples["ID"],
            target_region=regions["target"],
        ),
        expand(
            "results/bugfix/{ID}/table/simulations/background/{SAMPLE}-{target_region}_background_sim_COV.csv",
            SAMPLE=samples["sample"],
            ID=samples["ID"],
            target_region=regions["target"],
        ),
        expand(
            "results/bugfix/{ID}/table/simulations/background/{SAMPLE}-{target_region}_background_sim_STARTS.csv",
            SAMPLE=samples["sample"],
            ID=samples["ID"],
            target_region=regions["target"],
        ),
        #expand(
        #    "results/plots/overlays/{ID}/{target_region}.{SAMPLE}_overlays_sim.pdf",
        #    SAMPLE=samples["sample"],
        #    ID=samples["ID"],
        #    target_region=regions["target"],
        #    kmer=config["simulations"]["kmer"],
        #),

checkpoint exclude_blacklist:
    input:
        Region_file = lambda wildcards: regions["path"][wildcards.target_region],
        blacklist = config["blacklist_GRCH37"]
    output:
        "results/bugfix/{ID}/regions/target_region/{target_region}_blacklist-excluded.bed"
    conda: "workflow/envs/read_preprocessing.yml"
    shell:
        """
        bedtools intersect -v -a {input.Region_file} -b {input.blacklist} > {output}
        """

#rule generate_random_background:
checkpoint generate_random_background:
    input:
        #Region_file = lambda wildcards: regions["path"][wildcards.target_region],
        region="results/bugfix/{ID}/regions/target_region/{target_region}_blacklist-excluded.bed",
        genome=config["GRCh37_genome"],
        universal_blacklist=config["universal_blacklist_GRCH37"],
    output:
        "results/bugfix/{ID}/regions/background_region/{target_region}_background_regions.bed",
    params:
        length = get_length, #lambda wildcards, input: get_length #(input.region),
        num = get_num #lambda wildcads, input: get_num #(input.region)
    conda:"workflow/envs/background.yml"
    shell:
        """
        bedtools random -n {params.num} -l {params.length} -g {input.genome} | \
        bedtools shuffle -i stdin -g {input.genome} -excl {input.universal_blacklist} -noOverlapping \
        > {output}
        """


rule combine_target_background_region:
    input:
        target="results/bugfix/{ID}/regions/target_region/{target_region}_blacklist-excluded.bed",
        background = "results/bugfix/{ID}/regions/background_region/{target_region}_background_regions.bed"
    output: 
        "results/bugfix/{ID}/regions/{target_region}_combined.bed"
    shell:
        "cat {input.target} {input.background} > {output}"


rule simulate_reads:
    input:
        lenDist=config["simulation_values"]["lenDist"],
        fwdKMerGenome=config["simulation_values"]["fwdKMerGenome"],
        fwdPKMers=config["simulation_values"]["fwdPKMers"],
        fwdMKMers=config["simulation_values"]["fwdMKMers"],
        revPKMers=config["simulation_values"]["revPKMers"],
        revMKMers=config["simulation_values"]["revMKMers"],
        regionfile_chrom = "results/bugfix/{ID}/regions/{target_region}_combined.bed",
    output: 
        "results/bugfix/{ID}/simulations/sample/{SAMPLE}-{target_region}_sim.bam" #Do I need the sample ID here?
    params:
        sampling_repeats = config["simulations"]["sampling_repeats"],
        fasta = config["GRCh37"]
    conda: "workflow/envs/cfDNA_sim_chrom.yml"
    shell:
        """
        workflow/scripts/simulations/simulate_reads_bed.py -d {input.lenDist} \
        --fwdKMerGenome={input.fwdKMerGenome} \
        --revKMerGenome={input.fwdKMerGenome} \
        --fwdPKMers={input.fwdPKMers} \
        --fwdMKMers={input.fwdMKMers} \
        --revPKMers={input.revPKMers} \
        --revMKMers={input.revMKMers} \
        -f {params.fasta} \
        -s {params.sampling_repeats} -r {input.regionfile_chrom} \
        -o {output};
        """


rule bam_sorting:
    input:
        "results/bugfix/{ID}/simulations/sample/{SAMPLE}-{target_region}_sim.bam"
    output:
        "results/bugfix/{ID}/simulations/sample/{SAMPLE}-{target_region}_sim.sorted.bam"
    conda: "workflow/envs/read_preprocessing.yml"
    threads: config["threads"]
    shell:
        """samtools sort -T /fast/users/roeners_c/scratch/{wildcards.SAMPLE} -o {output} -O bam --threads {threads} {input}"""

rule index_and_stats:
    input:
        "results/bugfix/{ID}/simulations/sample/{SAMPLE}-{target_region}_sim.sorted.bam"
    output:
        "results/bugfix/{ID}/simulations/sample/{SAMPLE}-{target_region}_sim.sorted.bam.bai",
        "results/bugfix/{ID}/simulations/sample/{SAMPLE}-{target_region}_sim.sorted.bam_stats",
        "results/bugfix/{ID}/simulations/sample/{SAMPLE}-{target_region}_sim.sorted.byChromCounts.txt"
    params:
        allChrom = "results/bugfix/{ID}/simulations/sample/{SAMPLE}-{target_region}_sim.sorted"
    conda: "workflow/envs/read_preprocessing.yml"#"workflow/envs/cfDNA.yml"
    threads: config["threads"]
    shell:
        """
        samtools index -@ {threads} {input}
        samtools flagstat -@ {threads} {input} > {params.allChrom}.bam_stats
        samtools view -@ {threads} -F _2 {input} | cut -f 3 | uniq -c > {params.allChrom}.byChromCounts.txt
        """


rule extract_counts:
    input:
        #target=lambda wildcards: regions["path"][wildcards.target_region],
        target = "results/bugfix/{ID}/regions/target_region/{target_region}_blacklist-excluded.bed",
        BAMFILE=lambda wildcards: samples["path"][wildcards.SAMPLE],
    output:
        WPS="results/bugfix/{ID}/table/sample/target/{SAMPLE}-{target_region}_target_WPS.csv",
        COV="results/bugfix/{ID}/table/sample/target/{SAMPLE}-{target_region}_target_COV.csv",
        STARTS="results/bugfix/{ID}/table/sample/target/{SAMPLE}-{target_region}_target_STARTS.csv",
    params:
        minRL=config["minRL"],
        maxRL=config["maxRL"],
        out_pre="results/bugfix/{ID}/table/sample/target/{SAMPLE}-{target_region}_target_%s.csv",
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
        #target=lambda wildcards: regions["path"][wildcards.target_region],
        target = "results/bugfix/{ID}/regions/background_region/{target_region}_background_regions.bed",
        BAMFILE=lambda wildcards: samples["path"][wildcards.SAMPLE],
    output:
        WPS="results/bugfix/{ID}/table/sample/background/{SAMPLE}-{target_region}_background_WPS.csv",
        COV="results/bugfix/{ID}/table/sample/background/{SAMPLE}-{target_region}_background_COV.csv",
        STARTS="results/bugfix/{ID}/table/sample/background/{SAMPLE}-{target_region}_background_STARTS.csv",
    params:
        minRL=config["minRL"],
        maxRL=config["maxRL"],
        out_pre="results/bugfix/{ID}/table/sample/background/{SAMPLE}-{target_region}_background_%s.csv",
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

rule extract_counts_sim:
    input:
        #target=lambda wildcards: regions["path"][wildcards.target_region],
        target = "results/bugfix/{ID}/regions/target_region/{target_region}_blacklist-excluded.bed",
        BAMFILE="results/bugfix/{ID}/simulations/sample/{SAMPLE}-{target_region}_sim.sorted.bam",
        BAM_index="results/bugfix/{ID}/simulations/sample/{SAMPLE}-{target_region}_sim.sorted.bam.bai",
    output:
        WPS="results/bugfix/{ID}/table/simulations/target/{SAMPLE}-{target_region}_target_sim_WPS.csv",
        COV="results/bugfix/{ID}/table/simulations/target/{SAMPLE}-{target_region}_target_sim_COV.csv",
        STARTS="results/bugfix/{ID}/table/simulations/target/{SAMPLE}-{target_region}_target_sim_STARTS.csv",
    params:
        minRL=config["minRL"],
        maxRL=config["maxRL"],
        out_pre="results/bugfix/{ID}/table/simulations/target/{SAMPLE}-{target_region}_target_sim_%s.csv",
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


rule extract_counts_background_sim:
    input:
        #target=lambda wildcards: regions["path"][wildcards.target_region],
        target = "results/bugfix/{ID}/regions/background_region/{target_region}_background_regions.bed",
        BAMFILE="results/bugfix/{ID}/simulations/sample/{SAMPLE}-{target_region}_sim.sorted.bam",
        BAM_index="results/bugfix/{ID}/simulations/sample/{SAMPLE}-{target_region}_sim.sorted.bam.bai",
    output:
        WPS="results/bugfix/{ID}/table/simulations/background/{SAMPLE}-{target_region}_background_sim_WPS.csv",
        COV="results/bugfix/{ID}/table/simulations/background/{SAMPLE}-{target_region}_background_sim_COV.csv",
        STARTS="results/bugfix/{ID}/table/simulations/background/{SAMPLE}-{target_region}_background_sim_STARTS.csv",
    params:
        minRL=config["minRL"],
        maxRL=config["maxRL"],
        out_pre="results/bugfix/{ID}/table/simulations/background/{SAMPLE}-{target_region}_background_sim_%s.csv",
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

# add options for ALL and NONE -> plot all samples, only target sample -> get_WPS

rule plot_overlays:
    input:
        WPS="results/bugfix/{ID}/table/sample/target/{SAMPLE}-{target_region}_target_WPS.csv",
        WPS_ref=lambda wildcards: get_WPS_ref(wildcards.SAMPLE),
        WPS_back="results/bugfix/{ID}/table/sample/background/{SAMPLE}-{target_region}_background_WPS.csv",
        WPS_back_ref=lambda wildcards: get_WPS_background_ref(wildcards.SAMPLE),
        WPS_sim="results/bugfix/{ID}/table/simulations/target/{SAMPLE}-{target_region}_target_sim_WPS.csv",
        WPS_sim_ref=lambda wildcards: get_WPS_sim_ref(wildcards.SAMPLE),
        WPS_sim_back="results/bugfix/{ID}/table/simulations/background/{SAMPLE}-{target_region}_background_sim_WPS.csv",
        WPS_sim_back_ref=lambda wildcards: get_WPS_sim_background_ref(wildcards.SAMPLE),
        COV="results/bugfix/{ID}/table/sample/target/{SAMPLE}-{target_region}_target_COV.csv",
        COV_ref=lambda wildcards: get_COV_ref(wildcards.SAMPLE),
        COV_back="results/bugfix/{ID}/table/sample/background/{SAMPLE}-{target_region}_background_COV.csv",
        COV_back_ref=lambda wildcards: get_COV_background_ref(wildcards.SAMPLE),
        COV_sim="results/bugfix/{ID}/table/simulations/target/{SAMPLE}-{target_region}_target_sim_COV.csv",
        COV_sim_ref=lambda wildcards: get_COV_sim_ref(wildcards.SAMPLE),
        COV_sim_back="results/bugfix/{ID}/table/simulations/background/{SAMPLE}-{target_region}_background_sim_COV.csv",
        COV_sim_back_ref=lambda wildcards: get_COV_sim_background_ref(wildcards.SAMPLE),
    output:
        "results/plots/overlays/{ID}/{target_region}.{SAMPLE}_overlays_sim.pdf",
    params:
        target="{target_region}",
        sample="{SAMPLE}",
        ref_IDs=lambda wildcards: samples["ref_samples"][wildcards.SAMPLE].split(","),
    conda:
        "workflow/envs/overlays.yml"
    script:
        "workflow/scripts/simulations/overlays_sim.py"