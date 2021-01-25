from snakemake.utils import validate
import pandas as pd
import glob

configfile: "config/config_components.yml"

#ruleorder: exclude_blacklist > combine_chromosome_files

validate(config, schema="workflow/schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
samples.index.names = ["sample_id"]
validate(samples, schema="workflow/schemas/samples.schema.yaml")
regions = pd.read_csv(config["regions"], sep="\t").set_index("target", drop=False)
regions.index.names = ["region_id"]
validate(regions, schema="workflow/schemas/regions.schema.yaml")

def get_chroms(ID, target_region):
    BEDCOLS = ["chrom", "chromStart", "chromEnd",
           "name", "score", "strand",
           "thickStart", "thickEnd", "itemRGB",
           "blockCount", "blockSizes", "blockStarts"]
    #print(target_region)
    bedfile=f"results/intermediate/{ID}/regions/{target_region}_blacklist-excluded.bed"
    #print(bedfile)
    #bedfile=regions["path"][target_region]
    bed_df = pd.read_csv(bedfile, sep="\t", header=None)
    colnames = []
    for i in range(len(bed_df.columns)):
        colnames.append(BEDCOLS[i])
    bed_df.columns = colnames
    chroms = bed_df["chrom"].unique().tolist()
    #f_list = expand("results/intermediate/{{ID}}/simulations/tmp/{{SAMPLE}}.{{target_region}}_{{kmer}}mer_sim_{chrom}.bam",#"results/intermediate/{{ID}}/simulations/regions/{{target_region}}/{{target_region}}_{chrom}.bed",
    #                chrom=chroms)
    #print(f_list)
    return expand("results/intermediate/{{ID}}/simulations/tmp/{{SAMPLE}}.{{target_region}}_{{kmer}}mer_sim_{chrom}.bam",chrom=chroms)#"results/intermediate/{{ID}}/simulations/regions/{{target_region}}/{{target_region}}_{chrom}.bed",

def get_chroms_back(ID, target_region):
    BEDCOLS = ["chrom", "chromStart", "chromEnd",
           "name", "score", "strand",
           "thickStart", "thickEnd", "itemRGB",
           "blockCount", "blockSizes", "blockStarts"]
    #print(target_region)
    bedfile=f"results/intermediate/{ID}/regions/{target_region}_blacklist-excluded.bed"
    #print(bedfile)
    #bedfile=regions["path"][target_region]
    bed_df = pd.read_csv(bedfile, sep="\t", header=None)
    colnames = []
    for i in range(len(bed_df.columns)):
        colnames.append(BEDCOLS[i])
    bed_df.columns = colnames
    chroms = bed_df["chrom"].unique().tolist()
    #f_list = expand("results/intermediate/{{ID}}/simulations/tmp/{{SAMPLE}}.{{target_region}}_{{kmer}}mer_sim_{chrom}.bam",#"results/intermediate/{{ID}}/simulations/regions/{{target_region}}/{{target_region}}_{chrom}.bed",
    #                chrom=chroms)
    #print(f_list)
    return expand("results/intermediate/{{ID}}/simulations/tmp/background_region/background_{{SAMPLE}}.{{target_region}}_{{kmer}}mer_sim_{chrom}.bam",chrom=chroms)#"results/intermediate/{{ID}}/simulations/regions/{{target_region}}/{{target_region}}_{chrom}.bed",



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
        "results/intermediate/{{ID}}/background_region/table/sim_{{kmer}}mer_{{target_region}}.{ref_SAMPLE}_WPS.csv",
        ref_SAMPLE=samples["ref_samples"][sample].split(","),
    )


def get_COV_background_ref(sample):
    return expand(
        "results/intermediate/{{ID}}/background_region/table/sim_{{kmer}}mer_{{target_region}}.{ref_SAMPLE}_COV.csv",
        ref_SAMPLE=samples["ref_samples"][sample].split(","),
    )


def get_STARTS_background_ref(sample):
    return expand(
        "results/intermediate/{{ID}}/background_region/table/sim_{{kmer}}mer_{{target_region}}.{ref_SAMPLE}_STARTS.csv",
        ref_SAMPLE=samples["ref_samples"][sample].split(","),
    )



def get_length(input):
    df = pd.read_csv(input, sep="\t", header=None)
    length = df[2] - df[1]
    return length[0]


rule all:
    input:
        expand(
            "results/intermediate/{ID}/regions/{target_region}_blacklist-excluded.bed",
            ID=samples["ID"],
            target_region=regions["target"],
        ),
        expand("results/intermediate/{ID}/simulations/{SAMPLE}.{target_region}_{kmer}mer_allChrom.bam",
            SAMPLE=samples["sample"],
            ID=samples["ID"],
            target_region=regions["target"],
            kmer=config["simulations"]["kmer"],
        ),
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
            "results/intermediate/{ID}/table/sim_{kmer}mer_{target_region}.{SAMPLE}_WPS.csv",
            SAMPLE=samples["sample"],
            ID=samples["ID"],
            target_region=regions["target"],
            kmer=config["simulations"]["kmer"],
        ),
        expand(
            "results/intermediate/{ID}/table/background_region/sim_{kmer}mer_{target_region}.{SAMPLE}_WPS.csv",
            SAMPLE=samples["sample"],
            ID=samples["ID"],
            target_region=regions["target"],
            kmer=config["simulations"]["kmer"],
        ),
        ##expand(
        ##    "results/plots/overlays/{ID}/{target_region}.{SAMPLE}_overlays_sim_{kmer}mer.pdf",
        ##    SAMPLE=samples["sample"],
        ##    ID=samples["ID"],
        #    target_region=regions["target"],
        #    kmer=config["simulations"]["kmer"],
        #),

rule exclude_blacklist:
    input:
        Region_file = lambda wildcards: regions["path"][wildcards.target_region],
        blacklist = config["blacklist_GRCH37"]
    output:
        "results/intermediate/{ID}/regions/{target_region}_blacklist-excluded.bed"
    conda: "workflow/envs/read_preprocessing.yml"
    shell:
        """
        bedtools intersect -v -a {input.Region_file} -b {input.blacklist} > {output}
        """

rule generate_random_background:
    input:
        region="results/intermediate/{ID}/regions/{target_region}_blacklist-excluded.bed",
        genome=config["GRCh37_genome"],
        gap=config["UCSC_gap_GRCH37"],
    output:
        "results/intermediate/{ID}/background_region/{target_region}_background_regions.bed",
    params:
        length = lambda wildcards, input: get_length(input.region)
    conda:"workflow/envs/background.yml"
    shell:
        """
        bedtools random -n 1000 -l {params.length} -g {input.genome} | \
        bedtools shuffle -i stdin -g {input.genome} -excl {input.gap} -noOverlapping \
        > {output}
        """

rule target_regions_by_chrom:
    input:
        Region_file ="results/intermediate/{ID}/regions/{target_region}_blacklist-excluded.bed",
        blregion = "results/intermediate/{ID}/regions/{target_region}_blacklist-excluded.bed"
    output:
        "results/intermediate/{ID}/simulations/regions/{target_region}/tmp/{target_region}_{chrom}.bed"
    params:
        chrom= lambda wc: wc.chrom
    conda: "workflow/envs/overlays.yml"
    script:
        "workflow/scripts/simulations/split_region_by_chrom.py"

rule simulate_reads:
    input:
        lenDist=config["lenDist"],
        fwdKMerGenome=config["fwdKMerGenome"],
        fwdPKMers=config["fwdPKMers"],
        fwdMKMers=config["fwdMKMers"],
        revPKMers=config["revPKMers"],
        revMKMers=config["revMKMers"],
        regionfile_chrom = "results/intermediate/{ID}/simulations/regions/{target_region}/tmp/{target_region}_{chrom}.bed",
    output: 
        "results/intermediate/{ID}/simulations/tmp/{SAMPLE}.{target_region}_{kmer}mer_sim_{chrom}.bam"
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


rule combine_chromosome_files:
    input:
        files = lambda wildcards: get_chroms(wildcards.ID ,wildcards.target_region),
        #files = lambda wc: get_chroms("results/intermediate/{wc.ID}/regions/{wc.target_region}_blacklist-excluded.bed")
    output:
        "results/intermediate/{ID}/simulations/{SAMPLE}.{target_region}_{kmer}mer_allChrom.bam"
    params:
        chr_sample = "results/intermediate/{ID}/simulations/{SAMPLE}.{target_region}_{kmer}mer_sim_*.bam",
        allChrom = "results/intermediate/{ID}/simulations/{SAMPLE}.{target_region}_{kmer}mer_allChrom"
    conda: "workflow/envs/cfDNA.yml"
    threads: config["threads"]
    shell:
        """
        workflow/scripts/samtools merge -u {input.files} | workflow/scripts/samtools sort - {params.allChrom}
        """
        #"""
        #samtools merge -@ {threads} -u {input.files} | samtools sort -T /fast/users/roeners_c/scratch/tmp/{wildcards.SAMPLE} -O BAM -o {output} -@ {threads} - 
        #"""
        #"""
        #workflow/scripts/samtools merge -u {params.chr_sample} | workflow/scripts/samtools sort - {params.allChrom}
        #"""

rule index_and_stats:
    input:
        "results/intermediate/{ID}/simulations/{SAMPLE}.{target_region}_{kmer}mer_allChrom.bam"
    output:
        "results/intermediate/{ID}/simulations/{SAMPLE}.{target_region}_{kmer}mer_allChrom.bam.bai",
        "results/intermediate/{ID}/simulations/{SAMPLE}.{target_region}_{kmer}mer_allChrom.bam_stats",
        "results/intermediate/{ID}/simulations/{SAMPLE}.{target_region}_{kmer}mer_allChrom.byChromCounts.txt"
    params:
        allChrom = "results/intermediate/{ID}/simulations/{SAMPLE}.{target_region}_{kmer}mer_allChrom"
    conda: "workflow/envs/read_preprocessing.yml"#"workflow/envs/cfDNA.yml"
    threads: config["threads"]
    shell:
        """
        samtools index -@ {threads} {input}
        samtools flagstat -@ {threads} {input} > {params.allChrom}.bam_stats
        samtools view -@ {threads} -F _2 {input} | cut -f 3 | uniq -c > {params.allChrom}.byChromCounts.txt
        """



rule target_regions_by_chrom_background:
    input:
        Region_file ="results/intermediate/{ID}/background_region/{target_region}_background_regions.bed",
        blregion = "results/intermediate/{ID}/background_region/{target_region}_background_regions.bed"
    output:
        "results/intermediate/{ID}/simulations/background_region/regions/{target_region}/tmp/{target_region}_{chrom}.bed"
    params:
        chrom= lambda wc: wc.chrom
    conda: "workflow/envs/overlays.yml"
    script:
        "workflow/scripts/simulations/split_region_by_chrom.py"

rule simulate_reads_background:
    input:
        lenDist=config["lenDist"],
        fwdKMerGenome=config["fwdKMerGenome"],
        fwdPKMers=config["fwdPKMers"],
        fwdMKMers=config["fwdMKMers"],
        revPKMers=config["revPKMers"],
        revMKMers=config["revMKMers"],
        regionfile_chrom = "results/intermediate/{ID}/simulations/background_region/regions/{target_region}/tmp/{target_region}_{chrom}.bed",
    output: 
        "results/intermediate/{ID}/simulations/tmp/background_region/background_{SAMPLE}.{target_region}_{kmer}mer_sim_{chrom}.bam"
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


rule combine_chromosome_files_background:
    input:
        files = lambda wildcards: get_chroms_back(wildcards.ID ,wildcards.target_region),
        #files = lambda wc: get_chroms("results/intermediate/{wc.ID}/regions/{wc.target_region}_blacklist-excluded.bed")
    output:
        "results/intermediate/{ID}/simulations/background_region/{SAMPLE}.{target_region}_{kmer}mer_allChrom.bam"
    params:
        chr_sample = "results/intermediate/{ID}/simulations/background_region/{SAMPLE}.{target_region}_{kmer}mer_sim_*.bam",
        allChrom = "results/intermediate/{ID}/simulations/background_region/{SAMPLE}.{target_region}_{kmer}mer_allChrom"
    conda: "workflow/envs/cfDNA.yml"
    threads: config["threads"]
    shell:
        """
        workflow/scripts/samtools merge -u {input.files} | workflow/scripts/samtools sort - {params.allChrom}
        """
        #"""
        #samtools merge -@ {threads} -u {input.files} | samtools sort -T /fast/users/roeners_c/scratch/tmp/{wildcards.SAMPLE} -O BAM -o {output} -@ {threads} - 
        #"""
        #"""
        #workflow/scripts/samtools merge -u {params.chr_sample} | workflow/scripts/samtools sort - {params.allChrom}
        #"""

rule index_and_stats_background:
    input:
        "results/intermediate/{ID}/simulations/background_region/{SAMPLE}.{target_region}_{kmer}mer_allChrom.bam"
    output:
        "results/intermediate/{ID}/simulations/background_region/{SAMPLE}.{target_region}_{kmer}mer_allChrom.bam.bai",
        "results/intermediate/{ID}/simulations/background_region/{SAMPLE}.{target_region}_{kmer}mer_allChrom.bam_stats",
        "results/intermediate/{ID}/simulations/background_region/{SAMPLE}.{target_region}_{kmer}mer_allChrom.byChromCounts.txt"
    params:
        allChrom = "results/intermediate/{ID}/simulations/background_region/{SAMPLE}.{target_region}_{kmer}mer_allChrom"
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
        target = "results/intermediate/{ID}/regions/{target_region}_blacklist-excluded.bed",
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
        background="results/intermediate/{ID}/background_region/{target_region}_background_regions.bed",
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

rule extract_counts_sim:
    input:
        target = "results/intermediate/{ID}/regions/{target_region}_blacklist-excluded.bed",
        BAMFILE="results/intermediate/{ID}/simulations/{SAMPLE}.{target_region}_{kmer}mer_allChrom.bam",
    output:
        WPS="results/intermediate/{ID}/table/sim_{kmer}mer_{target_region}.{SAMPLE}_WPS.csv",
        COV="results/intermediate/{ID}/table/sim_{kmer}mer_{target_region}.{SAMPLE}_COV.csv",
        STARTS="results/intermediate/{ID}/table/sim_{kmer}mer_{target_region}.{SAMPLE}_STARTS.csv",
    params:
        minRL=config["minRL"],
        maxRL=config["maxRL"],
        out_pre="results/intermediate/{ID}/table/sim_{kmer}mer_{target_region}.{SAMPLE}_%s.csv",
    conda:"workflow/envs/cfDNA.yml"
    shell:
        """
        workflow/scripts/WPS/extractFromBAM_RegionBed_WPS_Cov.py \
        --minInsert={params.minRL} \
        --maxInsert={params.maxRL} \
        -i {input.target} \
        -o {params.out_pre} {input.BAMFILE}
        """

rule extract_counts_back_sim:
    input:
        target = "results/intermediate/{ID}/regions/{target_region}_blacklist-excluded.bed",
        BAMFILE="results/intermediate/{ID}/simulations/background_region/{SAMPLE}.{target_region}_{kmer}mer_allChrom.bam",
    output:
        WPS="results/intermediate/{ID}/table/background_region/sim_{kmer}mer_{target_region}.{SAMPLE}_WPS.csv",
        COV="results/intermediate/{ID}/table/background_region/sim_{kmer}mer_{target_region}.{SAMPLE}_COV.csv",
        STARTS="results/intermediate/{ID}/table/background_region/sim_{kmer}mer_{target_region}.{SAMPLE}_STARTS.csv",
    params:
        minRL=config["minRL"],
        maxRL=config["maxRL"],
        out_pre="results/intermediate/{ID}/table/sim_{kmer}mer_{target_region}.{SAMPLE}_%s.csv",
    conda:"workflow/envs/cfDNA.yml"
    shell:
        """
        workflow/scripts/WPS/extractFromBAM_RegionBed_WPS_Cov.py \
        --minInsert={params.minRL} \
        --maxInsert={params.maxRL} \
        -i {input.target} \
        -o {params.out_pre} {input.BAMFILE}
        """

rule plot_overlays:
    input:
        WPS="results/intermediate/{ID}/table/{target_region}.{SAMPLE}_WPS.csv",
        WPS_ref=lambda wildcards: get_WPS_ref(wildcards.SAMPLE),
        COV="results/intermediate/{ID}/table/{target_region}.{SAMPLE}_COV.csv",
        COV_ref=lambda wildcards: get_COV_ref(wildcards.SAMPLE),
        WPS_back="results/intermediate/{ID}/background_region/table/sim_{kmer}mer_{target_region}.{SAMPLE}_WPS.csv",
        WPS_back_ref=lambda wildcards: get_WPS_background_ref(wildcards.SAMPLE),
        COV_back="results/intermediate/{ID}/background_region/table/sim_{kmer}mer_{target_region}.{SAMPLE}_COV.csv",
        COV_back_ref=lambda wildcards: get_COV_background_ref(wildcards.SAMPLE),
    output:
        "results/plots/overlays/{ID}/{target_region}.{SAMPLE}_overlays_sim_{kmer}mer.pdf",
    params:
        target="{target_region}",
        sample="{SAMPLE}",
        ref_IDs=lambda wildcards: samples["ref_samples"][wildcards.SAMPLE].split(","),
    conda:
        "workflow/envs/overlays.yml"
    script:
        "workflow/scripts/WPS/overlays.py"