from snakemake.utils import validate
import pandas as pd

#import csv

configfile: "config/config_simulations.yml"

#validate(config, schema="workflow/schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
samples.index.names = ["sample_id"]
#validate(samples, schema="workflow/schemas/samples.schema.yaml")
regions = pd.read_csv(config["regions"], sep="\t").set_index("target", drop=False)
regions.index.names = ["region_id"]
#validate(regions, schema="workflow/schemas/regions.schema.yaml")

def get_chroms(target_region):
    BEDCOLS = ["chrom", "chromStart", "chromEnd",
           "name", "score", "strand",
           "thickStart", "thickEnd", "itemRGB",
           "blockCount", "blockSizes", "blockStarts"]
    bedfile=regions["path"][target_region]
    bed_df = pd.read_csv(bedfile, sep="\t", header=None)
    colnames = []
    for i in range(len(bed_df.columns)):
        colnames.append(BEDCOLS[i])
    bed_df.columns = colnames
    chroms = bed_df["chrom"].unique().tolist()
    f_list = expand("results/intermediate/{{ID}}/simulations/{{SAMPLE}}.{{target_region}}_{{kmer}}mer_sim_{chrom}.bam",#"results/intermediate/{{ID}}/simulations/regions/{{target_region}}/{{target_region}}_{chrom}.bed",
                    chrom=chroms)
    return f_list


rule all:
    input:
        expand("results/intermediate/{ID}/simulations/refKmer_grch37_regChroms_{kmer}mers.tsv",
                kmer=config["simulations"]["kmer"],
                ID=samples["ID"],
                ),
        expand("results/intermediate/{ID}/simulations/{SAMPLE}_lenDist.tsv",
                SAMPLE=samples["sample"],
                ID=samples["ID"],
                ),
        expand("results/intermediate/{ID}/simulations/{SAMPLE}_left_{kmer}mer_f.tsv", # remove after implementing chrom combine
                SAMPLE=samples["sample"],
                ID=samples["ID"],
                kmer=config["simulations"]["kmer"],
                ),
        expand("results/intermediate/{ID}/simulations/{SAMPLE}.{target_region}_{kmer}mer_allChrom.bam",
                SAMPLE=samples["sample"],
                ID=samples["ID"],
                target_region=regions["target"],
                kmer=config["simulations"]["kmer"],
                ),
        expand("results/intermediate/{ID}/simulations/{SAMPLE}.{target_region}_{kmer}mer_allChrom.bam.bai",
                SAMPLE=samples["sample"],
                ID=samples["ID"],
                target_region=regions["target"],
                kmer=config["simulations"]["kmer"],
                ),
        expand("results/intermediate/{ID}/simulations/{SAMPLE}.{target_region}_{kmer}mer_allChrom.bam_stats",
                SAMPLE=samples["sample"],
                ID=samples["ID"],
                target_region=regions["target"],
                kmer=config["simulations"]["kmer"],
                ),
        expand("results/intermediate/{ID}/simulations/{SAMPLE}.{target_region}_{kmer}mer_allChrom.byChromCounts.txt",
                SAMPLE=samples["sample"],
                ID=samples["ID"],
                target_region=regions["target"],
                kmer=config["simulations"]["kmer"],
                ),
        expand("results/intermediate/{ID}/simulations/{SAMPLE}.{target_region}_{kmer}mer_allChrom.fragpatterns.txt",
                SAMPLE=samples["sample"],
                ID=samples["ID"],
                target_region=regions["target"],
                kmer=config["simulations"]["kmer"],
                ),
        expand("results/intermediate/{ID}/simulations/{SAMPLE}.{target_region}_{kmer}mer_allChrom.length.tsv",
                SAMPLE=samples["sample"],
                ID=samples["ID"],
                target_region=regions["target"],
                kmer=config["simulations"]["kmer"],
                ),


rule referenceKMers:
    input:
        GRCh37=config["GRCh37"]
    output:
        "results/intermediate/{ID}/simulations/refKmer_grch37_regChroms_{kmer}mers.tsv"
    conda: "workflow/envs/cfDNA.yml"
    shell:
        "workflow/scripts/simulations/referenceKMers.py -f {input.GRCh37} -k {wildcards.kmer} -r '' -o {output}"


rule Extact_length_distribution:
    input:
        BAMFILE=lambda wildcards: samples["path"][wildcards.SAMPLE]
    output:
        "results/intermediate/{ID}/simulations/{SAMPLE}_lenDist.tsv"
    conda: "workflow/envs/cfDNA.yml"
    shell:
        """
        workflow/scripts/simulations/BAM_RG_Length.py --noRG -p tmp_outfile {input.BAMFILE}
        head -n 501 tmp_outfile.tsv > {output}
        """


rule BAM2FragKMers:
    input:
        BAMFILE=lambda wildcards: samples["path"][wildcards.SAMPLE]
    output:
        "results/intermediate/{ID}/simulations/{SAMPLE}_left_{kmer}mer_f.tsv",
        "results/intermediate/{ID}/simulations/{SAMPLE}_left_{kmer}mer_r.tsv",
        "results/intermediate/{ID}/simulations/{SAMPLE}_right_{kmer}mer_f.tsv",
        "results/intermediate/{ID}/simulations/{SAMPLE}_right_{kmer}mer_r.tsv"
    params:
        outfileLE="results/intermediate/{ID}/simulations/{SAMPLE}_left_{kmer}mer",
        outfileRE="results/intermediate/{ID}/simulations/{SAMPLE}_right_{kmer}mer"
    conda: "workflow/envs/cfDNA.yml"
    shell:
        """
        workflow/scripts/simulations/BAM2FragKMers.py --kmerLE={wildcards.kmer} --kmerRE={wildcards.kmer} \
        -v -r 'ALL' {input.BAMFILE} \
        --outfileLE={params.outfileLE} --outfileRE={params.outfileRE}"""

rule target_regions_by_chrom:
    input:
        Region_file = lambda wildcards: regions["path"][wildcards.target_region]
    output:
        "results/intermediate/{ID}/simulations/regions/{target_region}/{target_region}_{chrom}.bed"
    params:
        chrom= lambda wc: wc.chrom
    conda: "workflow/envs/overlays.yml"
    script:
        "workflow/scripts/simulations/split_region_by_chrom.py"

rule simulate_reads:
    input:
        lenDist="results/intermediate/{ID}/simulations/{SAMPLE}_lenDist.tsv",
        fwdKMerGenome="results/intermediate/{ID}/simulations/refKmer_grch37_regChroms_{kmer}mers.tsv",
        fwdPKMers="results/intermediate/{ID}/simulations/{SAMPLE}_left_{kmer}mer_f.tsv",
        fwdMKMers="results/intermediate/{ID}/simulations/{SAMPLE}_left_{kmer}mer_r.tsv",
        revPKMers="results/intermediate/{ID}/simulations/{SAMPLE}_right_{kmer}mer_f.tsv",
        revMKMers="results/intermediate/{ID}/simulations/{SAMPLE}_right_{kmer}mer_r.tsv",
        regionfile_chrom = "results/intermediate/{ID}/simulations/regions/{target_region}/{target_region}_{chrom}.bed",
    output: 
        "results/intermediate/{ID}/simulations/{SAMPLE}.{target_region}_{kmer}mer_sim_{chrom}.bam"
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
        files = lambda wildcards: get_chroms(wildcards.target_region),
    output:
        "results/intermediate/{ID}/simulations/{SAMPLE}.{target_region}_{kmer}mer_allChrom.bam"
    params:
        chr_sample = "results/intermediate/{ID}/simulations/{SAMPLE}.{target_region}_{kmer}mer_sim_*.bam",
        allChrom = "results/intermediate/{ID}/simulations/{SAMPLE}.{target_region}_{kmer}mer_allChrom"
    conda: "workflow/envs/cfDNA.yml"
    shell:
        """
        workflow/scripts/samtools merge -u {params.chr_sample} | workflow/scripts/samtools sort - {params.allChrom}
        """

rule index_and_stats:
    input:
        "results/intermediate/{ID}/simulations/{SAMPLE}.{target_region}_{kmer}mer_allChrom.bam"
    output:
        "results/intermediate/{ID}/simulations/{SAMPLE}.{target_region}_{kmer}mer_allChrom.bam.bai",
        "results/intermediate/{ID}/simulations/{SAMPLE}.{target_region}_{kmer}mer_allChrom.bam_stats",
        "results/intermediate/{ID}/simulations/{SAMPLE}.{target_region}_{kmer}mer_allChrom.byChromCounts.txt"
    params:
        allChrom = "results/intermediate/{ID}/simulations/{SAMPLE}.{target_region}_{kmer}mer_allChrom"
    conda: "workflow/envs/cfDNA.yml"
    shell:
        """
        workflow/scripts/samtools index {input}
        workflow/scripts/samtools flagstat {input} > {params.allChrom}.bam_stats
        workflow/scripts/samtools view -F _2 {input} | cut -f 3 | uniq -c > {params.allChrom}.byChromCounts.txt
        """

rule BAM2FragmentationPatterns_sim:
    input:
        "results/intermediate/{ID}/simulations/{SAMPLE}.{target_region}_{kmer}mer_allChrom.bam"
    output:
        "results/intermediate/{ID}/simulations/{SAMPLE}.{target_region}_{kmer}mer_allChrom.fragpatterns.txt"
    params:
        fasta = config["GRCh37"],
        allChrom = "results/intermediate/{ID}/simulations/{SAMPLE}.{target_region}_{kmer}mer_allChrom"
    conda: "workflow/envs/cfDNA.yml"
    shell:
        """
        workflow/scripts/simulations/BAM2FragmentationPatterns.py -f {params.fasta} -r ALL -o {params.allChrom}.fragpatterns.txt {input}
        """

rule Bam_RG_Length_sim:
    input:
        "results/intermediate/{ID}/simulations/{SAMPLE}.{target_region}_{kmer}mer_allChrom.bam"
    output:
        "results/intermediate/{ID}/simulations/{SAMPLE}.{target_region}_{kmer}mer_allChrom.length.tsv"
    params:
        allChrom = "results/intermediate/{ID}/simulations/{SAMPLE}.{target_region}_{kmer}mer_allChrom"
    conda: "workflow/envs/cfDNA.yml"
    shell:
        """
        workflow/scripts/simulations/BAM_RG_Length.py --noRG -p {params.allChrom}.length {input}
        """
#        