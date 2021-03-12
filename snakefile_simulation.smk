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


def get_region_list(bedfile):
    BEDCOLS = ["chrom", "chromStart", "chromEnd",
           "name", "score", "strand",
           "thickStart", "thickEnd", "itemRGB",
           "blockCount", "blockSizes", "blockStarts"]
    bed_df = pd.read_csv(bedfile, sep="\t", header=None)
    colnames = []
    for i in range(len(bed_df.columns)):
        colnames.append(BEDCOLS[i])
    bed_df.columns = colnames
    bed_df.head()
    bed_df["chrom"] = bed_df["chrom"].str.lstrip("chr")
    regions = list()
    for index, row in bed_df.iterrows():
        coord = "{chrom}:{start}-{end}".format(chrom=str(row["chrom"]),start=str(row["chromStart"]),end=str(row["chromEnd"]))
        regions.append(coord)
    return regions

def get_reg(target):
    bedfile = regions["path"][target]
    region_list = get_region_list(bedfile)
    return expand("results/intermediate/{{ID}}/simulations/{{SAMPLE}}.{{target_region}}_{{kmer}}mer_sim_{reg}.bam",
            reg=region_list)


print(config["simulations"]["kmer"])

#
#def get_region_dict(fai, chr_num):
#    cdict = dict()
#    try:
#        with open(fai,"r") as file:
#            count=0
#            fai_reader=csv.reader(file, delimiter="\t")
#            while count<chr_num:
#                line = next(fai_reader)
#                cdict[line[0]] = f"{line[0]}:1-{line[1]}"
#                count+=1
#    except FileNotFoundError as fnf_error:
#        print(fnf_error)
#    else:
#        return cdict
#
#print(config["GRCh37"]+".fai")
#print(samples)
#
#region_dict = get_region_dict(config["GRCh37"]+".fai", 24)
#
#print(region_dict)

rule all:
    input:
        expand("results/intermediate/{ID}/simulations/refKmer_grch37_regChroms_{kmer}mers.tsv",
                kmer=config["simulations"]["kmer"],
                ID=samples["ID"],
                target_region=regions["target"]),
        expand("results/intermediate/{ID}/simulations/{SAMPLE}_lenDist.tsv",
                SAMPLE=samples["sample"],
                ID=samples["ID"],
                target_region=regions["target"]),
        expand("results/intermediate/{ID}/simulations/{SAMPLE}_left_{kmer}mer_f.tsv", # remove after implementing chrom combine
                SAMPLE=samples["sample"],
                ID=samples["ID"],
                kmer=config["simulations"]["kmer"],
                target_region=regions["target"]),
        expand("results/tables/{ID}/simulations/{SAMPLE}.{target_region}_{kmer}mer_allChrom.length.tsv",
                SAMPLE=samples["sample"],
                ID=samples["ID"],
                kmer=config["simulations"]["kmer"],
                target_region=regions["target"])


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


rule simulate_reads:
    input:
        lenDist="results/intermediate/{ID}/simulations/{SAMPLE}_lenDist.tsv",
        fwdKMerGenome="results/intermediate/{ID}/simulations/refKmer_grch37_regChroms_{kmer}mers.tsv",
        fwdPKMers="results/intermediate/{ID}/simulations/{SAMPLE}_left_{kmer}mer_f.tsv",
        fwdMKMers="results/intermediate/{ID}/simulations/{SAMPLE}_left_{kmer}mer_r.tsv",
        revPKMers="results/intermediate/{ID}/simulations/{SAMPLE}_right_{kmer}mer_f.tsv",
        revMKMers="results/intermediate/{ID}/simulations/{SAMPLE}_right_{kmer}mer_r.tsv"
    output:
        temp("results/intermediate/{ID}/simulations/{SAMPLE}.{target_region}_{kmer}mer_sim_{reg}.bam")
    params:
        region = lambda wc: wc.reg,
        sampling_repeats = config["simulations"]["sampling_repeats"],
        fasta = config["GRCh37"]
    conda: "workflow/envs/cfDNA_sim.yml"
    shell:
        """
        workflow/scripts/simulations/simulate_reads.py -v -d {input.lenDist} \
        --fwdKMerGenome={input.fwdKMerGenome} \
        --revKMerGenome={input.fwdKMerGenome} \
        --fwdPKMers={input.fwdPKMers} \
        --fwdMKMers={input.fwdMKMers} \
        --revPKMers={input.revPKMers} \
        --revMKMers={input.revMKMers} \
        -f {params.fasta} \
        -s {params.sampling_repeats} -r {params.region} \
        -v \
        -o {output};
        """


rule combine_chromosome_files:
    input:
        files = lambda wildcards: get_reg(wildcards.target_region),
        #expand("results/intermediate/{ID}/simulations/{SAMPLE}.{target_region}_{kmer}mer_sim_{reg}.bam",
        #        SAMPLE=config["simulations"]["samples"],
        #        kmer=config["simulations"]["kmer"],
        #        reg=get_region_list(regions, target_region))
    output:
        "results/tables/{ID}/simulations/{SAMPLE}.{target_region}_{kmer}mer_allChrom.length.tsv"
    params:
        #chr_sample = "results/tables/{ID}/simulations/{SAMPLE}.{target_region}_{kmer}mer_sim_*.bam",
        chr_sample = "results/intermediate/{ID}/simulations/{SAMPLE}.{target_region}_{kmer}mer_sim_*.bam",
        allChrom = "results/tables/{ID}/simulations/{SAMPLE}.{target_region}_{kmer}mer_allChrom"
    conda: "workflow/envs/cfDNA.yml"
    shell:
        """
        workflow/scripts/samtools merge -u {params.chr_sample} | workflow/scripts/samtools sort - {params.allChrom}
        workflow/scripts/samtools index {params.allChrom}.bam
        workflow/scripts/samtools flagstat {params.allChrom}.bam > {params.allChrom}.bam_stats
        workflow/scripts/samtools view -F _2 {params.allChrom}.bam | cut -f 3 | uniq -c > {params.allChrom}.byChromCounts.txt
        workflow/scripts/simulations/BAM2FragmentationPatterns.py -r ALL -o {params.allChrom}.fragpatterns.txt {params.allChrom}.bam
        workflow/scripts/simulations/BAM_RG_Length.py --noRG -p {params.allChrom}.length {params.allChrom}.bam
        """