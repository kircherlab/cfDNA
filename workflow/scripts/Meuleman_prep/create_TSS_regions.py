import pandas as pd

input_ensemble = snakemake.input["ensemble"]
input_meuleman = snakemake.input["meuleman"]
output_transcript_all = snakemake.output["transcript_all"]
output_transcript_pattern = snakemake.params["output_transcript_pattern"]
output_TSS_all = snakemake.output["TSS_all"]
output_TSS_pattern = snakemake.params["output_TSS_pattern"]
chroms = snakemake.params["ref_chroms"]

#header = pd.read_csv(input_header).columns.to_list()
#df = pd.read_csv(input_file, sep="\t", names=header)

gene_cols = ["chrom", "source", "feature", "start", "end", "score", "strand", "frame", "gene_id", "gene_version", "gene_name", "gene_source", "gene_biotype"]
transcript_cols = ["chrom", "source", "feature", "start", "end", "score", "strand", "frame", 
             "gene_id", "gene_version", "transcript_id", "transcript_version","gene_name", "gene_source", "gene_biotype",
            "transcript_name", "transcript_source", "transcript_biotype"
            ]

ensemble_raw = pd.read_csv(input_ensemble, sep="\t", header=None)

### extract by gene ID ###

# filter for gene
ensemble_genes = ensemble_raw.loc[ensemble_raw[2] == "gene"]
# separate combined fields
ensemble_genes_cleaned = pd.concat([ensemble_genes.iloc[:,0:8], ensemble_genes[8].str.split(';', expand=True).iloc[:,0:-1]], axis=1)
# add column names
ensemble_genes_cleaned.columns = gene_cols
# clean fields
ensemble_genes_cleaned["gene_id"] = ensemble_genes_cleaned["gene_id"].str.replace("gene_id","").str.replace('"', '').str.strip()
ensemble_genes_cleaned["gene_version"] = ensemble_genes_cleaned["gene_version"].str.replace("gene_version","").str.replace('"', '').str.strip()
ensemble_genes_cleaned["gene_name"] = ensemble_genes_cleaned["gene_name"].str.replace("gene_name","").str.replace('"', '').str.strip()
ensemble_genes_cleaned["gene_source"] = ensemble_genes_cleaned["gene_source"].str.replace("gene_source","").str.replace('"', '').str.strip()
ensemble_genes_cleaned["gene_biotype"] = ensemble_genes_cleaned["gene_biotype"].str.replace("gene_biotype","").str.replace('"', '').str.strip()
ensemble_genes_cleaned["chrom"] = "chr" + ensemble_genes_cleaned["chrom"].astype(str)
# save to file
#ensemble_genes_cleaned.to_csv("Human.GRCh38.genes.tsv", sep="\t", index=None)

### extract by transcript ID ###

# filter for transcript
ensemble_transcript = ensemble_raw.loc[ensemble_raw[2] == "transcript"]
# separate combined fields
ensemble_transcripts_cleaned = pd.concat([ensemble_transcript.iloc[:,0:8], ensemble_transcript[8].str.split(';', expand=True).iloc[:,:-6]], axis=1).iloc[:,0:18]
# add column names
ensemble_transcripts_cleaned.columns = transcript_cols
# clean fields
ensemble_transcripts_cleaned["chrom"] = "chr" + ensemble_transcripts_cleaned["chrom"].astype(str)
ensemble_transcripts_cleaned["gene_id"] = ensemble_transcripts_cleaned["gene_id"].str.replace("gene_id","").str.replace('"', '').str.strip()
ensemble_transcripts_cleaned["gene_version"] = ensemble_transcripts_cleaned["gene_version"].str.replace("gene_version","").str.replace('"', '').str.strip()
ensemble_transcripts_cleaned["gene_name"] = ensemble_transcripts_cleaned["gene_name"].str.replace("gene_name","").str.replace('"', '').str.strip()
ensemble_transcripts_cleaned["gene_source"] = ensemble_transcripts_cleaned["gene_source"].str.replace("gene_source","").str.replace('"', '').str.strip()
ensemble_transcripts_cleaned["gene_biotype"] = ensemble_transcripts_cleaned["gene_biotype"].str.replace("gene_biotype","").str.replace('"', '').str.strip()
ensemble_transcripts_cleaned["transcript_id"] = ensemble_transcripts_cleaned["transcript_id"].str.replace("transcript_id","").str.replace('"', '').str.strip()
ensemble_transcripts_cleaned["transcript_version"] = ensemble_transcripts_cleaned["transcript_version"].str.replace("transcript_version","").str.replace('"', '').str.strip()
ensemble_transcripts_cleaned["transcript_name"] = ensemble_transcripts_cleaned["transcript_name"].str.replace("transcript_name","").str.replace('"', '').str.strip()
ensemble_transcripts_cleaned["transcript_source"] = ensemble_transcripts_cleaned["transcript_source"].str.replace("transcript_source","").str.replace('"', '').str.strip()
ensemble_transcripts_cleaned["transcript_biotype"] = ensemble_transcripts_cleaned["transcript_biotype"].str.replace("transcript_biotype","").str.replace('"', '').str.strip()
# save to file
#ensemble_transcripts_cleaned.to_csv("Human.GRCh38.transcripts.tsv", sep="\t", index=None)

### intersect with Meuleman data ###

# filter for target chroms
ensemble_transcripts_filtered = ensemble_transcripts_cleaned.loc[ensemble_transcripts_cleaned["chrom"].isin(chroms)]

# read meuleman supplement table 2
meuleman_df = pd.read_excel(input_meuleman, sheet_name=None)

# extract transcript_anno data for Gene expression workflow
Top1000_transcript_all = pd.DataFrame()
for key in meuleman_df.keys():
    name=key.replace(" & ","_").replace(".","").replace(" ","_").lower()
    outfile = output_transcript_pattern.format(name)
    sample=meuleman_df[key]
    sample=sample.loc[sample["category"] == 'protein_coding'].sort_values(by="top_comp_rank")
    candidates = ensemble_transcripts_filtered.merge(sample, how="inner", left_on="gene_name", right_on="name")
    candidates = candidates.sort_values(by="top_comp_rank")
    candidates = candidates[["gene_id","chrom","start_x", "end_x", "strand_x"]].drop_duplicates(subset=['gene_id'])
    candidates["chrom"] = candidates["chrom"].str.replace("chr","")
    Top1000_transcript_all = pd.concat([Top1000_transcript_all,candidates.iloc[:1000]])
    candidates.iloc[:1000].to_csv(outfile, sep="\t", header=None, index=None)

Top1000_transcript_all.to_csv(output_transcript_all, sep="\t", header=None, index=None)

# extract TSS data

Top1000_TSS_all = pd.DataFrame()
for key in meuleman_df.keys():
    name=key.replace(" & ","_").replace(".","").replace(" ","_").lower()
    outfile = output_TSS_pattern.format(name)
    sample=meuleman_df[key]
    sample=sample.loc[sample["category"] == 'protein_coding'].sort_values(by="top_comp_rank")
    candidates = ensemble_transcripts_filtered.merge(sample, how="inner", left_on="gene_name", right_on="name")
    candidates = candidates.sort_values(by="top_comp_rank").drop_duplicates(subset=['gene_id'])
    regions = pd.DataFrame()
    regions["chrom"]=candidates["chrom"] 
    regions["start"] = candidates["start_x"]
    regions["end"] = candidates["start_x"]+1
    regions["name"] = candidates["gene_id"]
    regions["score"] = candidates["score"]
    regions["strand"] = candidates["strand_x"]
    Top1000_TSS_all = pd.concat([Top1000_TSS_all,regions.iloc[:1000]])
    regions.iloc[:1000].to_csv(outfile, sep="\t", header=None, index=None)

Top1000_TSS_all.to_csv(output_TSS_all, sep="\t", header=None, index=None)