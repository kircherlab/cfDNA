# Documentation for config files

## config.yml

The config.yml file configures values that should stay constant between samples.

```yml
samples: "config/samples.tsv" # .tsv file containing sample names and locations
regions: "config/regions.tsv" # .tsv file containing bed files with regions of interest

tissue: ["CACO.2", "MCF7", "PC.3"] # proteinAtlas tissues for generating plots in GE workflow
refSample: "NPH001" # reference sample for rank correlation comparison
minRL: 120 # minimum read length for calculating WPS
maxRL: 180 # maximum read length for calculating WPS
bpProtection: 120 # value for WPS window
```

## samples.tsv

The samples.tsv contains a header with four columns:

```bash
ID	sample	path	ref_samples
experimentID	testsample1	"/path/to/testsample1.bam"	testsample2,testsample3
experimentID	testsample2	"/path/to/testsample2.bam"	testsample1,testsample3
experimentID	testsample3	"/path/to/testsample3.bam"	testsample1,testsample2
```

- **ID** - ID for a certain analysis to create identifiable directories and/or filenames
- **sample** - sample name used to identify files
- **path** - path to input file
- **ref_sample** - Reference sample for some        visualizations/calculations. ref_samples are comma separated, must be in present in the sample column and every sample needs a ref_sample (e.g. itself).

## regions.tsv

The regions.tsv contains a header with two columns:

```text
target  path
gene1   /path/to/gene1.bed
TF1     /pat/to/TF1BS.bed
```

- **target** - describes the targets defined in the correspoding .bed file
- **path** - path to input .bed file containing coordinates of interest (all coordinates should be centered around a specific feature and of same length)

**Note:** .bed has to contain the first 6 fields (chrom, chromStart, chromEnd, name, value, strand), even though name and value are not actively used.
