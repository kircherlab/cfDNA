# cfDNA-Workflow Changelog

## develop

### refactor

- changed config.yml to example.config.yml
- updated schemas
- updated Documentation

### new files

- updated blacklists for both genome builds
- updated genome files for both genome builds

### new feature

- snakefile_GE_analysis:
    - updated ProteinAtlas
        - 3 types: Blood, Tissue, Tissue+cell-lines
    - added support for GRCh38

- snakefile_WPS:
    - added support for normalization by random background sequences
        - generate random background sequences not overlapping target regions
        - calculate WPS,COV,STARTS table for background sequences

- extractFromBAM_RegionBed_WPS_Cov.py
    - added output for fragment endpoints (STARTS output)
    - added strand specificity

- overlays.py:
    - integrated normalization by background sequences

## v0.1.2

### refactor

- overlays.py
    - adding samples is now wrapped in a function
    - flanking regions are now determined dynamically

### documentation

- added README.md in config dir

- updated README.md in root dir
    - updated input sections for workflows documentation
    - added links to .tsv files and README.md in config dir