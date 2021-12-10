# cfDNA-Workflow Changelog

## develop

- general changes:
    - samples.schema.yaml
        - added optional info field for misc. info
        - added optional status field for sample status (e.g. disease status)


- snakefile_WPS:
    - added dynamic .gz support
    - most intermediate files are .gz compressed



## v0.3.2

- Hotfix: fixed bug in get_refsample

- update .gitignore

## v0.3.1

- Hotfix: unexpected behavior in conda environments

## v0.3.0

- updated documentation

- comfort:
    - added geneIDs to WPS/COV/STARTs and FFT tables

- snakefile_GE_analysis:
    - added background normalization
    - ported FFT calculation to python
    - separated multiple steps for increased maintainability

- snakefile_WPS:
    - restructured file paths
    - moved window extension to separate
    - moved WPS padding correction to extractFromBAM_RegionBed_WPS_Cov.py

- GE_unsupervised:
    - added rules for generating clusterings and heatmaps between all samples of an experiment using FFT intensities
        - kmeans clustering
        - HDBSCAN clustering
        - heatmaps
        - clustermaps (heatmap with hierarchical clustering)

## v0.2.1

- bugfix paths
- simplified path handling in plotting scripts
- improved readability of label files

## v0.2.0

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
