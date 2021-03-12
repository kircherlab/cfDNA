component_list = [
    "tissue_invariant",
    "stromal_a",
    "stromal_b",
    "primitive_embryonic",
    "renal_cancer",
    "lymphoid",
    "pulmonary_devel",
    "cardiac",
    "musculoskeletal",
    "myeloid_erythroid",
    "placental_trophoblast",
    "digestive",
    "organ_devel_renal",
    "neural",
    "vascular_endothelial",
    "cancer_epithelial",
]

configfile: "config/config_components.yml"

rule all:
    input:
        "resources/candidate_regions2/Meuleman/GRCh38/DHS/DHS_regions_all_summit_GRCh38.tsv.gz",
        "resources/candidate_regions2/Meuleman/GRCh37/DHS/DHS_regions_all_summit_GRCh37.tsv",
        "resources/candidate_regions2/Meuleman/GRCh38/DHS/DHS_mean_signal_matrix_GRCh38.tsv.gz",
        "resources/candidate_regions2/Meuleman/GRCh38/DHS/DHS_mean_signal_matrix_binary_GRCh38.tsv.gz",
        "resources/candidate_regions2/Meuleman/GRCh37/DHS/DHS_mean_signal_matrix_GRCh37.tsv.gz",
        "resources/candidate_regions2/Meuleman/GRCh37/DHS/DHS_mean_signal_matrix_binary_GRCh37.tsv.gz",
        expand("resources/candidate_regions2/Meuleman/GRCh38/DHSAnno/components/Top1000_{comp}_DHS_summitAnno_GRCh38.tsv.gz", comp=component_list),
        expand("resources/candidate_regions2/Meuleman/GRCh38/DHS/components/Top1000_{comp}_DHS_summit_GRCh38.bed.gz", comp=component_list),
        "resources/candidate_regions2/Meuleman/GRCh38/DHSAnno/Top1000_all_DHS_summitAnno_GRCh38.tsv.gz",
        "resources/candidate_regions2/Meuleman/GRCh38/DHS/Top1000_all_DHS_summit_GRCh38.tsv.gz",
        expand("resources/candidate_regions2/Meuleman/GRCh37/DHSAnno/components/Top1000_{comp}_DHS_summitAnno_GRCh37.tsv.gz", comp=component_list),
        expand("resources/candidate_regions2/Meuleman/GRCh37/DHS/components/Top1000_{comp}_DHS_summit_GRCh37.bed.gz", comp=component_list),
        "resources/candidate_regions2/Meuleman/GRCh37/DHSAnno/Top1000_all_DHS_summitAnno_GRCh37.tsv.gz",
        "resources/candidate_regions2/Meuleman/GRCh37/DHS/Top1000_all_DHS_summit_GRCh37.tsv.gz",


        #"resources/candidate_regions2/Meuleman/GRCh37/DHSAnno/Top1000_all_DHS_summitAnno_GRCh37.tsv.gz",
        #"resources/candidate_regions2/Meuleman/GRCh38/DHSAnno/Top1000_all_DHS_summitAnno_GRCh38.tsv.gz",
        #expand("resources/candidate_regions2/Meuleman/GRCh37/DHSAnno/components/Top1000_{comp}_DHS_summitAnno_GRCh37.tsv.gz", comp=component_list),
        #expand("resources/candidate_regions2/Meuleman/GRCh38/DHSAnno/components/Top1000_{comp}_DHS_summitAnno_GRCh38.tsv.gz", comp=component_list),
        
        expand("resources/candidate_regions2/Meuleman/GRCh37/transcriptAnno/components/Top1000_{comp}_transcriptAnno_GRCh37.tsv.gz", comp=component_list),
        expand("resources/candidate_regions2/Meuleman/GRCh37/TSS/components/Top1000_{comp}_TSS_GRCh37.bed.gz", comp=component_list),
        "resources/candidate_regions2/Meuleman/GRCh37/transcriptAnno/Top1000_all_transcriptAnno_GRCh37.tsv.gz",
        "resources/candidate_regions2/Meuleman/GRCh37/TSS/Top1000_all_TSS_GRCh37.bed.gz",
        expand("resources/candidate_regions2/Meuleman/GRCh38/transcriptAnno/components/Top1000_{comp}_transcriptAnno_GRCh38.tsv.gz", comp=component_list),
        expand("resources/candidate_regions2/Meuleman/GRCh38/TSS/components/Top1000_{comp}_TSS_GRCh38.bed.gz", comp=component_list),
        "resources/candidate_regions2/Meuleman/GRCh38/transcriptAnno/Top1000_all_transcriptAnno_GRCh38.tsv.gz",
        "resources/candidate_regions2/Meuleman/GRCh38/TSS/Top1000_all_TSS_GRCh38.bed.gz",


rule extract_DHS_summit:
    input:
        DHS_Index="resources/candidate_regions2/Meuleman/input/DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz",
    output:
        DHS_summit="resources/candidate_regions2/Meuleman/GRCh38/DHS/DHS_regions_all_summit_GRCh38.tsv.gz",
        DHS_summit_header="resources/candidate_regions2/Meuleman/GRCh38/DHS/DHS_regions_all_summit_GRCh38.tsv.gz.header",
    conda:
        "../../workflow/envs/preparations.yml"
    script:
        "../../workflow/scripts/Meuleman_prep/extract_dhs_summit.py"


rule liftover_DHS_summit_to_GRCh37:
    input:
        regions="resources/candidate_regions2/Meuleman/GRCh38/DHS/DHS_regions_all_summit_GRCh38.tsv.gz",
        header="resources/candidate_regions2/Meuleman/GRCh38/DHS/DHS_regions_all_summit_GRCh38.tsv.gz.header",
        chain="resources/common/liftover/hg38ToHg19.over.chain.gz",
    output:
        regions="resources/candidate_regions2/Meuleman/GRCh37/DHS/DHS_regions_all_summit_GRCh37.tsv",
        header="resources/candidate_regions2/Meuleman/GRCh37/DHS/DHS_regions_all_summit_GRCh37.tsv.header",
        unmapped="resources/candidate_regions2/Meuleman/GRCh37/DHS/DHS_regions_all_summit_GRCh37.unmapped",
    conda:
        "../../workflow/envs/preparations.yml"
    shell:
        """
        cp {input.header} {output.header}; \
        liftOver -bedPlus=3 -tab {input.regions} {input.chain} {output.regions} {output.unmapped}
        """


rule create_DHS_matrices_GRCh38:
    input:
        DHS_summit="resources/candidate_regions2/Meuleman/GRCh38/DHS/DHS_regions_all_summit_GRCh38.tsv.gz",
        header="resources/candidate_regions2/Meuleman/GRCh38/DHS/DHS_regions_all_summit_GRCh38.tsv.gz.header",
    output:
        DHS_continuous="resources/candidate_regions2/Meuleman/GRCh38/DHS/DHS_mean_signal_matrix_GRCh38.tsv.gz",
        DHS_binary="resources/candidate_regions2/Meuleman/GRCh38/DHS/DHS_mean_signal_matrix_binary_GRCh38.tsv.gz",
    conda:
        "../../workflow/envs/preparations.yml"
    script:
        "../../workflow/scripts/Meuleman_prep/create_DHS_matrices.py"


rule create_DHS_matrices_GRCh37:
    input:
        DHS_summit="resources/candidate_regions2/Meuleman/GRCh37/DHS/DHS_regions_all_summit_GRCh37.tsv",
        header="resources/candidate_regions2/Meuleman/GRCh37/DHS/DHS_regions_all_summit_GRCh37.tsv.header",
    output:
        DHS_continuous="resources/candidate_regions2/Meuleman/GRCh37/DHS/DHS_mean_signal_matrix_GRCh37.tsv.gz",
        DHS_binary="resources/candidate_regions2/Meuleman/GRCh37/DHS/DHS_mean_signal_matrix_binary_GRCh37.tsv.gz",
    conda:
        "../../workflow/envs/preparations.yml"
    script:
        "../../workflow/scripts/Meuleman_prep/create_DHS_matrices.py"


rule create_Top1000_TSS_regions_GRCh37:
    input:
        ensemble = "resources/common/ensemble/Homo_sapiens.GRCh37.87.gtf.noheader.gz",
        meuleman = "resources/candidate_regions2/Meuleman/input/41586_2020_2559_MOESM4_ESM.xlsx"
    output:
        expand("resources/candidate_regions2/Meuleman/GRCh37/transcriptAnno/components/Top1000_{comp}_transcriptAnno_GRCh37.tsv.gz", comp=component_list),
        expand("resources/candidate_regions2/Meuleman/GRCh37/TSS/components/Top1000_{comp}_TSS_GRCh37.bed.gz", comp=component_list),
        transcript_all = "resources/candidate_regions2/Meuleman/GRCh37/transcriptAnno/Top1000_all_transcriptAnno_GRCh37.tsv.gz",
        TSS_all = "resources/candidate_regions2/Meuleman/GRCh37/TSS/Top1000_all_TSS_GRCh37.bed.gz",
    params:
        output_transcript_pattern = lambda wildcards: "resources/candidate_regions2/Meuleman/GRCh37/transcriptAnno/components/Top1000_{}_transcriptAnno_GRCh37.tsv.gz",
        output_TSS_pattern = lambda wildcards: "resources/candidate_regions2/Meuleman/GRCh37/TSS/components/Top1000_{}_TSS_GRCh37.bed.gz",
        ref_chroms = [f"chr{i}" for i in config["reference_chromosomes"]]
    conda:"../../workflow/envs/preparations.yml"
    script:
        "../../workflow/scripts/Meuleman_prep/create_TSS_regions.py"

rule create_Top1000_TSS_regions_GRCh38:
    input:
        ensemble = "resources/common/ensemble/Homo_sapiens.GRCh38.102.gtf.noheader.gz",
        meuleman = "resources/candidate_regions2/Meuleman/input/41586_2020_2559_MOESM4_ESM.xlsx"
    output:
        expand("resources/candidate_regions2/Meuleman/GRCh38/transcriptAnno/components/Top1000_{comp}_transcriptAnno_GRCh38.tsv.gz", comp=component_list),
        expand("resources/candidate_regions2/Meuleman/GRCh38/TSS/components/Top1000_{comp}_TSS_GRCh38.bed.gz", comp=component_list),
        transcript_all = "resources/candidate_regions2/Meuleman/GRCh38/transcriptAnno/Top1000_all_transcriptAnno_GRCh38.tsv.gz",
        TSS_all = "resources/candidate_regions2/Meuleman/GRCh38/TSS/Top1000_all_TSS_GRCh38.bed.gz",
    params:
        output_transcript_pattern = lambda wildcards: "resources/candidate_regions2/Meuleman/GRCh38/transcriptAnno/components/Top1000_{}_transcriptAnno_GRCh38.tsv.gz",
        output_TSS_pattern = lambda wildcards: "resources/candidate_regions2/Meuleman/GRCh38/TSS/components/Top1000_{}_TSS_GRCh38.bed.gz",
        ref_chroms = [f"chr{i}" for i in config["reference_chromosomes"]]
    conda:"../../workflow/envs/preparations.yml"
    script:
        "../../workflow/scripts/Meuleman_prep/create_TSS_regions.py"


rule create_Top1000_DHS_regions_GRCh37:
    input:
        DHS_summit="resources/candidate_regions2/Meuleman/GRCh37/DHS/DHS_regions_all_summit_GRCh37.tsv",
        header="resources/candidate_regions2/Meuleman/GRCh37/DHS/DHS_regions_all_summit_GRCh37.tsv.header",
    output:
        expand("resources/candidate_regions2/Meuleman/GRCh37/DHSAnno/components/Top1000_{comp}_DHS_summitAnno_GRCh37.tsv.gz", comp=component_list),
        expand("resources/candidate_regions2/Meuleman/GRCh37/DHS/components/Top1000_{comp}_DHS_summit_GRCh37.bed.gz", comp=component_list),
        summitAnno_all="resources/candidate_regions2/Meuleman/GRCh37/DHSAnno/Top1000_all_DHS_summitAnno_GRCh37.tsv.gz",
        DHS_all="resources/candidate_regions2/Meuleman/GRCh37/DHS/Top1000_all_DHS_summit_GRCh37.tsv.gz",
    params:
        summitAnno_pattern = lambda wildcards: "resources/candidate_regions2/Meuleman/GRCh37/DHSAnno/components/Top1000_{}_DHS_summitAnno_GRCh37.tsv.gz",
        DHS_pattern = lambda wildcards: "resources/candidate_regions2/Meuleman/GRCh37/DHS/components/Top1000_{}_DHS_summit_GRCh37.bed.gz"
    conda:"../../workflow/envs/preparations.yml"
    script:
        "../../workflow/scripts/Meuleman_prep/create_DHS_regions.py"


rule create_Top1000_DHS_regions_GRCh38:
    input:
        DHS_summit="resources/candidate_regions2/Meuleman/GRCh38/DHS/DHS_regions_all_summit_GRCh38.tsv.gz",
        header="resources/candidate_regions2/Meuleman/GRCh38/DHS/DHS_regions_all_summit_GRCh38.tsv.gz.header",
    output:
        expand("resources/candidate_regions2/Meuleman/GRCh38/DHSAnno/components/Top1000_{comp}_DHS_summitAnno_GRCh38.tsv.gz", comp=component_list),
        expand("resources/candidate_regions2/Meuleman/GRCh38/DHS/components/Top1000_{comp}_DHS_summit_GRCh38.bed.gz", comp=component_list),
        summitAnno_all="resources/candidate_regions2/Meuleman/GRCh38/DHSAnno/Top1000_all_DHS_summitAnno_GRCh38.tsv.gz",
        DHS_all="resources/candidate_regions2/Meuleman/GRCh38/DHS/Top1000_all_DHS_summit_GRCh38.tsv.gz",
    params:
        summitAnno_pattern = lambda wildcards: "resources/candidate_regions2/Meuleman/GRCh38/DHSAnno/components/Top1000_{}_DHS_summitAnno_GRCh38.tsv.gz",
        DHS_pattern = lambda wildcards: "resources/candidate_regions2/Meuleman/GRCh38/DHS/components/Top1000_{}_DHS_summit_GRCh38.bed.gz"
    conda:"../../workflow/envs/preparations.yml"
    script:
        "../../workflow/scripts/Meuleman_prep/create_DHS_regions.py"
