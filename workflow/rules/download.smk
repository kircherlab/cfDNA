

rule all:
    input:
        "resources/candidate_regions2/Meuleman/input/DHS_Index_and_Vocabulary_hg19_WM20190703.txt.gz",
        "resources/candidate_regions2/Meuleman/input/DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz",
        "resources/candidate_regions2/Meuleman/input/DHS_Index_and_Vocabulary_legend.txt",
        "resources/candidate_regions2/Meuleman/input/41586_2020_2559_MOESM4_ESM.xlsx",
        "resources/common/ensemble/Homo_sapiens.GRCh38.102.gtf.noheader.gz",
        "resources/common/ensemble/Homo_sapiens.GRCh37.87.gtf.noheader.gz",
        "resources/common/liftover/hg38ToHg19.over.chain.gz"

rule download_gtf:
    output:
        "resources/common/ensemble/Homo_sapiens.GRCh38.102.gtf.gz",
        "resources/common/ensemble/Homo_sapiens.GRCh37.87.gtf.gz"
    shell:
        """
        wget -P resources/common/ensemble/ ftp://ftp.ensembl.org/pub/release-102/gtf/homo_sapiens/Homo_sapiens.GRCh38.102.gtf.gz
        wget -P resources/common/ensemble/ ftp://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz
        """

rule remove_gtf_header:
    input:
        GRCh37 = "resources/common/ensemble/Homo_sapiens.GRCh37.87.gtf.gz",
        GRCh38 = "resources/common/ensemble/Homo_sapiens.GRCh38.102.gtf.gz",
    output:
        GRCh37 = "resources/common/ensemble/Homo_sapiens.GRCh37.87.gtf.noheader.gz",
        GRCh38 = "resources/common/ensemble/Homo_sapiens.GRCh38.102.gtf.noheader.gz",
    conda: "../../workflow/envs/read_preprocessing.yml"
    shell:
        """
        zcat {input.GRCh37} | grep -E -v "^#" | bgzip -c > {output.GRCh37};
        zcat {input.GRCh38} | grep -E -v "^#" | bgzip -c > {output.GRCh38};
        """


rule download_Meuleman_data:
    output:
        "resources/candidate_regions2/Meuleman/input/DHS_Index_and_Vocabulary_hg19_WM20190703.txt.gz",
        "resources/candidate_regions2/Meuleman/input/DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz",
        "resources/candidate_regions2/Meuleman/input/DHS_Index_and_Vocabulary_legend.txt",
        "resources/candidate_regions2/Meuleman/input/41586_2020_2559_MOESM4_ESM.xlsx",
    shell:
        """
        wget -P resources/candidate_regions2/Meuleman/input/ https://zenodo.org/record/3838751/files/DHS_Index_and_Vocabulary_hg19_WM20190703.txt.gz
        wget -P resources/candidate_regions2/Meuleman/input/ https://zenodo.org/record/3838751/files/DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz
        wget -P resources/candidate_regions2/Meuleman/input/ https://zenodo.org/record/3838751/files/DHS_Index_and_Vocabulary_legend.txt
        wget -P resources/candidate_regions2/Meuleman/input/ https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-020-2559-3/MediaObjects/41586_2020_2559_MOESM4_ESM.xlsx
        """

rule download_liftover_chains:
    output:
        "resources/common/liftover/hg38ToHg19.over.chain.gz"
    shell:
        """
        wget --timestamping \
            'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz' \
            -O resources/common/liftover/hg38ToHg19.over.chain.gz
        """