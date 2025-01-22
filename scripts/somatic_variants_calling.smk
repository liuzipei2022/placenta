## 2024-07-18
## only tumor mode for somatic variants calling
## require env vep109

Files = glob_wildcards("2.mapping/{sample}_sorted_markdup.bam")

rule all:
    input:
        expand("5.somatic/{sample}.vcf.gz", sample=Files.sample),
        expand("5.somatic/filter/{sample}_pileups.table", sample=Files.sample),
        expand("5.somatic/filter/{sample}_contamination.table",sample=Files.sample),
        expand("5.somatic/filter/{sample}_segments.tsv",sample=Files.sample),
        expand("5.somatic/filter/{sample}_filtered.vcf.gz",sample=Files.sample),
        expand("5.somatic/filter/{sample}_filtered.vcf",sample=Files.sample),
        expand("5.somatic/filter/{sample}.vep.maf",sample=Files.sample),
        expand("5.somatic/filter/{sample}_filtered.exons.vcf",sample=Files.sample),
        expand("5.somatic/filter/{sample}.vep.exons.maf",sample=Files.sample),


rule Mutect2:
    input:
        "2.mapping/{sample}_sorted_markdup.bam",
    output:
        "5.somatic/{sample}.vcf.gz"
    shell:
        """
        /picb/lilab/tools/gatk-4.2.6.1/gatk Mutect2 \
        -R /picb/lilab5/liuzipei/placenta/WES_1/gatk_ref/Homo_sapiens_assembly38.fasta \
        -I {input[0]} \
        --germline-resource /picb/lilab5/liuzipei/placenta/WES_1/pon/af-only-gnomad.hg38.vcf.gz \
        -O {output[0]}
        """

rule GetPileupSummaries:
    input:
        "2.mapping/{sample}_sorted_markdup.bam",
    output:
        "5.somatic/filter/{sample}_pileups.table"
    shell:
        """
        /picb/lilab/tools/gatk-4.2.6.1/gatk GetPileupSummaries \
        -I {input[0]} \
        -V /picb/lilab5/liuzipei/placenta/WES_1/somatic_pair/filter/small_exac_common_3.hg38.vcf.gz \
        -L /picb/lilab5/liuzipei/placenta/WES_1/somatic_pair/filter/small_exac_common_3.hg38.vcf.gz \
        -O {output[0]}
        """

rule CalculateContamination:
    input:
        "5.somatic/filter/{sample}_pileups.table",
    output:
        "5.somatic/filter/{sample}_contamination.table",
        "5.somatic/filter/{sample}_segments.tsv"
    shell:
        """
        /picb/lilab/tools/gatk-4.2.6.1/gatk CalculateContamination \
        -I {input[0]} \
        -segments {output[1]} \
        -O {output[0]}
        """

rule FilterMutectCalls:
    input:
        "5.somatic/{sample}.vcf.gz",
        "5.somatic/filter/{sample}_contamination.table",
        "5.somatic/filter/{sample}_segments.tsv",
    output:
        "5.somatic/filter/{sample}_filtered.vcf.gz"
    shell:
        """
        /picb/lilab/tools/gatk-4.2.6.1/gatk FilterMutectCalls \
        -R /picb/lilab5/liuzipei/placenta/WES_1/gatk_ref/Homo_sapiens_assembly38.fasta \
        -V {input[0]} \
        --contamination-table {input[1]} \
        --tumor-segmentation {input[2]} \
        -O {output[0]}
        """

rule gunzip:
    input:
        "5.somatic/filter/{sample}_filtered.vcf.gz"
    output:
        "5.somatic/filter/{sample}_filtered.vcf"
    shell:
        """
        gunzip -c {input[0]} > {output[0]}
        """

rule vcf2maf:
    input:
        "5.somatic/filter/{sample}_filtered.vcf",
    output:
        "5.somatic/filter/{sample}.vep.maf"
    shell:
        """
        perl /picb/lilab5/liuzipei/biosoft/mskcc-vcf2maf-f6d0c40/vcf2maf.pl \
        --input-vcf {input[0]} \
        --output-maf  {output[0]} \
        --vep-path /picb/lilab5/liuzipei/biosoft/miniconda3.bck1/envs/vep109/bin/  \
        --vep-data /picb/lilab5/annotation/VEP/cache_109 \
        --ref-fasta /picb/lilab5/liuzipei/placenta/WES_1/gatk_ref/Homo_sapiens_assembly38.fasta \
        --ncbi-build GRCh38 \
        --tumor-id {wildcards.sample}
        """


rule keep_only_exons_regions:
    input:
        "5.somatic/filter/{sample}_filtered.vcf.gz"
    output:
        "5.somatic/filter/{sample}_filtered.exons.vcf"
    shell:
        """
        bcftools view -f PASS -R /picb/lilab5/liuzipei/placenta/WES_1/5.somatic/filter/hg38.exons.bed \
        {input[0]} > {output[0]}
        """

rule vcf2maf_only_exons:
    input:
        "5.somatic/filter/{sample}_filtered.exons.vcf",
    output:
        "5.somatic/filter/{sample}.vep.exons.maf"
    shell:
        """
        perl /picb/lilab5/liuzipei/biosoft/mskcc-vcf2maf-f6d0c40/vcf2maf.pl \
        --input-vcf {input[0]} \
        --output-maf  {output[0]} \
        --vep-path /picb/lilab5/liuzipei/biosoft/miniconda3.bck1/envs/vep109/bin/  \
        --vep-data /picb/lilab5/annotation/VEP/cache_109 \
        --ref-fasta /picb/lilab5/liuzipei/placenta/WES_1/gatk_ref/Homo_sapiens_assembly38.fasta \
        --ncbi-build GRCh38 \
        --tumor-id {wildcards.sample}
        """
