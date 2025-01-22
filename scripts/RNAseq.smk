# RNAseq standard pipeline
# Author: Liu Zipei
# Data: 2024-04-14

Files = glob_wildcards("0.rawdata/{sample}_R1.fastq.gz")

rule all:
    input:
        "0.rawdata/fastqc/multiqc_report.html",
        "1.cleandata/fastqc/multiqc_report.html",
        "3.featurecount/count.tsv",
        "3.featurecount/tpm_t.csv",
        expand("0.rawdata/fastqc/{sample}_{reads}_fastqc.html", sample=Files.sample, reads=["R1", "R2"]),
        expand("1.cleandata/fastqc/{sample}_{reads}.clean_fastqc.html", sample=Files.sample, reads=["R1", "R2"]),
        expand("2.mapping/{sample}_Aligned.sortedByCoord.out.bam", sample=Files.sample),
        expand("2.mapping/{sample}_Aligned.sortedByCoord.out.bam.bai", sample=Files.sample),
        

rule fastqc_raw:
    input:
        "0.rawdata/{sample}_R1.fastq.gz",
        "0.rawdata/{sample}_R2.fastq.gz",
    output:
        "0.rawdata/fastqc/{sample}_R1_fastqc.zip",
        "0.rawdata/fastqc/{sample}_R2_fastqc.zip",
        "0.rawdata/fastqc/{sample}_R1_fastqc.html",
        "0.rawdata/fastqc/{sample}_R2_fastqc.html",
    threads: 10
    shell:
        "fastqc -t {threads} -o 0.rawdata/fastqc/ {input}"

rule multiqc_raw:
    priority: -999
    output:
        "0.rawdata/fastqc/multiqc_report.html"
    shell:
        "multiqc -o 0.rawdata/fastqc/ 0.rawdata/fastqc/"
    

rule fastp:
    input:
        "0.rawdata/{sample}_R1.fastq.gz",
        "0.rawdata/{sample}_R2.fastq.gz",
    output:
        "1.cleandata/{sample}_R1.clean.fastq.gz",
        "1.cleandata/{sample}_R2.clean.fastq.gz",
        "1.cleandata/{sample}.html",
        "1.cleandata/{sample}.json",
    threads: 10
    shell:
        """
        fastp --detect_adapter_for_pe -w {threads} \
        -i {input[0]} -I {input[1]} -o {output[0]} -O {output[1]} \
        -h {output[2]} -j {output[3]}
        """

rule fastqc_clean:
    input:
        "1.cleandata/{sample}_R1.clean.fastq.gz",
        "1.cleandata/{sample}_R2.clean.fastq.gz",
    output:
        "1.cleandata/fastqc/{sample}_R1.clean_fastqc.zip",
        "1.cleandata/fastqc/{sample}_R2.clean_fastqc.zip",
        "1.cleandata/fastqc/{sample}_R1.clean_fastqc.html",
        "1.cleandata/fastqc/{sample}_R2.clean_fastqc.html",
    threads: 10
    shell:
        "fastqc -t {threads} -o 1.cleandata/fastqc/ {input}"

rule multiqc_clean:
    priority: -999
    output:
        "1.cleandata/fastqc/multiqc_report.html"
    shell:
        "multiqc -o 1.cleandata/fastqc/ 1.cleandata/fastqc/"

rule STAR:
    input:
        "1.cleandata/{sample}_R1.clean.fastq.gz",
        "1.cleandata/{sample}_R2.clean.fastq.gz",
    output:
        "2.mapping/{sample}_Aligned.sortedByCoord.out.bam"
    threads: 10
    log: "2.mapping/log/{sample}.log"
    shell:
        """
        STAR --genomeDir /picb/lilab5/liuzipei/placenta/WES/gatk_ref/STAR_index \
        --runThreadN {threads} \
        --readFilesIn  {input[0]} {input[1]} \
        --readFilesCommand zcat \
        --outFileNamePrefix 2.mapping/{wildcards.sample}_ \
        --outSAMstrandField intronMotif \
        --outSAMattributes All \
        --outFilterIntronMotifs RemoveNoncanonical \
        --outSAMtype BAM SortedByCoordinate > {log} 2>&1
        """

rule bam_index:
    input:
        "2.mapping/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        "2.mapping/{sample}_Aligned.sortedByCoord.out.bam.bai"
    threads: 10
    shell:
        "samtools index -@ {threads} {input}"


rule feature_count:
    priority: -999
    output:
        "3.featurecount/res.featureCounts"
    threads: 30
    params: strand=0
    log: "3.featurecount/featureCounts.log"
    shell:
        """
        featureCounts --primary -t exon -g gene_id -Q 20 -s {params.strand} \
        -p --countReadPairs -T {threads} \
        -a /picb/lilab5/share_data/ref/hg38/gencode.v44.annotation.gtf \
        -o {output} 2.mapping/*bam > {log} 2>&1
        """

rule get_count:
    input: 
        "3.featurecount/res.featureCounts"
    output: 
        "3.featurecount/count.tsv"
    shell: 
        """
        cut -f 1,7- {input} | sed 's|2.mapping/||g' | grep -v "#" > {output}
        """

rule TPM_FPKM:
    input: 
        "3.featurecount/count.tsv"
    output: 
        "3.featurecount/tpm_t.csv",
        "3.featurecount/fpkm_t.csv",
    shell: 
        """
        sh /home/liuzipei/scripts/RNAseq/tidy_count.sh
        """









