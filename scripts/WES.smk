# WES pipeline
# Data: 2024-05-10
# Require fastqc multiqc fastp bwa samtools R

Files = glob_wildcards("2.mapping/{sample}_sorted_markdup.bam")

rule all:
    input:
        expand("3.bam_qc/{sample}_depth.txt", sample=Files.sample),
        expand("3.bam_qc/{sample}_chimeric.txt", sample=Files.sample),
        expand("3.bam_qc/{sample}_insertsize.txt", sample=Files.sample),
        #expand("4.snp/gvcf/{sample}.HC.g.vcf.gz", sample=Files.sample),
        "3.bam_qc/all_samples_summary.txt",
        #"4.snp/vcf/merge_vcf/allsample.HC.vcf.gz",
        #"4.snp/vqsr/allsample.HC.VQSR.vcf",

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

rule BWA:
    input:
        "/picb/lilab5/liuzipei/placenta/WES_1/gatk_ref/Homo_sapiens_assembly38.fasta",
        "1.cleandata/{sample}_R1.clean.fastq.gz",
        "1.cleandata/{sample}_R2.clean.fastq.gz",
    output:
        "2.mapping/{sample}_sorted.bam"
    threads: 20
    log: "2.mapping/log/{sample}.log"
    shell:
        """
        bwa mem -t {threads} -R "@RG\\tID:{wildcards.sample}\\tLB:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:illumina\\tPU:{wildcards.sample}" \
        {input[0]} {input[1]} {input[2]} \
        | samtools view -@ {threads} -bS - \
        | samtools sort -@ {threads} - -o {output} > {log} 2>&1
        """

rule bam_index:
    input:
        "2.mapping/{sample}_sorted.bam"
    output:
        "2.mapping/{sample}_sorted.bam.bai"
    threads: 20
    shell:
        "samtools index -@ {threads} {input}"

rule bam_filter:
    input:
        "2.mapping/{sample}_sorted.bam"
    output:
        "2.mapping/{sample}_sorted.filter.bam"
    threads: 20
    shell:
        "samtools view -@ {threads} -h -f 0x2 {input} -o {output}"

rule markdup:
    input:
        "2.mapping/{sample}_sorted.filter.bam"
    output:
        "2.mapping/{sample}_sorted_markdup.bam",
        "2.mapping/{sample}_sorted_markdup_metrics.txt"
    threads: 20
    shell:
        """
        java -Xms60g -Xmx60g -XX:ParallelGCThreads={threads} \
        -Djava.io.tmpdir=/picb/lilab5/liuzipei/tmp_dir/ \
        -jar /home/liuxubing/software/picard.jar MarkDuplicates \
        -R /picb/lilab5/liuzipei/placenta/WES_1/gatk_ref/Homo_sapiens_assembly38.fasta \
        -I {input} -O {output[0]} -M {output[1]}
        """

rule markdup_bam_index:
    input:
        "2.mapping/{sample}_sorted_markdup.bam"
    output:
        "2.mapping/{sample}_sorted_markdup.bam.bai"
    threads: 20
    shell:
        "samtools index -@ {threads} {input}"

rule bam_qc_depth:
    input:
        "2.mapping/{sample}_sorted_markdup.bam"
    output:
        "3.bam_qc/{sample}_depth.txt"
    threads: 20
    shell:
        """
        samtools depth -@ {threads} -b /picb/lilab5/liuzipei/placenta/WES_1/5.somatic/filter/hg38.exons.refgene.bed\
        {input} > {output}
	"""

rule bam_qc_chimeric_reads:
    input:
        "2.mapping/{sample}_sorted_markdup.bam"
    output:
        "3.bam_qc/{sample}_chimeric.txt"
    threads: 20
    shell:
        """
        java -Xms60g -Xmx60g -XX:ParallelGCThreads={threads} \
        -Djava.io.tmpdir=/picb/lilab5/liuzipei/tmp_dir/ \
        -jar /home/liuxubing/software/picard.jar CollectJumpingLibraryMetrics \
        -R /picb/lilab5/liuzipei/placenta/WES_1/gatk_ref/Homo_sapiens_assembly38.fasta \
        -I {input} -O {output}  
        """

rule bam_qc_insert_size:
    input:
        "2.mapping/{sample}_sorted_markdup.bam"
    output:
        "3.bam_qc/{sample}_insertsize.txt",
        "3.bam_qc/{sample}_Histogram.pdf"
    threads: 20
    shell:
        """
        java -Xms60g -Xmx60g -XX:ParallelGCThreads={threads} \
        -Djava.io.tmpdir=/picb/lilab5/liuzipei/tmp_dir/ \
        -jar /home/liuxubing/software/picard.jar CollectInsertSizeMetrics \
        -I {input} -O {output[0]} \
        -R /picb/lilab5/liuzipei/placenta/WES_1/gatk_ref/Homo_sapiens_assembly38.fasta -H {output[1]}
        """

rule bam_qc_sum:
    priority: -999
    output:
        "3.bam_qc/all_samples_summary.txt"
    script:
        "/picb/lilab5/liuzipei/placenta/WES_1/scripts/WES_bam_qc_sum.sh"

rule gatk:
    input:
        "2.mapping/{sample}_sorted_markdup.bam",
        "2.mapping/{sample}_sorted_markdup.bam.bai"
    output:
        "4.snp/gvcf/{sample}.HC.g.vcf.gz",
        "4.snp/realign/{sample}_realign.bam"
    threads: 10
    shell:
        """
        /picb/lilab/tools/gatk-4.2.6.1/gatk --java-options -Xmx100G HaplotypeCaller \
        --tmp-dir /picb/lilab5/liuzipei/tmp_dir/ \
        -I {input[0]} \
        -R /picb/lilab5/liuzipei/placenta/WES_1/gatk_ref/Homo_sapiens_assembly38.fasta \
        -O {output[0]} \
        -bamout {output[1]} \
        --native-pair-hmm-threads {threads} \
        -ERC GVCF;
        """

rule gvcf_to_vcf:
    priority: -999
    output: "4.snp/vcf/merge_vcf/allsample.HC.vcf.gz",
    script:
        "/picb/lilab5/liuzipei/placenta/WES_1/scripts/gvcf_to_vcf.sh"

rule VQSR:
    input: "4.snp/vcf/merge_vcf/allsample.HC.vcf.gz",
    output: "4.snp/vqsr/allsample.HC.VQSR.vcf",
    script:
        "/picb/lilab5/liuzipei/placenta/WES_1/scripts/VQSR.sh"
