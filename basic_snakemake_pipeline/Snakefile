# RNAseq snakemake pipeline

"""
Rules for analysing fasta files with FastQC, Trimmomatic, hisat2 and htseq2.
For usage, include this in your workflow.
"""

__author__ = "Theo Killian (theodore.killian@uclouvain.be)"
__license__ = "Artistic-2.0"


directories = ["reports", "reports/fastqc1", "reports/fastqc2", "reports/trimmo", "reports/hisat2", "reports/htseq"]

fastqc_input1 = ["data/010419E18_R1.fastq.gz", "data/010419E18_R2.fastq.gz"]

fastqc_input2 = ["reports/trimmo/010419E18_R1_paired.fastq.gz", "reports/trimmo/010419E18_R2_paired.fastq.gz"]

fastqc_output1 = ["reports/fastqc1/010419E18_R1_fastqc.html", "reports/fastqc1/010419E18_R2_fastqc.html"]

fastqc_output2 = ["reports/fastqc2/010419E18_R1_paired_fastqc.html", "reports/fastqc2/010419E18_R2_paired_fastqc.html"]

trimmo_output = ["reports/trimmo/010419E18_R1_paired.fastq.gz", "reports/trimmo/010419E18_R2_paired.fastq.gz",
    "reports/trimmo/010419E18_R1_unpaired.fastq.gz", "reports/trimmo/010419E18_R2_unpaired.fastq.gz"]

hisat2_output = ["reports/hisat2/010419E18.sam"]

sam2bam_output = ["data/010419E18.bam"]

sort_bam_output = ["data/010419E18_sorted.bam"]

index_bam_output = ["data/010419E18_sorted.bam.bai"]

htseq_count_output = ["reports/htseq/010419E18.tsv"]

count_matrix_output = ["reports/htseq/010419E18_counts"]

ASSEMBLY_VERSION = "genome_snp_tran"

THREADS = 8

rule all:
    input:
        directories,
        fastqc_output1,
        trimmo_output,
        fastqc_output2,
        hisat2_output,
        sam2bam_output,
        sort_bam_output,
        index_bam_output,
        htseq_count_output,
        count_matrix_output

rule make_directories:
    shell:
        """
        mkdir reports
        mkdir reports/fastqc1
        mkdir reports/fastqc2
        mkdir reports/trimmo
        mkdir reports/hisat2
        mkdir reports/htseq
        """

rule fastqc_1:
    input:
        fastqc_input1
    output:
        html="reports/fastqc1/{sample}_fastqc.html",
        zip="reports/fastqc1/{sample}_fastqc.zip"
    params: ""
    threads: THREADS
    log:
        "reports/fastqc1/{sample}.log"
    wrapper:
        "0.36.0/bio/fastqc"

rule trimmomatic_pe:
    input:
        r1="data/{sample}_R1.fastq.gz",
        r2="data/{sample}_R2.fastq.gz"
    output:
        r1="reports/trimmo/{sample}_R1_paired.fastq.gz",
        r2="reports/trimmo/{sample}_R2_paired.fastq.gz",
        r1_unpaired="reports/trimmo/{sample}_R1_unpaired.fastq.gz",
        r2_unpaired="reports/trimmo/{sample}_R2_unpaired.fastq.gz"
    log:
        "logs/trimmomatic/{sample}.log"
    params:
        # list of trimmers (see manual)
        trimmer=["ILLUMINACLIP:adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:10:20 MINLEN:50"],
        # optional parameters
        extra="",
        compression_level="-9"
    threads: THREADS
    wrapper:
        "0.36.0/bio/trimmomatic/pe"

rule fastqc_2:
    input:
        fastqc_input2
    output:
        html="reports/fastqc2/{sample}_paired_fastqc.html",
        zip="reports/fastqc2/{sample}_paired_fastqc.zip"
    params: ""
    threads: THREADS
    log:
        "reports/fastqc2/{sample}_paired.log"
    wrapper:
        "0.36.0/bio/fastqc"

rule hisat2_align:
    input:
        r1="reports/trimmo/{sample}_R1_paired.fastq.gz",
        r2="reports/trimmo/{sample}_R2_paired.fastq.gz"
    output:
        sam="reports/hisat2/{sample}.sam",
        log="reports/hisat2/{sample}.log"
    threads: THREADS
    shell:
        """
        hisat2 -p 6 --remove-chrname --dta --summary-file {output.log} -x INDEXES/{ASSEMBLY_VERSION} -1 {input.r1} -2 {input.r2} -S {output.sam}
        """

rule sam2bam:
    input:
        "reports/hisat2/{sample}.sam"
    output:
        "data/{sample}.bam"
    threads: THREADS
    shell:
        """
        samtools view -S -b {input} > {output} 
        """

rule sort_bam:
    input:
        "data/{sample}.bam"
    output:
        "data/{sample}_sorted.bam"
    threads: THREADS
    shell:
        """
        samtools sort {input} > {output} 
        """

rule index_bam:
    input:
        "data/{sample}_sorted.bam"
    output:
        "data/{sample}_sorted.bam.bai"
    threads: THREADS
    shell:
        """
        samtools index {input}
        """

rule htseq_count:
    input:
        sam = "reports/hisat2/{sample}.sam",
        gff = "GTFs/Mus_musculus.GRCm38.94.gtf.gz"
    output:
        #"reports/htseq/{sample}.log"
        "reports/htseq/{sample}.tsv"
    threads: THREADS
    shell:
        """
        htseq-count -q -s no {input.sam} {input.gff} > {output}
        """

rule count_matrix:
    input:
        "reports/htseq/{sample}.tsv"
    output:
        counts="reports/htseq/{sample}_counts"
    threads: THREADS
    shell:
        """
        cat {input} | egrep "^ENS" | cut -f2 >> {output.counts}
        """

