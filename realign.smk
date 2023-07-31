configfile: "config.yaml"


# Define a rule to convert BAM to FASTQ using bedtools bamtofastq
rule bam_to_fastq:
    input:
        bam = config["BAM_DIR"] + "{sample}.bam"
    output:
        fastq = config["FASTQ_DIR"] + "{sample}.fq"
    log:
        "logs/bam_to_fastq_{wildcards.sample}.log"
    conda: "envs/bed.yaml"
    shell:
        "bedtools bamtofastq -i {input.bam} -fq {output.fastq} &> {log}"

# Define a rule to align the FASTQ files to the old reference using bwa mem
rule bwa_mem_align:
    input:
        fastq = config["FASTQ_DIR"] + "{sample}.fq",
        ref = config["REF_FASTA"]
    output:
        bam = "realigned_bams/{sample}_GRCh38.bam"
    log:
        "logs/bwa_mem_align_{wildcards.sample}.log"
    conda: "envs/bwa.yaml"
    shell:
        """
        mkdir -p realigned_bams
        bwa mem {input.ref} {input.fastq} | samtools sort -o {output.bam} &> {log}
        """n

# Define the 'all' rule to run the entire workflow
rule all:
    input:
        expand(config["FASTQ_DIR"] + "{sample}.fq", sample=config["SAMPLES"]),
        expand("realigned_bams/{sample}_GRCh38.bam", sample=read_bam_files())
