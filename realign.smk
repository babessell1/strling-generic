configfile: "config.yaml"

# Define a rule to generate the list of BAM files in the input directory
rule get_bam_files:
    output:
        "bam_files.txt"
    run:
        import os
        with open(output[0], "w") as f:
            bam_files = [f for f in os.listdir(config["BAM_DIR"]) if f.endswith(".bam")]
            f.write("\n".join(bam_files))

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
        bam = "str-results/{sample}_GRCh38.bam"
    log:
        "logs/bwa_mem_align_{wildcards.sample}.log"
    conda: "envs/bwa.yaml"
    shell:
        "bwa mem {input.ref} {input.fastq} | samtools sort -o {output.bam} &> {log}"

# Create a rule to create the output directories
rule create_output_dirs:
    output:
        directory(config["FASTQ_DIR"]),
        directory("str-results")

# Define a workflow to run all the steps
workflow realign_bam_to_new_ref:
    # First, create the output directories
    rule create_output_dirs

    # Second, get the list of BAM files in the directory
    rule get_bam_files

    # Third, convert BAM to FASTQ for each BAM file
    rule bam_to_fastq

    # Fourth, align the FASTQ files to the old reference for each sample
    rule bwa_mem_align

# Define the 'all' rule to run the entire workflow
rule all:
    input:
        expand("str-results/{sample}_GRCh38.bam", sample=read_bam_files())
