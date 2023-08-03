#!/bin/bash

#SBATCH --account=remills1
#SBATCH --job-name=realign
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem-per-cpu=8G
#SBATCH --time=48:00:00
#SBATCH --output=logs/realign.out
#SBATCH --error=logs/realign.err

eval "$(conda shell.bash hook)"
conda init bash
conda activate snake
module load Bioinformatics
module load samtools

snakemake -s realign.smk --unlock
snakemake -s realign.smk --rerun-incomplete --cores 12