#!/bin/bash

#SBATCH --account=remills1
#SBATCH --job-name=realign
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem-per-cpu=4G
#SBATCH --time=24:00:00

eval "$(conda shell.bash hook)"
conda init bash
conda activate snake
module load Bioinformatics
module load samtools

snakemake -s realign.smk --unlock
snakemake -s realign.smk --rerun-incomplete --cores 12