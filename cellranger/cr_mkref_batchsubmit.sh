#!/bin/bash

#SBATCH --job-name=cellranger_mkref

# Set the file to write the stdout and stderr to (if -e is not set; -o or --output).
#SBATCH --output=logs/%x-%j.log

# Set the number of cores (-n or --ntasks).
#SBATCH --ntasks=16

# Set the memory per CPU. Units can be given in T|G|M|K.
#SBATCH --mem-per-cpu=8G

# Set the partition to be used (-p or --partition).
#SBATCH --partition=medium
 
# Set the expected running time of your job (-t or --time).
# Formats are MM:SS, HH:MM:SS, Days-HH, Days-HH:MM, Days-HH:MM:SS
#SBATCH --time=12:00:00

cellranger mkref \
--genome=hg38_chikgfp_genome \
--fasta=hg38.103_chikvgfp.fa \
--genes=hg38.102_chikgfp.gtf
