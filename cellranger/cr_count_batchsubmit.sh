#!/bin/bash

#SBATCH --job-name=cellranger_count

# Set the file to write the stdout and stderr to (if -e is not set; -o or --output).
#SBATCH --output=logs/%x-%j.log

# Set the number of cores (-n or --ntasks).
#SBATCH --ntasks=16

# Set the memory per CPU. Units can be given in T|G|M|K.
#SBATCH --mem-per-cpu=10G

# Set the partition to be used (-p or --partition).
#SBATCH --partition=medium
 
# Set the expected running time of your job (-t or --time).
# Formats are MM:SS, HH:MM:SS, Days-HH, Days-HH:MM, Days-HH:MM:SS
#SBATCH --time=12:00:00


#1 = name for run or name to give sample, e.g. 
#2 = sample number in format 001, 002 etc.


echo "----------------------------------------------------------------------------"
echo This is the processing of
echo $1
echo "Sample number"
echo $2
echo "----------------------------------------------------------------------------"

# get sample number and append to CG_FP_ to complete path
# e.g. CG_FP + 001
num=$2
sample="CG_FP_"
sample+=$num

cellranger count --id=$1 \
--fastqs=/fast/scratch/users/postmusd_c/scCHIK/P891/$sample/HW533DMXX/L002/ \
--transcriptome=/fast/work/users/postmusd_c/yard/genomes/hg38.chikvgfp/hg38_chikgfp_genome/ \
  --jobmode=slurm \
  --maxjobs=100 \
  --jobinterval=1000
