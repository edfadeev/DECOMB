#!/bin/bash
#
#SBATCH --job-name=fastqc
#SBATCH --cpus-per-task=4
#SBATCH --mem=3GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --output=/proj/DECOMB/analysis/process_metaG/Log/%x-%j.out
#SBATCH --error=/proj/DECOMB/analysis/process_metaG/Log/%x-%j.err

#load module
module load fastqc

#Set up the path to the working directory
WORKDIR=/proj/DECOMB/analysis/process_metaG/
cd $WORKDIR

#run the program
fastqc -t 4 --outdir ./QC/ CDT3KANXX_${SLURM_ARRAY_TASK_ID}.fastq 
