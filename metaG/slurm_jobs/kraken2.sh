#!/bin/bash
#
#SBATCH --job-name=kraken2-classification
#SBATCH --cpus-per-task=12
#SBATCH --mem=48GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --output=/proj/DECOMB/analysis/process_metaG/Log/%x-%j.out
#SBATCH --error=/proj/DECOMB/analysis/process_metaG/Log/%x-%j.err

#load module
module load kraken/2.0.8-beta

#Set up the path to the working directory
WORKDIR=/proj/DECOMB/analysis/process_metaG/
cd $WORKDIR

#run the program
kraken2 --db ./kraken2/krakendb --threads 12 --output ./kraken2/CDT3KANXX_${SLURM_ARRAY_TASK_ID} ./TRIM/CDT3KANXX_${SLURM_ARRAY_TASK_ID}_QC.fastq.gz
