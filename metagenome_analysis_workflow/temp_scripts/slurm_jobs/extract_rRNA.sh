#!/bin/bash
#
#SBATCH --job-name=extract_16S_rRNA
#SBATCH --cpus-per-task=4
#SBATCH --mem=8GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --output=/proj/DECOMB/analysis/process_metaG/Log/%x-%j.out
#SBATCH --error=/proj/DECOMB/analysis/process_metaG/Log/%x-%j.err

#load module
module load bbmap/37.61

#Set up the path to the working directory
WORKDIR=/proj/DECOMB/analysis/process_metaG/
cd $WORKDIR

#run the program
bbduk.sh in=./TRIM/CDT3KANXX_${SLURM_ARRAY_TASK_ID}_QC.fastq.gz outm=./16S_rRNA/CDT3KANXX_${SLURM_ARRAY_TASK_ID}_16S.fastq outu=CDT3KANXX_${SLURM_ARRAY_TASK_ID}_non_rRNA.fastq k=31 ref=./16S_rRNA/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz
