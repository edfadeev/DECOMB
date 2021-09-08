#!/bin/bash
#
#SBATCH --job-name=trimming-fastp
#SBATCH --cpus-per-task=8
#SBATCH --mem=12GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --output=/proj/DECOMB/analysis/process_metaG/Log/%x-%j.out
#SBATCH --error=/proj/DECOMB/analysis/process_metaG/Log/%x-%j.err

#load module
module load fastp

#Set up the path to the working directory
WORKDIR=/proj/DECOMB/analysis/process_metaG/
cd $WORKDIR

#run the program
fastp --interleaved_in --thread 8 --qualified_quality_phred 20 --unqualified_percent_limit 20 --length_required 50 -i ./fastq/CDT3KANXX_${SLURM_ARRAY_TASK_ID}.fastq --out1 ./TRIM/CDT3KANXX_${SLURM_ARRAY_TASK_ID}_QC_R1.fastq.gz --out2 ./TRIM/CDT3KANXX_${SLURM_ARRAY_TASK_ID}_QC_R2.fastq.gz -h ./TRIM/CDT3KANXX_${SLURM_ARRAY_TASK_ID}_QC.html -j ./TRIM/CDT3KANXX_${SLURM_ARRAY_TASK_ID}_QC.json --failed_out ./TRIM/CDT3KANXX_${SLURM_ARRAY_TASK_ID}_QC_failed.fastq.gz

