#!/bin/bash
#
#SBATCH --job-name=assembly_SPAdes
#SBATCH --cpus-per-task=16
#SBATCH --mem=250GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --output=/proj/DECOMB/analysis/process_metaG/Log/%x-%j.out
#SBATCH --error=/proj/DECOMB/analysis/process_metaG/Log/%x-%j.err

#load module
module load spades/3.14.1

#Set up the path to the working directory
WORKDIR=/proj/DECOMB/analysis/process_metaG/
cd $WORKDIR

#run the program
spades.py --meta --tmp-dir $TMPDIR/tmp --threads 16 -1 ./TRIM/CDT3KANXX_${SLURM_ARRAY_TASK_ID}_QC_R1.fastq.gz -2 ./TRIM/CDT3KANXX_${SLURM_ARRAY_TASK_ID}_QC_R2.fastq.gz -o ./assembly/CDT3KANXX_${SLURM_ARRAY_TASK_ID}

mv $TMPDIR/tmp ./assembly/corrected/tmp
