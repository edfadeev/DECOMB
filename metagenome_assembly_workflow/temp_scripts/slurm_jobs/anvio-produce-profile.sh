#!/bin/bash
#
#SBATCH --job-name=produce_anvio_profile
#SBATCH --cpus-per-task=20
#SBATCH --mem=20GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --output=/proj/DECOMB/analysis/process_metaG/Log/%x-%j.out
#SBATCH --error=/proj/DECOMB/analysis/process_metaG/Log/%x-%j.err

#load module

#Set up the path to the working directory
WORKDIR=/proj/DECOMB/analysis/process_metaG/
cd $WORKDIR

#generate anvio profile and map the reads from each metaG
anvi-profile -i $WORKDIR/mapping/CDT3KANXX_${SLURM_ARRAY_TASK_ID}/scaffolds_mapped.bam \
-c $WORKDIR/anvio/merged_metaG.db \
--num-threads 20 \
--sample-name CDT3KANXX_${SLURM_ARRAY_TASK_ID} \
--output-dir $WORKDIR/anvio/CDT3KANXX_${SLURM_ARRAY_TASK_ID}
