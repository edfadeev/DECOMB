#!/bin/bash
#
#SBATCH --job-name=anvio_taxonomy
#SBATCH --cpus-per-task=12
#SBATCH --mem=20GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --output=/proj/DECOMB/analysis/process_metaG/Log/%x-%j.out
#SBATCH --error=/proj/DECOMB/analysis/process_metaG/Log/%x-%j.err


#Set up the path to the working directory
WORKDIR=/proj/DECOMB/analysis/process_metaG/
cd $WORKDIR

#extimate taxonomy
anvi-estimate-scg-taxonomy -c $WORKDIR/anvio/CDT3KANXX_${SLURM_ARRAY_TASK_ID}.db \
--num-parallel-processes 12 \
-p $WORKDIR/anvio/CDT3KANXX_${SLURM_ARRAY_TASK_ID}/PROFILE.db \
--num-threads 12 --metagenome-mode --compute-scg-coverages --update-profile-db-with-taxonomy
