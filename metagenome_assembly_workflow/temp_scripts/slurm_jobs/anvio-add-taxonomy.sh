#!/bin/bash
#
#SBATCH --job-name=anvio_taxonomy
#SBATCH --cpus-per-task=12
#SBATCH --mem=20GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --output=/proj/DECOMB/analysis/process_metaG/Log/%x-%j.out
#SBATCH --error=/proj/DECOMB/analysis/process_metaG/Log/%x-%j.err

#load module
module load diamond/0.9.36

#Set up the path to the working directory
WORKDIR=/proj/DECOMB/analysis/process_metaG/
cd $WORKDIR

#identify single copy genes
anvi-run-scg-taxonomy -c $WORKDIR/anvio/CDT3KANXX_${SLURM_ARRAY_TASK_ID}.db --num-threads 4 --num-parallel-processes 3


#estimate taxonomy
anvi-estimate-scg-taxonomy -c $WORKDIR/anvio/CDT3KANXX_${SLURM_ARRAY_TASK_ID}.db \
-p $WORKDIR/anvio/CDT3KANXX_${SLURM_ARRAY_TASK_ID}/PROFILE.db \
--num-threads 12 --metagenome-mode --compute-scg-coverages --update-profile-db-with-taxonomy \
--output-file $WORKDIR/anvio/CDT3KANXX_${SLURM_ARRAY_TASK_ID}/SCG_taxonomy.txt
