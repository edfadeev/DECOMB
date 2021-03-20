#!/bin/bash
#
#SBATCH --job-name=import_to_anvio
#SBATCH --cpus-per-task=12
#SBATCH --mem=20GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --output=/proj/DECOMB/analysis/process_metaG/Log/%x-%j.out
#SBATCH --error=/proj/DECOMB/analysis/process_metaG/Log/%x-%j.err

#load module
module load anvio/7

#Set up the path to the working directory
WORKDIR=/proj/DECOMB/analysis/process_metaG/
cd $WORKDIR

#reformat headers
anvi-script-reformat-fasta ${WORKDIR}/assembly/CDT3KANXX_${SLURM_ARRAY_TASK_ID}/scaffolds.fasta \
-o $WORKDIR/anvio/CDT3KANXX_${SLURM_ARRAY_TASK_ID}_scaffolds-fixed.fasta \
-l 1000 --simplify-names

#generate db
anvi-gen-contigs-database -f $WORKDIR/anvio/CDT3KANXX_${SLURM_ARRAY_TASK_ID}_scaffolds-fixed.fasta \
-o $WORKDIR/anvio/CDT3KANXX_${SLURM_ARRAY_TASK_ID}.db --project-name CDT3KANXX_${SLURM_ARRAY_TASK_ID}

#annotate using HMMs
anvi-run-hmms -c $WORKDIR/anvio/CDT3KANXX_${SLURM_ARRAY_TASK_ID}.db --num-threads 12

#annotate COGs
anvi-run-ncbi-cogs -c $WORKDIR/anvio/CDT3KANXX_${SLURM_ARRAY_TASK_ID}.db \
--num-threads 12 --cog-data-dir /proj/DECOMB/source/anvio-COG

#Profiling BAM files
#sort BAM files
anvi-init-bam --num-threads 12 $WORKDIR/mapping/CDT3KANXX_${SLURM_ARRAY_TASK_ID}/scaffolds_mapped-RAW.bam \
-o $WORKDIR/mapping/CDT3KANXX_${SLURM_ARRAY_TASK_ID}/scaffolds_mapped.bam

#generate anvio profile
anvi-profile -i $WORKDIR/mapping/CDT3KANXX_${SLURM_ARRAY_TASK_ID}/scaffolds_mapped.bam \
-c $WORKDIR/anvio/CDT3KANXX_${SLURM_ARRAY_TASK_ID}.db \
--num-threads 12 \
--sample-name CDT3KANXX_${SLURM_ARRAY_TASK_ID} \
--output-dir $WORKDIR/anvio/CDT3KANXX_${SLURM_ARRAY_TASK_ID}
