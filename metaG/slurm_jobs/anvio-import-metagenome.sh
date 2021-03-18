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
module load samtools/1.11

#Set up the path to the working directory
WORKDIR=/proj/DECOMB/analysis/process_metaG/
cd $WORKDIR

#reformat headers
anvi-script-reformat-fasta ${WORKDIR}/assembly/CDT3KANXX_${SLURM_ARRAY_TASK_ID}/scaffolds.fasta \
-o $WORK_DIR/anvio/CDT3KANXX_${SLURM_ARRAY_TASK_ID}_scaffolds-fixed.fasta \
-l 1000 --simplify-names

#generate db
anvi-gen-contigs-database -f $WORK_DIR/anvio/CDT3KANXX_${SLURM_ARRAY_TASK_ID}_scaffolds-fixed.fasta \
-o $WORK_DIR/anvio/CDT3KANXX_${SLURM_ARRAY_TASK_ID}.db --project-name CDT3KANXX_${SLURM_ARRAY_TASK_ID}

#annotate using HMMs
anvi-run-hmms -c $WORK_DIR/anvio/CDT3KANXX_${SLURM_ARRAY_TASK_ID}.db --num-threads 12

#annotate COGs
anvi-run-ncbi-cogs -c $WORK_DIR/anvio/CDT3KANXX_${SLURM_ARRAY_TASK_ID}.db --num-threads 12

#Profiling BAM files
#sort BAM files
anvi-init-bam --num-threads 12 $WORK_DIR/mapping/$WORK_DIR/mapping/CDT3KANXX_${SLURM_ARRAY_TASK_ID}/scaffolds_mapped-RAW.bam \
-o $WORK_DIR/mapping/$WORK_DIR/mapping/CDT3KANXX_${SLURM_ARRAY_TASK_ID}/scaffolds_mapped.bam

#generate anvio profile
anvi-profile -i $WORK_DIR/mapping/$WORK_DIR/mapping/CDT3KANXX_${SLURM_ARRAY_TASK_ID}/scaffolds_mapped.bam \
-c $WORK_DIR/anvio/CDT3KANXX_${SLURM_ARRAY_TASK_ID}.db \
--num-threads 12 \
--sample-name ${SLURM_ARRAY_TASK_ID} \
--output-dir $WORK_DIR/anvio/CDT3KANXX_${SLURM_ARRAY_TASK_ID} \
