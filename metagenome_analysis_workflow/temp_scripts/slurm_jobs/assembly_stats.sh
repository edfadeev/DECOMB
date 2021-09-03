#!/bin/bash
#
#SBATCH --job-name=assembly_stats_Quast
#SBATCH --cpus-per-task=8
#SBATCH --mem=20GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --output=/proj/DECOMB/analysis/process_metaG/Log/%x-%j.out
#SBATCH --error=/proj/DECOMB/analysis/process_metaG/Log/%x-%j.err

#load module
module load quast/5.0.2
module load circos/0.69.9
module load ncbiblastplus/2.10.0

#Set up the path to the working directory
WORKDIR=/proj/DECOMB/analysis/process_metaG/
cd $WORKDIR

#run QUAST
metaquast.py $WORKDIR/assembly/CDT3KANXX_${SLURM_ARRAY_TASK_ID}/contigs.fasta \
--threads 8 --circos --glimmer --rna-finding \
--blast-db $WORKDIR/quast/silva_db/silva_db \
-1 $WORKDIR/assembly/CDT3KANXX_${SLURM_ARRAY_TASK_ID}/corrected/CDT3KANXX_${SLURM_ARRAY_TASK_ID}_QC_R1.fastq.00.0_0.cor.fastq.gz \
-2 $WORKDIR/assembly/CDT3KANXX_${SLURM_ARRAY_TASK_ID}/corrected/CDT3KANXX_${SLURM_ARRAY_TASK_ID}_QC_R2.fastq.00.0_0.cor.fastq.gz \
-o $WORKDIR/quast/CDT3KANXX_${SLURM_ARRAY_TASK_ID}
