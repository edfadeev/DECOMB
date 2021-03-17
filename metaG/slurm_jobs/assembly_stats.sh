#!/bin/bash
#
#SBATCH --job-name=assembly_stats_Quast
#SBATCH --cpus-per-task=8
#SBATCH --mem=20GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --output=/proj/DECOMB/analysis/process_metaG/Log/%x-%j.out
#SBATCH --error=/proj/DECOMB/analysis/process_metaG/Log/%x-%j.err

#load module
module unload python3
module load quast/5.0.2
module load circos/0.69.9

#Set up the path to the working directory
WORKDIR=/proj/DECOMB/analysis/process_metaG/
cd $WORKDIR

#run QUAST
metaquast.py ./assembly/CDT3KANXX_${SLURM_ARRAY_TASK_ID}/contigs.fasta \
--threads 8 --circos --glimmer --rna-finding \
/mirror/silva/current/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz
-1 ./assembly/CDT3KANXX_${SLURM_ARRAY_TASK_ID}/corrected/CDT3KANXX_${SLURM_ARRAY_TASK_ID}_QC_R1.fastq.00.0_0.cor.fastq.gz \
-2 ./assembly/CDT3KANXX_${SLURM_ARRAY_TASK_ID}/corrected/CDT3KANXX_${SLURM_ARRAY_TASK_ID}_QC_R2.fastq.00.0_0.cor.fastq.gz \
-o ./quast/CDT3KANXX_${SLURM_ARRAY_TASK_ID}
