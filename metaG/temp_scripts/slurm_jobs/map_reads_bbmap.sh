#!/bin/bash
#
#SBATCH --job-name=mapping-bbmap
#SBATCH --cpus-per-task=12
#SBATCH --mem=20GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --output=/proj/DECOMB/analysis/process_metaG/Log/%x-%j.out
#SBATCH --error=/proj/DECOMB/analysis/process_metaG/Log/%x-%j.err

#load module
module load bbmap/37.61
module load samtools/1.11

#Set up the path to the working directory
WORKDIR=/proj/DECOMB/analysis/process_metaG/
cd $WORKDIR

mkdir $WORKDIR/mapping/CDT3KANXX_${SLURM_ARRAY_TASK_ID}/

#reformat headers
anvi-script-reformat-fasta ${WORKDIR}/assembly/CDT3KANXX_${SLURM_ARRAY_TASK_ID}/scaffolds.fasta \
-o $WORKDIR/mapping/CDT3KANXX_${SLURM_ARRAY_TASK_ID}/scaffolds-fixed.fasta \
-l 1000 --simplify-names

#run the program
bbmap.sh in1=$WORKDIR/TRIM/CDT3KANXX_${SLURM_ARRAY_TASK_ID}_QC_R1.fastq.gz \
in2=$WORKDIR/TRIM/CDT3KANXX_${SLURM_ARRAY_TASK_ID}_QC_R2.fastq.gz \
ref=$WORKDIR/mapping/CDT3KANXX_${SLURM_ARRAY_TASK_ID}/scaffolds-fixed.fasta \
out=$WORKDIR/mapping/CDT3KANXX_${SLURM_ARRAY_TASK_ID}/scaffolds_mapped.sam \
statsfile=$WORKDIR/mapping/CDT3KANXX_${SLURM_ARRAY_TASK_ID}/scaffolds_mapped.stats \
covstats=$WORKDIR/mapping/CDT3KANXX_${SLURM_ARRAY_TASK_ID}/scaffolds_mapped.covstats \
rpkm=$WORKDIR/mapping/CDT3KANXX_${SLURM_ARRAY_TASK_ID}/scaffolds_mapped.rpkm \
nodisk=t \
Xmx=20G

#change SAM to BAM mapping files
samtools view -bS --threads 12 -o $WORKDIR/mapping/CDT3KANXX_${SLURM_ARRAY_TASK_ID}/scaffolds_mapped-RAW.bam \
$WORKDIR/mapping/CDT3KANXX_${SLURM_ARRAY_TASK_ID}/scaffolds_mapped.sam

#sort BAM files
anvi-init-bam --num-threads 12 $WORKDIR/mapping/CDT3KANXX_${SLURM_ARRAY_TASK_ID}/scaffolds_mapped-RAW.bam \
-o $WORKDIR/mapping/CDT3KANXX_${SLURM_ARRAY_TASK_ID}/scaffolds_mapped.bam
