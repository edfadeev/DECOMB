#!/bin/bash
#
#SBATCH --job-name=mapping-bbmap_spades
#SBATCH --cpus-per-task=20
#SBATCH --mem=100GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --output=/proj/DECOMB/analysis/metaG_anvio/Log/%x-%j.out
#SBATCH --error=/proj/DECOMB/analysis/metaG_anvio/Log/%x-%j.err

#load module
module load bbmap/37.61
module load samtools/1.11

#Set up the path to the working directory
WORKDIR=/proj/DECOMB/analysis/metaG_anvio
cd $WORKDIR

mkdir $WORKDIR/04_MAPPING/spades/CDT3KANXX_${SLURM_ARRAY_TASK_ID}/

#run the program
bbmap.sh in1=$WORKDIR/01_QC/CDT3KANXX_${SLURM_ARRAY_TASK_ID}_QC_R1.fastq.gz \
in2=$WORKDIR/01_QC/CDT3KANXX_${SLURM_ARRAY_TASK_ID}_QC_R2.fastq.gz \
ref=$WORKDIR/03_CONTIGS/spades-contigs-fixed.fasta \
out=$WORKDIR/04_MAPPING/SPAdes/CDT3KANXX_${SLURM_ARRAY_TASK_ID}/scaffolds_mapped.sam \
statsfile=$WORKDIR/04_MAPPING/SPAdes/CDT3KANXX_${SLURM_ARRAY_TASK_ID}/scaffolds_mapped.stats \
covstats=$WORKDIR/04_MAPPING/SPAdes/CDT3KANXX_${SLURM_ARRAY_TASK_ID}/scaffolds_mapped.covstats \
rpkm=$WORKDIR/04_MAPPING/SPAdes/CDT3KANXX_${SLURM_ARRAY_TASK_ID}/scaffolds_mapped.rpkm \
nodisk=t \
Xmx=100G

#change SAM to BAM mapping files
samtools view -bS --threads 20 -o $WORKDIR/04_MAPPING/SPAdes/CDT3KANXX_${SLURM_ARRAY_TASK_ID}/scaffolds_mapped-RAW.bam \
$WORKDIR/04_MAPPING/SPAdes/CDT3KANXX_${SLURM_ARRAY_TASK_ID}/scaffolds_mapped.sam

#sort BAM files
anvi-init-bam --num-threads 20 $WORKDIR/04_MAPPING/SPAdes/CDT3KANXX_${SLURM_ARRAY_TASK_ID}/scaffolds_mapped-RAW.bam \
-o $WORKDIR/04_MAPPING/SPAdes/CDT3KANXX_${SLURM_ARRAY_TASK_ID}/scaffolds_mapped.bam
