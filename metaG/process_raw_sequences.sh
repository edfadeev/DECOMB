#!/bin/bash
#
#SBATCH --job-name=process_raw_sequences
#SBATCH --cpus-per-task=8
#SBATCH --mem=12GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --output=/proj/DECOMB/analysis/metaG_anvio/Log/%x-%j.out
#SBATCH --error=/proj/DECOMB/analysis/metaG_anvio/Log/%x-%j.err

#load module
module load tophat/2.1.1
module load fastqc
module load fastp

#Set up the path to the working directory
WORKDIR=/proj/DECOMB/analysis/metaG_anvio/
cd $WORKDIR

#export sequences from BAM files
for file in ./unzipped/*.bam;

do bam2fastx --fastq --all -o ../$file.fastq $file

done 

#move the sequences to a new directory
mkdir $WORKDIR/01_QC
mv *.fastq $WORKDIR/01_QC/

#run fastQC
fastqc -t 8 --outdir $WORKDIR/01_QC/ $WORKDIR/01_QC/CDT3KANXX_${SLURM_ARRAY_TASK_ID}.fastq 

#run fastP
fastp --interleaved_in --thread 8 --qualified_quality_phred 20 --unqualified_percent_limit 20 --length_required 50 -i $WORKDIR/01_QC/CDT3KANXX_${SLURM_ARRAY_TASK_ID}.fastq --out1 $WORKDIR/01_QC/CDT3KANXX_${SLURM_ARRAY_TASK_ID}_QC_R1.fastq.gz --out2 $WORKDIR/01_QC/CDT3KANXX_${SLURM_ARRAY_TASK_ID}_QC_R2.fastq.gz -h $WORKDIR/01_QC/CDT3KANXX_${SLURM_ARRAY_TASK_ID}_QC.html -j $WORKDIR/01_QC/CDT3KANXX_${SLURM_ARRAY_TASK_ID}_QC.json --failed_out $WORKDIR/01_QC/CDT3KANXX_${SLURM_ARRAY_TASK_ID}_QC_failed.fastq.gz




