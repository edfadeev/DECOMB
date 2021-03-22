#!/bin/bash
#
#SBATCH --job-name=assembly_SPAdes
#SBATCH --cpus-per-task=20
#SBATCH --mem=250GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --output=/proj/DECOMB/analysis/process_metaG/Log/%x-%j.out
#SBATCH --error=/proj/DECOMB/analysis/process_metaG/Log/%x-%j.err

#load module
module load spades/3.14.1

#Set up the path to the working directory
WORKDIR=/proj/DECOMB/analysis/metaG_anvio
cd $WORKDIR

#run the program
spades.py --meta --tmp-dir $TMPDIR/tmp --threads 20 \
--pe1-1 $WORKDIR/01_QC/CDT3KANXX_105352_QC_R1.fastq.gz --pe1-2 $WORKDIR/01_QC/CDT3KANXX_105352_QC_R2.fastq.gz \
--pe2-1 $WORKDIR/01_QC/CDT3KANXX_105353_QC_R1.fastq.gz --pe2-2 $WORKDIR/01_QC/CDT3KANXX_105353_QC_R2.fastq.gz \
--pe3-1 $WORKDIR/01_QC/CDT3KANXX_105354_QC_R1.fastq.gz --pe3-2 $WORKDIR/01_QC/CDT3KANXX_105354_QC_R2.fastq.gz \
-o $WORKDIR/02_ASSEMBLY/SPAdes

mv $TMPDIR/tmp ./assembly/corrected/tmp
