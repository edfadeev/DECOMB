#!/bin/bash
#
#SBATCH --job-name=assembly_MEGAHIT
#SBATCH --cpus-per-task=20
#SBATCH --mem=250GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --output=/proj/DECOMB/analysis/metaG_anvio/Log/%x-%j.out
#SBATCH --error=/proj/DECOMB/analysis/metaG_anvio/Log/%x-%j.err

#load module
module load megahit/1.1.2

#Set up the path to the working directory
WORKDIR=/proj/DECOMB/analysis/metaG_anvio
cd $WORKDIR

#run the program
megahit -1 $WORKDIR/01_QC/CDT3KANXX_105352_QC_R1.fastq.gz,$WORKDIR/01_QC/CDT3KANXX_105353_QC_R1.fastq.gz,$WORKDIR/01_QC/CDT3KANXX_105354_QC_R1.fastq.gz \
-2 $WORKDIR/01_QC/CDT3KANXX_105352_QC_R2.fastq.gz,$WORKDIR/01_QC/CDT3KANXX_105353_QC_R2.fastq.gz,$WORKDIR/01_QC/CDT3KANXX_105354_QC_R2.fastq.gz \
-o $WORKDIR/02_ASSEMBLY/megahit -t 20
