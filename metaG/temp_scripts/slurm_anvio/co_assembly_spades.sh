#!/bin/bash
#
#SBATCH --job-name=co-assembly_SPAdes
#SBATCH --cpus-per-task=40
#SBATCH --mem=450GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --output=/proj/DECOMB/analysis/metaG_anvio/Log/%x-%j.out
#SBATCH --error=/proj/DECOMB/analysis/metaG_anvio/Log/%x-%j.err

#load module
module load spades/3.14.1

#Set up the path to the working directory
WORKDIR=/proj/DECOMB/analysis/metaG_anvio
cd $WORKDIR

#merge fastq libraries
cat $WORKDIR/01_QC/*_QC_R1.fastq.gz > $WORKDIR/01_QC/R1.fastq.gz
cat $WORKDIR/01_QC/*_QC_R2.fastq.gz > $WORKDIR/01_QC/R2.fastq.gz

#run the program
spades.py --meta --tmp-dir $TMPDIR/tmp --threads 40 \
-1 $WORKDIR/01_QC/R1.fastq.gz -2 $WORKDIR/01_QC/R2.fastq.gz \
--memory 450 \
-o $WORKDIR/02_ASSEMBLY/SPAdes 

mv $TMPDIR/tmp ./assembly/corrected/tmp
