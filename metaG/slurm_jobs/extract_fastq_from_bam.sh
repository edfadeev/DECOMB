#!/bin/bash
#
#SBATCH --job-name=bam2fastq
#SBATCH --cpus-per-task=4
#SBATCH --mem=3GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --output=log/%x-%j.out
#SBATCH --error=log/%x-%j.err

#load module
module load tophat/2.1.1

#Set up the path to the working directory
WORKDIR=/proj/DECOMB/data/process_metaG/
cd $WORKDIR

#run the program
bam2fastx --fastq --all -o CDT3KANXX_2.fastq ./unzipped/CDT3KANXX_2\#0.bam
