#!/bin/bash
#
#SBATCH --job-name=bam2fastq
#SBATCH --cpus-per-task=4
#SBATCH --mem=3GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --output=/proj/DECOMB/analysis/process_metaG/Log/%x-%j.out
#SBATCH --error=/proj/DECOMB/analysis/process_metaG/Log/%x-%j.err

#load module
module load tophat/2.1.1

#Set up the path to the working directory
WORKDIR=/proj/DECOMB/analysis/process_metaG/
cd $WORKDIR

#run the program
for file in ./unzipped/*.bam;

do bam2fastx --fastq --all -o ../$file.fastq $file

done 
