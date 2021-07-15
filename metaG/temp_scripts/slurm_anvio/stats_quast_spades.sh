#!/bin/bash
#
#SBATCH --job-name=SPAdes_stats
#SBATCH --cpus-per-task=20
#SBATCH --mem=20GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --output=/proj/DECOMB/analysis/metaG_anvio/Log/%x-%j.out
#SBATCH --error=/proj/DECOMB/analysis/metaG_anvio/Log/%x-%j.err

#load module
module load quast/5.0.2
module load circos/0.69.9
module load ncbiblastplus/2.10.0

#Set up the path to the working directory
WORKDIR=/proj/DECOMB/analysis/metaG_anvio/
cd $WORKDIR

#run QUAST
metaquast.py $WORKDIR/03_CONTIGS/contigs.fasta \
--threads 20 --circos --mgm --rna-finding \
--blast-db /proj/DECOMB/source/silva_db/silva_db \
-1 $WORKDIR/01_QC/R1.fastq.gz \
-2 $WORKDIR/01_QC/R2.fastq.gz \
-o $WORKDIR/03_CONTIGS/SPAdes
