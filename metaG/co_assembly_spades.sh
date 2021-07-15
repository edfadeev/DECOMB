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
module load quast/5.0.2
module load circos/0.69.9
module load ncbiblastplus/2.10.0


#Set up the path to the working directory
WORKDIR=/proj/DECOMB/analysis/metaG_anvio
cd $WORKDIR

#merge fastq libraries
cat $WORKDIR/01_QC/*_QC_R1.fastq.gz > $WORKDIR/01_QC/R1.fastq.gz
cat $WORKDIR/01_QC/*_QC_R2.fastq.gz > $WORKDIR/01_QC/R2.fastq.gz

#run the program
spades.py --meta --tmp-dir $TMPDIR/tmp --threads 20 \
-1 $WORKDIR/01_QC/R1.fastq.gz -2 $WORKDIR/01_QC/R2.fastq.gz \
-o $WORKDIR/02_ASSEMBLY/SPAdes

mv $TMPDIR/tmp ./assembly/corrected/tmp

#extract contigs
mkdir $WORKDIR/03_CONTIGS

#correct deflines for Anvio
anvi-script-reformat-fasta $WORKDIR/02_ASSEMBLY/SPAdes/contigs.fasta \
-o $WORKDIR/03_CONTIGS/spades-contigs-fixed.fasta \
-l 1000 --simplify-names

#run QUAST
metaquast.py $WORKDIR/03_CONTIGS/spades-contigs-fixed.fasta \
--threads 20 --circos --mgm --rna-finding \
--blast-db /proj/DECOMB/source/silva_db/silva_db \
-1 $WORKDIR/01_QC/R1.fastq.gz \
-2 $WORKDIR/01_QC/R2.fastq.gz \
-o $WORKDIR/03_CONTIGS/SPAdes

