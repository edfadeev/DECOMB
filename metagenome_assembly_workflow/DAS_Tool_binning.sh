#!/bin/bash
#
#SBATCH --job-name=combine_binning_spades
#SBATCH --cpus-per-task=40
#SBATCH --mem=100GB
#SBATCH --mail-user=eduard.fadeev@univie.ac.at
#SBATCH --output=/proj/DECOMB/analysis/metaG_anvio/Log/%x-%j.out
#SBATCH --error=/proj/DECOMB/analysis/metaG_anvio/Log/%x-%j.err

#load module
module unload R
module load R/4.0.3
module load ncbiblast/2.2.26
module load dastool/1.1.0

#Set up the path to the working directory
WORKDIR=/proj/DECOMB/analysis/metaG_anvio
cd $WORKDIR

DAS_Tool -i $WORKDIR/06_BINS/DAS_Tool/spades_metabat2_contig.txt,$WORKDIR/06_BINS/DAS_Tool/spades_concoct_contig.txt -l metabat,concoct \
-c $WORKDIR/05_ANVIO/spades_contigs.fa --proteins $WORKDIR/06_BINS/DAS_Tool/spades-AAs-for-binning.fasta \
-o $WORKDIR/06_BINS/DAS_Tool/spades_ --search_engine blast --threads 40
