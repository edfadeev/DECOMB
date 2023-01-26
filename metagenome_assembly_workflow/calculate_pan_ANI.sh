#!/bin/bash
#
#SBATCH --job-name=Nucleotide_identitiy_calculations
#SBATCH --cpus-per-task=20
#SBATCH --mem=40GB
#SBATCH --partition=basic
#SBATCH --mail-type=ALL
#SBATCH --output=/scratch/oceanography/efadeev/DECOMB/analysis/Pphenol_pangenomics/00_LOGS/%x-%j.out
#SBATCH --time=1-24:00:00

#Set up the path to the working directory
WORKDIR=/scratch/oceanography/efadeev/DECOMB/analysis/Pphenol_pangenomics/
cd $WORKDIR

#run calculation
anvi-compute-genome-similarity -e $WORKDIR/Pphenol_genomes.txt \
-o $WORKDIR/03_PAN/ANI \
-p $WORKDIR/03_PAN/P_phenolica_pangenome-PAN.db \
--program pyANI \
-T 20