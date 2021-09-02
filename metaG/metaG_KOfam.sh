#!/bin/bash
#
#SBATCH --job-name=metaG_KOfam
#SBATCH --cpus-per-task=40
#SBATCH --mem=100GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --output=/proj/DECOMB/analysis/metaG_anvio/Log/%x-%j.out
#SBATCH --error=/proj/DECOMB/analysis/metaG_anvio/Log/%x-%j.err

#Set up the path to the working directory
WORKDIR=/proj/DECOMB/analysis/metaG_anvio/
  cd $WORKDIR


anvi-run-kegg-kofams -c $WORKDIR/05_ANVIO/spades.db \
--kegg-data-dir /proj/DECOMB/source/KOfam -T 40 \
--just-do-it

#produce tables with KO hits and modules
anvi-estimate-metabolism -c $WORKDIR/05_ANVIO/spades.db \
-O $WORKDIR/05_ANVIO/spades-Kofam\
--metagenome-mode
--kegg-data-dir /proj/DECOMB/source/KOfam \
--kegg-output-modes kofam_hits,modules

