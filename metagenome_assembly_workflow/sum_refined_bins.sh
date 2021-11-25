#!/bin/bash
#
#SBATCH --job-name=summarize_refined_bins
#SBATCH --cpus-per-task=10
#SBATCH --mem=40GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --output=/proj/DECOMB/analysis/metaG_anvio/Log/%x-%j.out
#SBATCH --error=/proj/DECOMB/analysis/metaG_anvio/Log/%x-%j.err

#Set up the path to the working directory
WORKDIR=/proj/DECOMB/analysis/metaG_anvio
cd $WORKDIR

#summarize refined bins
anvi-summarize -p $WORKDIR/05_ANVIO/SPAdes/merged_profile/PROFILE.db -c $WORKDIR/05_ANVIO/spades.db -C DAS_Tool \
-o $WORKDIR/06_BINS/DAS_Tool_refined

#export each bin into a separate db
anvi-split -p $WORKDIR/05_ANVIO/SPAdes/merged_profile/PROFILE.db -c $WORKDIR/05_ANVIO/spades.db \
--collection-name DAS_Tool -o $WORKDIR/06_BINS/REFINED
