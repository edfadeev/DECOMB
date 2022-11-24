#!/bin/bash
#
#SBATCH --job-name=bins_checkM
#SBATCH --cpus-per-task=40
#SBATCH --mem=100GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --output=/scratch/oceanography/efadeev/DECOMB/analysis/metaG_anvio/Log/%x-%j.out
#SBATCH --error=/scratch/oceanography/efadeev/DECOMB/analysis/metaG_anvio/Log/%x-%j.err

#Set up the path to the working directory
WORKDIR=/scratch/oceanography/efadeev/DECOMB/analysis/metaG_anvio
  cd $WORKDIR

module load checkm/latest

checkm lineage_wf -t 20 -x fasta $WORKDIR/06_BINS/REFINED/ $WORKDIR/06_BINS/Refined_DAS_bins_checkM
