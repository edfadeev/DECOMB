#!/bin/bash
#
#SBATCH --job-name=bins_checkM
#SBATCH --cpus-per-task=40
#SBATCH --mem=100GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --output=/proj/DECOMB/analysis/metaG_anvio/Log/%x-%j.out
#SBATCH --error=/proj/DECOMB/analysis/metaG_anvio/Log/%x-%j.err

#Set up the path to the working directory
WORKDIR=/proj/DECOMB/analysis/metaG_anvio/
  cd $WORKDIR

module load checkm/latest

#checkm taxonomy_wf order Alteromonadales -x fa -t 40 $WORKDIR/06_BINS/Selected_bins_summary/Selected_bins_fasta $WORKDIR/06_BINS/Selected_bins_summary/CheckM_Alteromonadales

#checkm taxonomy_wf order Rhodobacterales -x fa -t 40 $WORKDIR/06_BINS/Selected_bins_summary/Selected_bins_fasta $WORKDIR/06_BINS/Selected_bins_summary/CheckM_Rhodobacterales

#checkm taxonomy_wf order Vibrionales -x fa -t 40 $WORKDIR/06_BINS/Selected_bins_summary/Selected_bins_fasta $WORKDIR/06_BINS/Selected_bins_summary/CheckM_Vibrionales

#checkm taxonomy_wf order Cellvibrionales -x fa -t 40 $WORKDIR/06_BINS/Selected_bins_summary/Selected_bins_fasta $WORKDIR/06_BINS/Selected_bins_summary/CheckM_Cellvibrionales

checkm taxonomy_wf order Rhodospirillales -x fa -t 40 $WORKDIR/06_BINS/Selected_bins_summary/Selected_bins_fasta $WORKDIR/06_BINS/Selected_bins_summary/CheckM_Rhodospirillales


