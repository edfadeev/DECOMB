#!/bin/bash
#
#SBATCH --job-name=bin_spades_co-assembly
#SBATCH --cpus-per-task=40
#SBATCH --mem=400GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --output=/proj/DECOMB/analysis/metaG_anvio/Log/%x-%j.out
#SBATCH --error=/proj/DECOMB/analysis/metaG_anvio/Log/%x-%j.err

#load module
module unload python3/3.9.0
module load python3/3.7.0

module load metabat/2.15
module load concoct/1.1.0-py37

#Set up the path to the working directory
WORKDIR=/proj/DECOMB/analysis/metaG_anvio
cd $WORKDIR

#binning with concoct
anvi-cluster-contigs -p $WORKDIR/05_ANVIO/SPAdes/merged_profile/PROFILE.db \
-c $WORKDIR/05_ANVIO/spades.db --collection-name CONCOCT --num-threads 40 \
--driver concoct --just-do-it

anvi-summarize -p $WORKDIR/05_ANVIO/SPAdes/merged_profile/PROFILE.db -c $WORKDIR/05_ANVIO/spades.db -C CONCOCT \
-o $WORKDIR/06_BINS/concoct_summary

#binning with metabat2
anvi-cluster-contigs -p $WORKDIR/05_ANVIO/SPAdes/merged_profile/PROFILE.db \
-c $WORKDIR/05_ANVIO/spades.db --collection-name metabat2 --num-threads 40 \
--driver metabat2 --just-do-it

anvi-summarize -p $WORKDIR/05_ANVIO/SPAdes/merged_profile/PROFILE.db -c $WORKDIR/05_ANVIO/spades.db -C metabat2 \
-o $WORKDIR/06_BINS/metabat2_summary
