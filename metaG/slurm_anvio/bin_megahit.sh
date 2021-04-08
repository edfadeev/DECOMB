#!/bin/bash
#
#SBATCH --job-name=bin_megahit_co-assembly
#SBATCH --cpus-per-task=40
#SBATCH --mem=200GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --output=/proj/DECOMB/analysis/metaG_anvio/Log/%x-%j.out
#SBATCH --error=/proj/DECOMB/analysis/metaG_anvio/Log/%x-%j.err

#load module
module load metabat/2.15
module load concoct/1.1.0

#Set up the path to the working directory
WORKDIR=/proj/DECOMB/analysis/metaG_anvio
cd $WORKDIR

#binning with concoct
anvi-cluster-contigs $WORKDIR/05_ANVIO/megahit/merged_profile/PROFILE.db \
-c $WORKDIR/05_ANVIO/megahit.db --collection-name CONCOCT --num-threads 40 \
--driver concoct 

#binning with metabat2
anvi-cluster-contigs $WORKDIR/05_ANVIO/megahit/merged_profile/PROFILE.db \
-c $WORKDIR/05_ANVIO/megahit.db --collection-name metabat2 --num-threads 40 \
--driver metabat2 
