#!/bin/bash
#
#SBATCH --job-name=bin_spades_co-assembly
#SBATCH --cpus-per-task=40
#SBATCH --mem=100GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --output=/home/project/oceanography/DECOMB/metaG_anvio/Log/%x-%j.out
#SBATCH --error=/home/project/oceanography/DECOMB/metaG_anvio/Log/%x-%j.err

#load module
module load conda
conda activate anvio_7.1

module unload python3/3.9.0
module load python3/3.7.15

module load metabat/2.15
module load concoct/1.1.0-py37

#Set up the path to the working directory
WORKDIR=/home/project/oceanography/DECOMB/metaG_anvio
cd $WORKDIR

#binning with concoct
anvi-cluster-contigs -p $WORKDIR/05_ANVIO/SPAdes/merged_profile/PROFILE.db \
-c $WORKDIR/05_ANVIO/spades.db --collection-name concoct_new --num-threads 40 \
--driver concoct --just-do-it

#binning with metabat2
anvi-cluster-contigs -p $WORKDIR/05_ANVIO/SPAdes/merged_profile/PROFILE.db \
-c $WORKDIR/05_ANVIO/spades.db --collection-name metabat2_new --num-threads 40 \
--driver metabat2 --just-do-it
