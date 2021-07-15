<<<<<<< HEAD
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
=======
#!/bin/bash
#
#SBATCH --job-name=bin_megahit_co-assembly
#SBATCH --cpus-per-task=20
#SBATCH --mem=100GB
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
anvi-cluster-contigs -p $WORKDIR/05_ANVIO/megahit/merged_profile/PROFILE.db \
-c $WORKDIR/05_ANVIO/megahit.db --collection-name CONCOCT --num-threads 20 \
--driver concoct --just-do-it 

#binning with metabat2
#anvi-cluster-contigs -p $WORKDIR/05_ANVIO/megahit/merged_profile/PROFILE.db \
#-c $WORKDIR/05_ANVIO/megahit.db --collection-name metabat2 --num-threads 20 \
#--driver metabat2 --just-do-it
>>>>>>> 92a014ab59d030391dd94e0174ecff8b34a348b6
