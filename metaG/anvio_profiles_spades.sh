#!/bin/bash
#
#SBATCH --job-name=merge_profiles_anvio_spades
#SBATCH --cpus-per-task=20
#SBATCH --mem=100GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --output=/proj/DECOMB/analysis/metaG_anvio/Log/%x-%j.out
#SBATCH --error=/proj/DECOMB/analysis/metaG_anvio/Log/%x-%j.err

#Set up the path to the working directory
WORKDIR=/proj/DECOMB/analysis/metaG_anvio/
cd $WORKDIR

for i in {105352..105354}
do 
  anvi-profile -i $WORKDIR/04_MAPPING/SPAdes/CDT3KANXX_${i}/scaffolds_mapped.bam \
  -c $WORKDIR/05_ANVIO/spades.db -T 20 -o $WORKDIR/05_ANVIO/SPAdes/CDT3KANXX_${i}
done

#merge profiles
anvi-merge $WORKDIR/05_ANVIO/SPAdes/*/PROFILE.db \
-o $WORKDIR/05_ANVIO/SPAdes/merged_profile \
-c $WORKDIR/05_ANVIO/spades.db
