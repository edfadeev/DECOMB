#!/bin/bash
#
#SBATCH --job-name=merge_profiles_anvio_megahit
#SBATCH --cpus-per-task=20
#SBATCH --mem=10GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --output=/proj/DECOMB/analysis/metaG_anvio/Log/%x-%j.out
#SBATCH --error=/proj/DECOMB/analysis/metaG_anvio/Log/%x-%j.err

#Set up the path to the working directory
WORKDIR=/proj/DECOMB/analysis/metaG_anvio/
cd $WORKDIR

for i in {105352..105354}
do 
anvi-db-info $WORKDIR/05_ANVIO/megahit/CDT3KANXX_${i}/PROFILE.db \
--self-key sample_id --self-value "CDT3KANXX_${i}" --just-do-it
done

#merge profiles
anvi-merge $WORKDIR/05_ANVIO/megahit/*/PROFILE.db \
-o $WORKDIR/05_ANVIO/megahit/merged_profile \
-c $WORKDIR/05_ANVIO/megahit.db
