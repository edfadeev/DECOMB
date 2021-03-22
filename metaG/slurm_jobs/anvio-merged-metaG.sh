#!/bin/bash
#
#SBATCH --job-name=merge_metaG
#SBATCH --cpus-per-task=12
#SBATCH --mem=20GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --output=/proj/DECOMB/analysis/process_metaG/Log/%x-%j.out
#SBATCH --error=/proj/DECOMB/analysis/process_metaG/Log/%x-%j.err

#load module
module load metabat/2.15

#Set up the path to the working directory
WORKDIR=/proj/DECOMB/analysis/process_metaG/
cd $WORKDIR

#merge contigs
anvi-merge --contigs-db $WORKDIR/anvio/CDT3KANXX_*.db --sample-name merged_metaGs \
-o $WORKDIR/anvio/merged_metaG/merged_contigs.db

#generate anvio profile and map the reads from each metaG
anvi-profile -i $WORKDIR/mapping/CDT3KANXX_*/scaffolds_mapped.bam \
-c $WORKDIR/anvio/merged_metaG/merged_contigs.db \
--num-threads 12 \
--output-dir $WORKDIR/anvio/merged_metaG/MERGED_SAMPLES
