#!/bin/bash
#
#SBATCH --job-name=merge_metaG
#SBATCH --cpus-per-task=20
#SBATCH --mem=20GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --output=/proj/DECOMB/analysis/process_metaG/Log/%x-%j.out
#SBATCH --error=/proj/DECOMB/analysis/process_metaG/Log/%x-%j.err

#load module
#module load metabat/2.15

#Set up the path to the working directory
WORKDIR=/proj/DECOMB/analysis/process_metaG/
cd $WORKDIR

#merge scaffolds to one file
cat $WORKDIR/mapping/CDT3KANXX_*/scaffolds-fixed.fasta >> $WORKDIR/mapping/merged_scaffolds-fixed.fasta

#generate db
anvi-gen-contigs-database --num-threads 20 -f $WORKDIR/mapping/merged_scaffolds-fixed.fasta \
-o $WORKDIR/anvio/merged_metaG.db --project-name merged_metaG

#annotate using HMMs
anvi-run-hmms -c $WORKDIR/anvio/merged_metaG.db --num-threads 20

#annotate COGs
anvi-run-ncbi-cogs -c $WORKDIR/anvio/merged_metaG.db \
--num-threads 20 --cog-data-dir /proj/DECOMB/source/anvio-COG

#generate anvio profile and map the reads from each metaG
anvi-profile -i $WORKDIR/mapping/CDT3KANXX_*/scaffolds_mapped.bam \
-c $WORKDIR/anvio/merged_metaG.db \
--num-threads 20 \
--output-dir $WORKDIR/anvio/merged_metaG/MERGED_SAMPLES

#sxtimate taxonomy
anvi-estimate-scg-taxonomy -c $WORKDIR/anvio/merged_metaG.db \
--num-parallel-processes 5 \
-p $WORKDIR/anvio/merged_metaG/MERGED_SAMPLES/PROFILE.db \
--num-threads 20 --metagenome-mode --compute-scg-coverages --update-profile-db-with-taxonomy
