#!/bin/bash
#
#SBATCH --job-name=import_to_anvio_megahit
#SBATCH --cpus-per-task=20
#SBATCH --mem=100GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --output=/proj/DECOMB/analysis/metaG_anvio/Log/%x-%j.out
#SBATCH --error=/proj/DECOMB/analysis/metaG_anvio/Log/%x-%j.err

#load module
#module load metabat/2.15
module load centrifuge/1.0.3-beta

#Set up the path to the working directory
WORKDIR=/proj/DECOMB/analysis/metaG_anvio/
cd $WORKDIR

#generate db
anvi-gen-contigs-database --num-threads 20 -f $WORKDIR/03_CONTIGS/megahit-contigs-fixed.fasta \
-o $WORKDIR/05_ANVIO/megahit.db --project-name megahit_assembly

#annotate using HMMs
anvi-run-hmms -c $WORKDIR/05_ANVIO/megahit.db --num-threads 20

#annotate COGs
anvi-run-ncbi-cogs -c $WORKDIR/05_ANVIO/megahit.db \
--num-threads 20 --cog-data-dir /proj/DECOMB/source/anvio-COG

#add gene taxonomies
anvi-get-sequences-for-gene-calls -c $WORKDIR/05_ANVIO/megahit.db -o $WORKDIR/05_ANVIO/megahit_gene_calls.fa
centrifuge -f -x /apps/centrifuge/1.0.3-beta/p+h+v_171103/p+h+v $WORKDIR/05_ANVIO/megahit_gene_calls.fa -S $WORKDIR/05_ANVIO/megahit_centrifuge_hits.tsv -p 20 \
--report-file $WORKDIR/05_ANVIO/megahit_centrifuge_reports.tsv

anvi-import-taxonomy-for-genes -c $WORKDIR/05_ANVIO/megahit.db \
-i $WORKDIR/05_ANVIO/megahit_centrifuge_reports.tsv $WORKDIR/05_ANVIO/megahit_centrifuge_hits.tsv \
-p centrifuge
