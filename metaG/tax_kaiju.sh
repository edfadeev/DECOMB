#!/bin/bash
#
#SBATCH --job-name=kaiju_annotation 
#SBATCH --cpus-per-task=40
#SBATCH --mem=100GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --output=/proj/DECOMB/analysis/metaG_anvio/Log/%x-%j.out
#SBATCH --error=/proj/DECOMB/analysis/metaG_anvio/Log/%x-%j.err

module load kaiju/1.7.3

#Set up the path to the working directory
WORKDIR=/proj/DECOMB/analysis/metaG_anvio/
  cd $WORKDIR

#kaiju -t /apps/kaiju/database/nodes.dmp \
#-f /apps/kaiju/database/refseq/kaiju_db_refseq.fmi \
#-i $WORKDIR/05_ANVIO/spades_gene_calls.fa  \
#-o $WORKDIR/05_ANVIO/spades_gene_calls_refseq_tax.out \
#-z 40 \
#-v

kaiju-addTaxonNames -t /apps/kaiju/database/nodes.dmp \
-n /apps/kaiju/database/names.dmp \
-i $WORKDIR/05_ANVIO/spades_gene_calls_refseq_tax.out \
-o $WORKDIR/05_ANVIO/spades_gene_calls_refseq_names.out \
-r superkingdom,phylum,order,class,family,genus,species

anvi-import-taxonomy-for-genes -i $WORKDIR/05_ANVIO/spades_gene_calls_refseq_names.out \
-c $WORKDIR/05_ANVIO/spades.db \
-p kaiju --just-do-it


