#!/bin/bash
#
#SBATCH --job-name=centrifuge_tax_annotation 
#SBATCH --cpus-per-task=40
#SBATCH --mem=100GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --output=/proj/DECOMB/analysis/metaG_anvio/Log/%x-%j.out
#SBATCH --error=/proj/DECOMB/analysis/metaG_anvio/Log/%x-%j.err

module load centrifuge/1.0.3-beta

#Set up the path to the working directory
WORKDIR=/proj/DECOMB/analysis/metaG_anvio/
  cd $WORKDIR

anvi-get-sequences-for-gene-calls -c $WORKDIR/05_ANVIO/spades.db -o $WORKDIR/05_ANVIO/spades_gene_calls.fa

centrifuge -f $WORKDIR/05_ANVIO/spades_gene_calls.fa -S $WORKDIR/05_ANVIO/spades_centrifuge_hits.tsv -p 40 \
--report-file $WORKDIR/05_ANVIO/spades_centrifuge_reports.tsv

anvi-import-taxonomy-for-genes -c $WORKDIR/05_ANVIO/spades.db \
-i $WORKDIR/05_ANVIO/spades_centrifuge_reports.tsv $WORKDIR/05_ANVIO/spades_centrifuge_hits.tsv \
-p centrifuge