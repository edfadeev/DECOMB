#!/bin/bash
#
#SBATCH --job-name=annotate_assembly
#SBATCH --cpus-per-task=20
#SBATCH --mem=20GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --output=/proj/DECOMB/analysis/metaG_anvio/Log/%x-%j.out
#SBATCH --error=/proj/DECOMB/analysis/metaG_anvio/Log/%x-%j.err

#load module
module load centrifuge/1.0.3-beta

#Set up the path to the working directory
WORKDIR=/proj/DECOMB/analysis/metaG_anvio/
cd $WORKDIR

#generate db
anvi-gen-contigs-database --num-threads 20 -f $WORKDIR/03_CONTIGS/spades-contigs-fixed.fasta \
-o $WORKDIR/05_ANVIO/spades.db --project-name spades_assembly

#annotate using HMMs
anvi-run-hmms -c $WORKDIR/05_ANVIO/spades.db --num-threads 20

#annotate genes with COGs
anvi-run-ncbi-cogs -c $WORKDIR/05_ANVIO/spades.db \
--num-threads 20 --cog-data-dir /proj/DECOMB/source/anvio-COG

#add annotations from interproscan
anvi-get-sequences-for-gene-calls -c $WORKDIR/05_ANVIO/spades.db --get-aa-sequences -o $WORKDIR/05_ANVIO/spades-AA-sequences.fa

/mirror/interpro/interproscan-5.51-85.0/interproscan.sh -i $WORKDIR/05_ANVIO/spades-AA-sequences.fa -f tsv --cpu 20 -o $WORKDIR/05_ANVIO/spades-interpro-output.tsv

/proj/DECOMB/source/InterProScanParser/iprs2anvio.sh -i $WORKDIR/05_ANVIO/spades-interpro-output.tsv -o $WORKDIR/05_ANVIO/spades-interpro-output -g -p -r

anvi-import-functions -c $WORKDIR/05_ANVIO/spades.db -i $WORKDIR/05_ANVIO/spades-interpro-output_iprs2anvio.tsv

#add gene taxonomy
anvi-get-sequences-for-gene-calls -c $WORKDIR/05_ANVIO/spades.db -o $WORKDIR/05_ANVIO/spades_gene_calls.fa

centrifuge -f $WORKDIR/05_ANVIO/spades_gene_calls.fa -S $WORKDIR/05_ANVIO/spades_centrifuge_hits.tsv -p 20 \
--report-file $WORKDIR/05_ANVIO/spades_centrifuge_reports.tsv

anvi-import-taxonomy-for-genes -c $WORKDIR/05_ANVIO/spades.db \
-i $WORKDIR/05_ANVIO/spades_centrifuge_reports.tsv $WORKDIR/05_ANVIO/spades_centrifuge_hits.tsv \
-p centrifuge
