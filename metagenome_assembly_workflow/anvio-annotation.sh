#!/bin/bash
#
#SBATCH --job-name=annotate_assembly
#SBATCH --cpus-per-task=40
#SBATCH --mem=100GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --output=/proj/DECOMB/analysis/metaG_anvio/Log/%x-%j.out
#SBATCH --error=/proj/DECOMB/analysis/metaG_anvio/Log/%x-%j.err

#load module
module load kaiju/1.7.3

#Set up the path to the working directory
WORKDIR=/proj/DECOMB/analysis/metaG_anvio/
cd $WORKDIR

#generate db
anvi-gen-contigs-database --num-threads 40 -f $WORKDIR/03_CONTIGS/spades-contigs-fixed.fasta \
-o $WORKDIR/05_ANVIO/spades.db --project-name spades_assembly

#annotate using HMMs
anvi-run-hmms -c $WORKDIR/05_ANVIO/spades.db --num-threads 40

#annotate genes with COGs
anvi-run-ncbi-cogs -c $WORKDIR/05_ANVIO/spades.db \
--num-threads 40 --cog-data-dir /proj/DECOMB/source/anvio-COG

#add annotations from interproscan
#anvi-get-sequences-for-gene-calls -c $WORKDIR/05_ANVIO/spades.db --get-aa-sequences -o $WORKDIR/05_ANVIO/spades-AA-sequences.fa

/proj/DECOMB/source/my_interproscan/interproscan-5.52-86.0/interproscan.sh -i $WORKDIR/05_ANVIO/spades-AA-sequences.fa -f tsv --cpu 40 --goterms --pathways --applications Hamap,SUPERFAMILY,TIGRFAM,SFLD,Pfam,PANTHER -o $WORKDIR/05_ANVIO/spades-interpro-output.tsv -dp 

/proj/DECOMB/source/InterProScanParser/iprs2anvio.sh -i $WORKDIR/05_ANVIO/spades-interpro-output.tsv -o $WORKDIR/05_ANVIO/spades-interpro-output -g -p -r

anvi-import-functions -c $WORKDIR/05_ANVIO/spades.db -i $WORKDIR/05_ANVIO/spades-interpro-output_iprs2anvio.tsv


#add annotations for KEGG modules
anvi-run-kegg-kofams -c $WORKDIR/05_ANVIO/spades.db \
--kegg-data-dir /proj/DECOMB/source/KOfam -T 40 \
--just-do-it

#produce tables with KO hits and modules
anvi-estimate-metabolism -c $WORKDIR/05_ANVIO/spades.db \
-O $WORKDIR/05_ANVIO/spades-Kofam \
--metagenome-mode \
--kegg-data-dir /proj/DECOMB/source/KOfam \
--kegg-output-modes kofam_hits,modules

#add gene taxonomic assignments to genes
kaiju -t /apps/kaiju/database/nodes.dmp \
-f /apps/kaiju/database/refseq/kaiju_db_refseq.fmi \
-i $WORKDIR/05_ANVIO/spades_gene_calls.fa  \
-o $WORKDIR/05_ANVIO/spades_gene_calls_refseq_tax.out \
-z 40 \
-v

kaiju-addTaxonNames -t /apps/kaiju/database/nodes.dmp \
-n /apps/kaiju/database/names.dmp \
-i $WORKDIR/05_ANVIO/spades_gene_calls_refseq_tax.out \
-o $WORKDIR/05_ANVIO/spades_gene_calls_refseq_names.out \
-r superkingdom,phylum,order,class,family,genus,species

anvi-import-taxonomy-for-genes -i $WORKDIR/05_ANVIO/spades_gene_calls_refseq_names.out \
-c $WORKDIR/05_ANVIO/spades.db \
-p kaiju --just-do-it

