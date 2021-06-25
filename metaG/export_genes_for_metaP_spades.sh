#!/bin/bash
#
#SBATCH --job-name=spades_assembly_genes_output
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

#annotate COGs
anvi-run-ncbi-cogs -c $WORKDIR/05_ANVIO/spades.db \
--num-threads 20 --cog-data-dir /proj/DECOMB/source/anvio-COG

anvi-export-functions -c 05_ANVIO/spades.db --annotation-sources COG20_FUNCTION -o 05_ANVIO/spades-functions.txt




#add gene taxonomy
anvi-get-sequences-for-gene-calls -c $WORKDIR/05_ANVIO/spades.db -o $WORKDIR/05_ANVIO/spades_gene_calls.fa

centrifuge -f $WORKDIR/05_ANVIO/spades_gene_calls.fa -S $WORKDIR/05_ANVIO/spades_centrifuge_hits.tsv -p 20 \
--report-file $WORKDIR/05_ANVIO/spades_centrifuge_reports.tsv

anvi-import-taxonomy-for-genes -c $WORKDIR/05_ANVIO/spades.db \
-i $WORKDIR/05_ANVIO/spades_centrifuge_reports.tsv $WORKDIR/05_ANVIO/spades_centrifuge_hits.tsv \
-p centrifuge

#export taxonomy of each gene
anvi-export-table --table genes_taxonomy -o spades-genes-taxonomy.txt spades.db
#export taxonomy table
anvi-export-table --table taxon_names -o spades-tax-names.txt spades.db



#generate reference database for metaP and functions table
sed -n '2,$p' 05_ANVIO/spades-functions.txt | sort > spades-functions-sorted.txt
sed -n '2,$p' 05_ANVIO/spades-gene-calls.txt | sort > spades-gene-calls-sorted.txt

#merge all data together
sed -e 's/ /_/g' 05_ANVIO/spades-functions.txt > 05_ANVIO/spades-functions-corrected.txt
join -11 -a1 --header 05_ANVIO/spades-gene-calls.txt 05_ANVIO/spades-functions-corrected.txt > 05_ANVIO/spades-genes-fun-merged.txt

#generate AAs fasta file
awk '{print ">"$2"\n"$10}' 05_ANVIO/spades-genes-fun-merged.txt > 05_ANVIO/spades-AA-ref-db.fasta
