#!/bin/bash
#
#SBATCH --job-name=summarize_refined_bins
#SBATCH --cpus-per-task=10
#SBATCH --mem=40GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --output=/scratch/oceanography/efadeev/DECOMB/analysis/metaG_anvio/Log/%x-%j.out
#SBATCH --error=/scratch/oceanography/efadeev/DECOMB/analysis/metaG_anvio/Log/%x-%j.err

#Set up the path to the working directory
WORKDIR=/scratch/oceanography/efadeev/DECOMB/analysis/metaG_anvio
cd $WORKDIR

#summarize refined bins
anvi-summarize -p $WORKDIR/05_ANVIO/SPAdes/merged_profile/PROFILE.db -c $WORKDIR/05_ANVIO/spades.db -C Refined_DAS_bins \
-o $WORKDIR/06_BINS/Refined_DAS_bins

#export each bin into a separate db
anvi-split -p $WORKDIR/05_ANVIO/SPAdes/merged_profile/PROFILE.db -c $WORKDIR/05_ANVIO/spades.db \
--collection-name Refined_DAS_bins -o $WORKDIR/06_BINS/REFINED

#extract gene calls as GFF3 files and contig sequences 
readarray -t BINS < $WORKDIR/06_BINS/Refined_bins.txt
for bin in ${BINS[@]}; do 
anvi-get-sequences-for-gene-calls --contigs-db $WORKDIR/06_BINS/REFINED/${bin}/CONTIGS.db \
--export-gff3 --output-file $WORKDIR/06_BINS/REFINED/${bin}.gff

anvi-export-contigs --contigs-db $WORKDIR/06_BINS/REFINED/${bin}/CONTIGS.db --output-file $WORKDIR/06_BINS/REFINED/${bin}.fasta

done


