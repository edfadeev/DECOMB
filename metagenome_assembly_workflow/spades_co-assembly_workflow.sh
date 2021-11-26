#login with port forwarding
ssh -L 5678:localhost:5678 slurm

#workflow for co-assembly of metagenomes using SPAdes and annotation using Anvio
conda activate /apps/anvio/7

#Set up the path to the working directory
WORKDIR=/proj/DECOMB/analysis/metaG_anvio/
cd $WORKDIR

#export sequences from bam files, run QC and trim adapters
sbatch -a 105352-105354 ../DECOMB/metaG/process_raw_sequences.sh

#run co-assembly and produce statistics on it
sbatch ../DECOMB/metaG/co_assembly_spades.sh

#map reads from each library on the co-assembly
sbatch ../DECOMB/metaG/map_reads_bbmap_spades.sh

#import to anvio and annotate the metagenome
sbatch ../DECOMB/metaG/anvio-annotation.sh

################################
#export data for further analysis
################################
#export taxonomy of each gene
anvi-export-table --table genes_taxonomy -o $WORKDIR/05_ANVIO/spades-genes-taxonomy.txt $WORKDIR/05_ANVIO/spades.db

#export taxonomy table
anvi-export-table --table taxon_names -o $WORKDIR/05_ANVIO/spades-tax-names.txt $WORKDIR/05_ANVIO/spades.db

#export AAs of all the gene calls
anvi-get-sequences-for-gene-calls -c $WORKDIR/05_ANVIO/spades.db --report-extended-deflines --get-aa-sequences -o $WORKDIR/05_ANVIO/

sed -n '2,$p' $WORKDIR/05_ANVIO/spades-gene-calls.txt | sort -n -k 1 > $WORKDIR/05_ANVIO/spades-gene-calls-sorted.txt

#generate AAs reference fasta file for the entire metagenome
awk '{print ">"$1"_"$2"\n"$10}' $WORKDIR/05_ANVIO/spades-gene-calls-sorted.txt > $WORKDIR/05_ANVIO/spades-AAs-ref-db.fasta

#export functions of each gene

DATABASES=("COG20_FUNCTION" "Pfam"  "GO" "InterPro" "Hamap")

for db in ${DATABASES[@]}; do
anvi-export-functions -c 05_ANVIO/spades.db --annotation-sources $db -o 05_ANVIO/spades-$db-functions.txt

#remove spaces
sed -e 's/ /_/g' 05_ANVIO/spades-$db-functions.txt > 05_ANVIO/spades-$db-functions-corrected.txt

#sort the table
sed -n '2,$p' 05_ANVIO/spades-$db-functions-corrected.txt | sort -n -k 1 > 05_ANVIO/spades-$db-functions-sorted.txt
done

################################
#Produce automatic bins
################################

mkdir 06_BINS

sbatch ../DECOMB/metaG/bin_co-assembly.sh

################################
#combine bins using DAS Tool
################################
mkdir $WORKDIR/06_BINS/DAS_Tool

#export list of contigs per bin from each binner
anvi-export-collection -C concoct_new \
                        -p $WORKDIR/05_ANVIO/SPAdes/merged_profile/PROFILE.db \
                        -O $WORKDIR/06_BINS/DAS_Tool/spades-concoct

anvi-export-collection -C metabat2_new \
                        -p $WORKDIR/05_ANVIO/SPAdes/merged_profile/PROFILE.db \
                        -O $WORKDIR/06_BINS/DAS_Tool/spades-metabat2
                        
#extract only the contig names without the splits                         
sed 's/_split_[0-9]*//g' $WORKDIR/06_BINS/DAS_Tool/spades-concoct.txt| uniq -u - > $WORKDIR/06_BINS/DAS_Tool/spades_concoct_contig.txt
sed 's/_split_[0-9]*//g' $WORKDIR/06_BINS/DAS_Tool/spades-metabat2.txt| uniq -u - > $WORKDIR/06_BINS/DAS_Tool/spades_metabat2_contig.txt

#generate protein fasta file for DAS Tool
awk '{print ">"$2"_"$1"\n"$10}' $WORKDIR/05_ANVIO/spades-gene-calls-sorted.txt > $WORKDIR/06_BINS/DAS_Tool/spades-AAs-for-binning.fasta

#combine bins using DAS Tool
sbatch ../DECOMB/metaG/combined_binning.sh

#import the bins from DAS
anvi-import-collection --collection-name DAS_Tool \
                        --pan-or-profile-db $WORKDIR/05_ANVIO/SPAdes/merged_profile/PROFILE.db \
                        --contigs-db $WORKDIR/05_ANVIO/spades.db \
                        --contigs-mode $WORKDIR/06_BINS/DAS_Tool/spades__DASTool_scaffolds2bin.txt

#explore bins of CONCOCT (they were better than the ones got out of metabat)
anvi-interactive -p $WORKDIR/05_ANVIO/SPAdes/merged_profile/PROFILE.db -c $WORKDIR/05_ANVIO/spades.db -C DAS_Tool --server-only -P 5678

################################
#Refine successful bins
################################
#refine selected bins

#Bins that were relatively more abundant in the Jelly metagenome

#Bin 2
#The bin was split into two separated bins of Alteromonas
#Bin_2_1 - Alteromonas - length 3.7Mbp (C78.9/R0)
#Bin_2_2 - Alteromonas - length 5.19Mbp (C80.3/R0)
anvi-refine -p $WORKDIR/05_ANVIO/SPAdes/merged_profile/PROFILE.db -c $WORKDIR/05_ANVIO/spades.db -C CONCOCT -b Bin_2 --server-only -P 5678


#Bin 134
#The bin was split into two separated bins of Rhodobacteraceae
#Bin_134_1 - Rhodobacteraceae - length 2.14Mbp (C87.3/R4.2)
#Bin_134_2 - Rhodobacteraceae - length 2.3Mbp (C80.3/R1.4)
anvi-refine -p $WORKDIR/05_ANVIO/SPAdes/merged_profile/PROFILE.db -c $WORKDIR/05_ANVIO/spades.db -C CONCOCT -b Bin_134 --server-only -P 5678

#Bin 5
#Based on coverage large part of the bin was discarded, however the remained bin also varies in coverage
#Bin_5_1 - Alteromonas - length 4.63Mbp (C59.2/R1.4) - not included due to low completness
#Bin_5_2 - family Nitrincolaceae - length 5.19Mbp (C100/R0)
#Bin_5_3 - Reinekea blandensis (Pseudomonadales) - length 4.08Mbp (C100/R2.8)
#Bin_5_4 - Pseudoalteromonas - length 4.87Mbp (C42.3/R0) - not included due to low completness
anvi-refine -p $WORKDIR/05_ANVIO/SPAdes/merged_profile/PROFILE.db -c $WORKDIR/05_ANVIO/spades.db -C CONCOCT -b Bin_5 --server-only -P 5678

#Bin 76
#Based on coverage large part of the bin was discarded.
#Bin_76_1 - Marinobacterium jannaschii (weak taxonomy, but for sure family Nitrincolaceae) - length 3.78Mbp (C100/R0)
anvi-refine -p $WORKDIR/05_ANVIO/SPAdes/merged_profile/PROFILE.db -c $WORKDIR/05_ANVIO/spades.db -C CONCOCT -b Bin_76 --server-only -P 5678

#Bin 84
#in the following step part of the contigs were removed due to coverage differences
#the removed conitgs had no taxonomic assignment
#Bin_84_1 - Pseudoalteromonas phenolica - length 3.92Mbp (C93/R0). 

anvi-refine -p $WORKDIR/05_ANVIO/SPAdes/merged_profile/PROFILE.db -c $WORKDIR/05_ANVIO/spades.db -C CONCOCT -b Bin_84 --server-only -P 5678

#Bin 38
#Based on coverage large part of the bin was discarded, however the remained bin also varies in coverage
#Bin_38_1 - Vibrio - length 1.7Mbp (C83/R0)
anvi-refine -p $WORKDIR/05_ANVIO/SPAdes/merged_profile/PROFILE.db -c $WORKDIR/05_ANVIO/spades.db -C CONCOCT -b Bin_38 --server-only -P 5678

#Bin 102
#Based on coverage large part of the bin was discarded, however the remained bin also varies in coverage
#Bin_102_1 - Kordiimonas lacus - length 2.18Mbp (C73/R0)
anvi-refine -p $WORKDIR/05_ANVIO/SPAdes/merged_profile/PROFILE.db -c $WORKDIR/05_ANVIO/spades.db -C CONCOCT -b Bin_102 --server-only -P 5678

#Merging Bin_75 and Bin_177 did not povide better binning of the Pseudoalteromnas genome
#anvi-merge-bins -p $WORKDIR/05_ANVIO/SPAdes/merged_profile/PROFILE.db -C CONCOCT --bin-names-list "Bin_75, Bin_177" --new-bin-name Bin_75_177
#anvi-refine -p $WORKDIR/05_ANVIO/SPAdes/merged_profile/PROFILE.db -c $WORKDIR/05_ANVIO/spades.db -C CONCOCT -b Bin_75_177 --server-only -P 5678


#Bins abundant in the control metagenome
#Bin 179
#Based on coverage large part of the bin was discarded, however the remained bin also varies in coverage
#Bin_179_1 - family Cellvibrionaceae - length 3.18Mbp (C93/R2.8)
anvi-refine -p $WORKDIR/05_ANVIO/SPAdes/merged_profile/PROFILE.db -c $WORKDIR/05_ANVIO/spades.db -C CONCOCT -b Bin_179 --server-only -P 5678

#Bin 115
#Based on coverage, the bin was split into two sub bins and part of the sequences were discarded.
#Bin_115_1 - Bermanella sp002683575 (Pseudomonadales) - length 4.86Mbp (C98.6/R0)
#Bin_115_2 - Glaciecola sp000155775 (Enterobacterales) - length 2.22Mbp (C100/R0)

anvi-refine -p $WORKDIR/05_ANVIO/SPAdes/merged_profile/PROFILE.db -c $WORKDIR/05_ANVIO/spades.db -C CONCOCT -b Bin_115 --server-only -P 5678

#Bin 12
#Based on coverage large part of the bin was discarded, however the remained bin also varies in coverage
#Bin_12_1 - Saccharospirillum sp003054965 - length 3.75Mbp (C95.8/R1.4)
anvi-refine -p $WORKDIR/05_ANVIO/SPAdes/merged_profile/PROFILE.db -c $WORKDIR/05_ANVIO/spades.db -C CONCOCT -b Bin_102 --server-only -P 5678


#Bins enriched in T0
#Bin 150
#The bin was refined twice to reduce redundancy
#Bin_150_1_1 - Rhodobacteraceae - length 2.6Mbp (C97.2/R1.4)
anvi-refine -p $WORKDIR/05_ANVIO/SPAdes/merged_profile/PROFILE.db -c $WORKDIR/05_ANVIO/spades.db -C CONCOCT -b Bin_150 --server-only -P 5678

#***for not no additional bins were refined for T0***


#summarize and export them for further analysis
mkdir $WORKDIR/06_BINS/REFINED

#export the bins collection
anvi-export-collection -p $WORKDIR/05_ANVIO/SPAdes/merged_profile/PROFILE.db -C CONCOCT -O $WORKDIR/06_BINS/concoct_refined_bins_list

#subset the bins of interest and create new collection
BINS=("Bin_2_1" "Bin_2_2" "Bin_134_1" "Bin_134_2" "Bin_5_2" "Bin_5_3" "Bin_84_1" "Bin_76_1" "Bin_38_1" "Bin_102_1" "Bin_179_1" "Bin_115_2" "Bin_115_1" "Bin_12_1")
for bin in ${BINS[@]}; do 
echo $bin | grep -w -F -f - $WORKDIR/06_BINS/concoct_refined_bins_list.txt >> $WORKDIR/06_BINS/collection_of_interest.txt
done

anvi-import-collection -C Selected_bins \
                        -p $WORKDIR/05_ANVIO/SPAdes/merged_profile/PROFILE.db \
                        -c $WORKDIR/05_ANVIO/spades.db \
                        $WORKDIR/06_BINS/collection_of_interest.txt

#produce summary of the selected bins
sbatch ../DECOMB/metaG/sum_refined_bins.sh


#explore the selected bins collection

#add taxonomy to each bin for visualization (still work in progress)
awk '{print $1,$11,$13}' 06_BINS/Selected_bins_summary/bins_summary.txt > 06_BINS/Selected_bins_tax.txt
sed -i "1s/.*/item_name categorical_1 text_layer_01/" 06_BINS/Selected_bins_tax.txt
sed -i 's/ /\t/g' 06_BINS/Selected_bins_tax.txt

anvi-interactive -p $WORKDIR/05_ANVIO/SPAdes/merged_profile/PROFILE.db -c $WORKDIR/05_ANVIO/spades.db -C Selected_bins --server-only -P 5678


#check bins completness with checkM
mkdir $WORKDIR/06_BINS/Selected_bins_summary/Selected_bins_fasta

find in $WORKDIR/06_BINS/Selected_bins_summary/ -type f -name '*.fa' -exec cp '{}' $WORKDIR/06_BINS/Selected_bins_summary/Selected_bins_fasta

sbatch ../DECOMB/metaG/bin_checkM.sh

################################
#Metabolic reconstruction of each bin
################################
mkdir 07_METABOLISM

#list of bins and their locations
echo -e "name\tcontigs_db_path" > $WORKDIR/07_METABOLISM/selected-bins.txt
for bin in ${BINS[@]}; do 
echo -e "${bin}\t${WORKDIR}/06_BINS/REFINED/${bin}/CONTIGS.db" >> $WORKDIR/07_METABOLISM/selected-bins.txt
done

#populate each bin with KOfam annotation and estimate metabolism
#make sure to set up a KOfam database (run once):
#anvi-setup-kegg-kofams --kegg-data-dir /proj/DECOMB/source/KOfam
sbatch ../DECOMB/metaG/bin_KOfam-annotation.sh

# explore the produces metabolic heatmap
anvi-interactive --manual-mode \
                 -d $WORKDIR/07_METABOLISM/Bins-completeness-MATRIX.txt \
                 -t $WORKDIR/07_METABOLISM/Bins-completeness-MATRIX.txt.newick \
                 -p $WORKDIR/07_METABOLISM/Bins_metabolism_PROFILE.db \
                 --title "Bins Metabolism Heatmap" \
                 --server-only -P 5678

#produce table of the different groups of bins
JELLY_BINS=("Bin_2_1" "Bin_2_2" "Bin_134_1" "Bin_134_2" "Bin_5_2" "Bin_5_3" "Bin_84_1" "Bin_76_1" "Bin_38_1" "Bin_102_1") 

CONTROL_BINS=("Bin_179_1" "Bin_115_2" "Bin_115_1" "Bin_12_1")

echo -e "sample\tgroup" > $WORKDIR/07_METABOLISM/bin-groups.txt
for bin in ${JELLY_BINS[@]}; do 
echo -e "${bin}\tJelly" >> $WORKDIR/07_METABOLISM/bin-groups.txt
done

for bin in ${CONTROL_BINS[@]}; do 
echo -e "${bin}\tControl" >> $WORKDIR/07_METABOLISM/bin-groups.txt
done

#load R version 4.0.2
module unload R
module load R/4.0.2

#calculate enriched metabolism
anvi-compute-functional-enrichment -M $WORKDIR/07_METABOLISM/Bins_modules.txt \
                                   -G $WORKDIR/07_METABOLISM/bin-groups.txt \
                                   -o $WORKDIR/07_METABOLISM/Bin_enriched_modules.txt
