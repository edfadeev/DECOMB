#login with port forwarding
ssh -L 5678:localhost:5678 slurm

#load conda
module load conda

#workflow for co-assembly of metagenomes using SPAdes and annotation using Anvio
conda activate anvio_7.1

#Set up the path to the working directory and the scripts directory
DECOMB_git=/scratch/oceanography/efadeev/20220302_DECOMB/DECOMB_git

WORKDIR=/home/project/oceanography/DECOMB/metaG_anvio
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
sbatch ../DECOMB/metaG/DAS_Tool_binning.sh

#import the bins from DAS
anvi-import-collection --collection-name DAS_Tool \
                        --pan-or-profile-db $WORKDIR/05_ANVIO/SPAdes/merged_profile/PROFILE.db \
                        --contigs-db $WORKDIR/05_ANVIO/spades.db \
                        --contigs-mode $WORKDIR/06_BINS/DAS_Tool/spades__DASTool_scaffolds2bin.txt

#explore bins of DAS Tool
anvi-interactive -p $WORKDIR/05_ANVIO/SPAdes/merged_profile/PROFILE.db -c $WORKDIR/05_ANVIO/spades.db -C DAS_Tool --server-only -P 5678

#summarize the DAS Tool bins
anvi-summarize -p $WORKDIR/05_ANVIO/SPAdes/merged_profile/PROFILE.db -c $WORKDIR/05_ANVIO/spades.db -C DAS_Tool \
-o $WORKDIR/06_BINS/DAS_Tool_summary


################################
#Refine successful bins
################################
#refine selected bins
#Bin 189
anvi-refine -p $WORKDIR/05_ANVIO/SPAdes/merged_profile/PROFILE.db -c $WORKDIR/05_ANVIO/spades.db -C DAS_Tool \
-b Bin_189 --server-only -P 5678
#Refined into - UBA8296 sp002338335 (Balneolales) - length 1.4Mbp (C94.4/R8.5)

#Bin 50_sub
anvi-refine -p $WORKDIR/05_ANVIO/SPAdes/merged_profile/PROFILE.db -c $WORKDIR/05_ANVIO/spades.db -C DAS_Tool \
-b Bin_189 --server-only -P 5678
#Refined into Bin_50_sub_1 UBA7446 sp002470745 (Flavobacteriales) - length 1.09Mbp (C71.8/R8.5)

#Bin METABAT_111
anvi-refine -p $WORKDIR/05_ANVIO/SPAdes/merged_profile/PROFILE.db -c $WORKDIR/05_ANVIO/spades.db -C DAS_Tool \
-b METABAT_111 --server-only -P 5678
#Refined into METABAT_111_1 GCA-2707915 sp004214065 (SAR86) - length 1.2Mbp (C63.4/R4.2)

#bin smaller than 1Mbp with completness below 80% were excluded, as well as bins with redundancy of >10% (according to anvio)
#in total 44 bins were selected
#store the list of bins in '$WORKDIR/06_BINS/Refined_bins.txt'
#Create a collection of the manually selected bins
anvi-export-collection -C DAS_Tool \
                        -p $WORKDIR/05_ANVIO/SPAdes/merged_profile/PROFILE.db \
                        -O $WORKDIR/06_BINS/DAS_Tool/DAS_Tool

sed 's/_split_[0-9]*//g' $WORKDIR/06_BINS/DAS_Tool/DAS_Tool.txt| uniq -u - > $WORKDIR/06_BINS/DAS_Tool/DAS_Tool_contig.txt
awk 'NR==FNR{a[$0];next} $NF in a' 06_BINS/Refined_bins.txt 06_BINS/DAS_Tool/DAS_Tool_contig.txt > 06_BINS/Refined_bins_collection.txt

#import the collection to Anvio
anvi-import-collection --collection-name Refined_DAS_bins \
                        --pan-or-profile-db $WORKDIR/05_ANVIO/SPAdes/merged_profile/PROFILE.db \
                        --contigs-db $WORKDIR/05_ANVIO/spades.db --contigs-mode \
                        $WORKDIR/06_BINS/Refined_bins_collection.txt
                        

#add taxonomy to each bin for visualization
awk '{print $1,$11}' $WORKDIR/06_BINS/Refined_DAS_bins/bins_summary.txt > $WORKDIR/06_BINS/Refined_DAS_bins_tax.txt
sed -i "1s/.*/item_name categorical_1/" $WORKDIR/06_BINS/Refined_DAS_bins_tax.txt
sed -i 's/ /\t/g' $WORKDIR/06_BINS/Refined_DAS_bins_tax.txt

#explore the selected bins collection
anvi-interactive -p $WORKDIR/05_ANVIO/SPAdes/merged_profile/PROFILE.db -c $WORKDIR/05_ANVIO/spades.db \
-C Refined_DAS_bins --server-only -P 5678 --additional-layers $WORKDIR/06_BINS/Refined_DAS_bins_tax.txt

#summarize and export the bins for further analysis
sbatch $DECOMB_git/metagenome_assembly_workflow/sum_refined_bins.sh

#the gff3 and fasta files were using for generating genbank files and annotation using RAST server

#check bins completness with checkM using the fasta files
sbatch $DECOMB_git/metagenome_assembly_workflow/bins_checkm.sh


################################
#Metabolic reconstruction of each bin
################################
mkdir 07_METABOLISM

#list of bins and their locations
#list of selected bins was manually produced and saved into:
# 'selected-bins.csv'
readarray -t BINS < $WORKDIR/06_BINS/Refined_bins.txt

echo -e "name\tcontigs_db_path" > $WORKDIR/07_METABOLISM/selected-bins-collection.txt
for bin in ${BINS[@]}; do 
echo -e "${bin}\t${WORKDIR}/06_BINS/REFINED/${bin}/CONTIGS.db" >> $WORKDIR/07_METABOLISM/selected-bins-collection.txt
done

#populate each bin with KOfam annotation and estimate metabolism
#make sure to set up a KOfam database (run once):
#anvi-setup-kegg-kofams --kegg-data-dir /proj/DECOMB/source/KOfam
sbatch $DECOMB_git/metagenome_assembly_workflow/Bin_KEGG_modules.sh

# explore the produced metabolic heatmap
anvi-interactive --manual-mode \
                 -d $WORKDIR/07_METABOLISM/Bins-completeness-MATRIX.txt \
                 -t $WORKDIR/07_METABOLISM/Bins-completeness-MATRIX.txt.newick \
                 -p $WORKDIR/07_METABOLISM/Bins_metabolism_PROFILE.db \
                 --title "Bins Metabolism Heatmap" \
                 --server-only -P 5678

#############################################################################
# draft
#############################################################################
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








