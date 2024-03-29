#login with port forwarding
ssh -L 5678:localhost:5678 slurm

#load conda
module load conda

#workflow for co-assembly of metagenomes using SPAdes and annotation using Anvio
conda activate anvio-dev

#Set up the path to the working directory and the scripts directory
DECOMB_git=/scratch/oceanography/efadeev/DECOMB/DECOMB_git

WORKDIR=/scratch/oceanography/efadeev/DECOMB/analysis/metaG_anvio
# currently at WORKDIR=/home/project/oceanography/DECOMB/analysis

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

#export all the gene calls
anvi-export-gene-calls -c $WORKDIR/05_ANVIO/spades.db --gene-caller prodigal -o $WORKDIR/05_ANVIO/spades-gene-calls.txt

#generate AAs reference fasta file for the entire metagenome
awk '{print ">"$1"_"$2"\n"$10}' $WORKDIR/05_ANVIO/spades-gene-calls.txt > $WORKDIR/05_ANVIO/spades-AAs-ref-db.fasta

#export functions of each gene
DATABASES=("COG20_FUNCTION" "Pfam"  "GO" "InterPro" "Hamap" "KEGG_Class" "KeggGhostKoala" "KOfam" "COG20_CATEGORY" "KEGG_Module" "COG20_PATHWAY")

for db in ${DATABASES[@]}; do
anvi-export-functions -c $WORKDIR/05_ANVIO/spades.db --annotation-sources $db -o $WORKDIR/05_ANVIO/spades-$db-functions.txt
#remove spaces
sed -i -e 's/ /_/g' $WORKDIR/05_ANVIO/spades-$db-functions.txt 

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
-b Bin_50_sub --server-only -P 5678
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
awk '{print $1,$10,$11,$12}' $WORKDIR/06_BINS/Refined_DAS_bins/bins_summary.txt > $WORKDIR/06_BINS/Refined_DAS_bins_tax.txt
sed -i "1s/.*/item_name Class Order Family/" $WORKDIR/06_BINS/Refined_DAS_bins_tax.txt
sed -i 's/ /\t/g' $WORKDIR/06_BINS/Refined_DAS_bins_tax.txt

#explore the selected bins collection
anvi-interactive -p $WORKDIR/05_ANVIO/SPAdes/merged_profile/PROFILE.db -c $WORKDIR/05_ANVIO/spades.db \
-C Refined_DAS_bins --server-only -P 5678 --additional-layers $WORKDIR/06_BINS/Refined_DAS_bins_tax.txt

#summarize and export the bins for further analysis
sbatch $DECOMB_git/metagenome_assembly_workflow/sum_refined_bins.sh

#generate table with gene calls per bin
readarray -t BINS < $WORKDIR/06_BINS/Refined_bins.txt

echo -e "Bin\tcontig\tgene_callers_id" > $WORKDIR/06_BINS/selected-bins-gene-calls.txt
for bin in ${BINS[@]}; do 
sed '1d' $WORKDIR/06_BINS/REFINED/$bin/$bin-gene-calls.txt | awk -v b=$bin '{print b"\t"$2"\t"$1}' >> $WORKDIR/06_BINS/selected-bins-gene-calls.txt
done

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



################################
#Generate phylogeny and pangenome for Pseudoalteromonas bin
################################
mkdir 08_BIN_PAN

#download all the complete pseudoalteromonas genomes
ncbi-genome-download bacteria \
                     --assembly-level chromosome,complete \
                     --genera Pseudoalteromonas \
                     --metadata $WORKDIR/08_BIN_PAN/NCBI-METADATA.txt

anvi-script-process-genbank-metadata -m $WORKDIR/08_BIN_PAN/NCBI-METADATA.txt \
                                     --output-dir $WORKDIR/08_BIN_PAN/Pseudoalt_genomes \
                                     --output-fasta-txt $WORKDIR/08_BIN_PAN/Pseudoalt_genomes.txt

#extract names and paths of the fastas
cut -f1,2 08_BIN_PAN/Pseudoalt_genomes.txt> 08_BIN_PAN/Pseudoalt_fasta.txt

#add manually the bins
cp 06_BINS/Refined_DAS_bins/bin_by_bin/Bin_84/Bin_84-contigs.fa 08_BIN_PAN/Pseudoalt_genomes/Bin_84.fa

anvi-run-workflow -w contigs --get-default-config 08_BIN_PAN/contig-config.json

anvi-run-workflow -w contigs \
-c $WORKDIR/08_BIN_PAN/contig-config.json \
--additional-params \
--cores 48 \
--cluster \
'sbatch --job-name=Pseudoalt_genomes2anvio \
        --mail-user=eduard.fadeev@univie.ac.at \
        --output=/scratch/oceanography/efadeev/DECOMB/analysis/metaG_anvio/Log/%x-%j.out \
        --error=/scratch/oceanography/efadeev/DECOMB/analysis/metaG_anvio/Log/%x-%j.out \
        --cpus-per-task=6 \
        --time=1-12:00:00 \
        --mem=20GB'

anvi-run-workflow -w phylogenomics --get-default-config 08_BIN_PAN/phylo-config.json

anvi-run-workflow -w phylogenomics \
-c $WORKDIR/08_BIN_PAN/phylo-config.json \
--additional-params \
--cores 10 \
--cluster \
'sbatch --job-name=Pseudoalt_phylogenomics \
        --mail-user=eduard.fadeev@univie.ac.at \
        --output=/scratch/oceanography/efadeev/DECOMB/analysis/metaG_anvio/Log/%x-%j.out \
        --error=/scratch/oceanography/efadeev/DECOMB/analysis/metaG_anvio/Log/%x-%j.out \
        --cpus-per-task=10 \
        --time=1-12:00:00 \
        --mem=40GB'

#inspect the tree
anvi-interactive --tree 08_BIN_PAN/PHYLOGENOMICS/Pseudoalteromonas_phylogenomics-proteins_GAPS_REMOVED.fa.contree \
-p 08_BIN_PAN/phylo-profile.db --manual --server-only -P 5678

################################
#Generate pangenome for P. phenolica
################################
#download all genomes of P. phenolica
ncbi-genome-download bacteria -t 161398,1315281 \
            --assembly-level all \
            --metadata $WORKDIR/08_BIN_PAN/Pphenol-METADATA.txt

anvi-script-process-genbank-metadata -m $WORKDIR/08_BIN_PAN/Pphenol-METADATA.txt \
                                     --output-dir $WORKDIR/08_BIN_PAN/Pphenol_genomes \
                                     --output-fasta-txt $WORKDIR/08_BIN_PAN/Pphenol_genomes.txt

#remove the annotation paths
cut -f1,2 $WORKDIR/08_BIN_PAN/Pphenol_genomes.txt> $WORKDIR/08_BIN_PAN/Pphenol_fasta.txt

#add the paths of the bins
tail -n 3 $WORKDIR/08_BIN_PAN/Pseudoalt_fasta.txt >> $WORKDIR/08_BIN_PAN/Pphenol_fasta.txt

#need to move the entire analysis to /tmp/
#otherwise the makeblastdb is not working
#run the pangenome workflow
screen -S anvio_snakemake

WORKDIR=/scratch/oceanography/efadeev/DECOMB/analysis/Pphenol_pangenomics/

anvi-run-workflow -w pangenomics \
-c Pphenol_pangenome.json \
--additional-params \
--jobs 10 \
--cluster \
'sbatch --job-name=Pphenol_pangenomics \
        --mail-user=ALL \
        --output=/scratch/oceanography/efadeev/DECOMB/analysis/Pphenol_pangenomics/00_LOGS/%x-%j.out \
        --time=1-24:00:00 \
        --mem=50GB'

##!!!!#######
#due to some unclear configuration of the cluster
#the makeblastdb command needs to be executed separately with output on tmp

#inspect the pangenome
anvi-display-pan -g $WORKDIR/03_PAN/P_phenolica_pangenome-GENOMES.db \
-p $WORKDIR/03_PAN/P_phenolica_pangenome-PAN.db --server-only -P 5678