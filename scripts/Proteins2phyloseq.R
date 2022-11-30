########################################
#This script merges the different gene annotations that were produced by Anvio.
#Then it combines the metaproteome and exoproteome tables with the annotation and saves is as a phyloseq object
########################################

#load libraries
require(dplyr)
require(tidyr)
require(phyloseq)

########################################
#import taxonomy and annotation of each gene in the reference metagenome
########################################
#genes list
gene_annotations_df <- read.csv("./data/metagenome/spades-gene-calls.txt",
                                sep="\t", h= T) %>% 
                        left_join(read.csv("./data/metagenome/Refined-bins-gene-calls.txt",sep="\t", h= T), 
                                    by = c("gene_callers_id","contig")) 

#import gene taxonomy
gene_taxa_df<- read.csv("./data/metagenome/spades-genes-taxonomy.txt",
                        sep="\t", h= T)  %>% 
                left_join(read.csv("./data/metagenome/spades-tax-names.txt",sep="\t", h= T), 
                          by = "taxon_id") %>% 
                dplyr::rename(Class = t_order ,Order = t_class,
                              Phylum = t_phylum, Family = t_family,
                              Genus = t_genus) #switch columns due to some taxonomy parsing bug in anvio


#import functional annotations
source_files <- list.files(path = "./data/metagenome", 
                           pattern = "-.*functions\\.txt",
                           full.names = TRUE)

annotations_list <- lapply(source_files, function(x) {
                            source_db <- unlist(strsplit(basename(x), split ="-"))[2]
                            dat <- read.csv(x,sep="\t",col.names = c("gene_callers_id","db",
                                                          "accession","function.","e_value")) %>% 
    #mutate_all(as.character) %>% 
                                   select(gene_callers_id, accession, function.) %>% 
                                   group_by(gene_callers_id)%>%
                                   summarise_each(funs(paste(unique(.), collapse='|')), matches('^\\D+$')) %>% 
                                   plyr::rename(replace= c("accession"=paste(source_db, "accession", sep ="_"), 
                                                            "function."=paste(source_db, "function", sep ="_"))) 
                            return(dat)
                            } )

#merge all dataframes
source_dbs <- gsub("spades-|-functions.txt","",basename(source_files))

annotations_df <- annotations_list %>% 
  purrr::reduce(full_join, by = "gene_callers_id") %>% 
  select(gene_callers_id,contains(source_dbs)) %>% 
  replace(is.na(.), "Unk")

gene_meta_df <- gene_annotations_df %>% 
  select(gene_callers_id, contig, start, stop, partial, aa_sequence, Bin) %>% 
  merge(annotations_df, by = "gene_callers_id") %>% 
  merge(gene_taxa_df, by = "gene_callers_id") %>% 
  mutate(prot_length = nchar(aa_sequence)) %>% 
  mutate(rows = gene_callers_id) %>% 
  tibble::column_to_rownames('rows')

########################################
#import samples metadata
########################################
#run in terminal to remove paths of the raw files
#sed "s/F:\\\\EFadeev\\\\DECOMB_raw_data\\\\DECOMB-all\\\\//g"
#DE-COMB_all_prot_InputFiles.txt > DE-COMB_all_prot_InputFiles_corrected.txt 

samples_df <- read.csv("./data/metaproteome/DE-COMB_all_prot_InputFiles_corrected.txt",
                                sep="\t", h= T) %>% 
                mutate(File.Name = case_when(File.Name =="T0_EL_200616145603.raw" ~ "T0_EL (2).raw",
                                             TRUE ~ File.Name)) %>% 
                mutate(Sample_name= gsub(" \\(2\\)","_2",gsub("\\.raw","",File.Name)),
                       Treatment = case_when(gsub("[1-3]_.*","",File.Name)=="C" ~ "Control",
                                             gsub("[1-3]_.*","",File.Name)=="J" ~ "Jelly",
                                             gsub("0_.*","",File.Name)=="T" ~ "Inoculum"),
                       Fraction = gsub(" ","", gsub("\\(.*|\\.raw","",gsub("C[1-3]_|J[1-3]_|T0_","",File.Name))),
                       Replicate = gsub("C|J","",gsub("_.*","",File.Name)),
                       Run = case_when(grepl("_2", Sample_name) ~ "B",
                                              TRUE ~ "A"),
                       SampleID = gsub("_2","",Sample_name)) %>% 
                select(File.ID,SampleID,Sample_name,Treatment,Fraction,Replicate,Run)

 

########################################
#import metaproteome data
########################################
#rename sample columns and filter out low confidence proteins and replace NA with 0
PD_df <- read.csv("data/metaproteome/DE-COMB_all_prot_Proteins.txt",
                  sep="\t", h= T)%>% 
  dplyr::rename(gene_callers_id = Accession) %>% 
  mutate(gene_callers_id = as.integer(gsub("_.*","",gene_callers_id))) %>% 
  select_at(vars(!contains("Found"))) %>% 
  rename_with(~gsub("Abundance\\.|\\.Sample","",.), everything()) %>% 
  rename_at(all_of(samples_df$File.ID), ~ samples_df$Sample_name) %>% 
  filter(Number.of.PSMs >=2 , Number.of.Unique.Peptides>=1)%>% 
  mutate_if(is.numeric, funs(replace_na(., 0))) %>% 
  mutate_if(is.numeric,as.integer) %>% 
  select(c("gene_callers_id",samples_df$Sample_name))

########################################
#for convenient analysis the metaproteome data is imported into phyloseq
########################################
#metadata table
meta_df<- sample_data(samples_df)
sample_names(meta_df)<-samples_df$Sample_name     

#produce protein counts table
metaP_prot_counts<- otu_table(data.frame(PD_df[, all_of(samples_df$Sample_name)], 
                                         row.names = PD_df$gene_callers_id), taxa_are_rows = TRUE)

#gene ids as taxonomy table
protein_annotation<- tax_table(as.matrix(gene_meta_df))                          
                                
#merge into phyloseq object
metaP_obj0<- phyloseq(metaP_prot_counts, protein_annotation, meta_df)

#remove proteins that were not observed 
metaP_obj0<- prune_taxa(taxa_sums(metaP_obj0)>0,metaP_obj0)

#save metaproteome phyloseq
saveRDS(metaP_obj0, "data/metaproteome/metaP_ps_raw.rds")

#print session info and clean the workspace
sessionInfo()
rm(list = ls())
gc()
