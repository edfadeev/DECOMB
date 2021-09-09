########################################
#This script merges the different gene annotations that were produced by Anvio.
#Then it combines the metaproteome and exoproteome tables with the annotation and saves is as a phyloseq object
########################################

#set working directory
#macOS
wd <- "/Users/eduardfadeev/Google Drive (dr.eduard.fadeev@gmail.com)/DECOMB/"

#Linux 
wd <- "~/Data/Postdoc-Vienna/DECOMB/"

#Windows
wd <- "D:/Postdoc-Vienna/DECOMB/"

#load libraries
require(dplyr)
require(tidyr)
require(phyloseq)

########################################
#import taxonomy and annotation of each gene in the reference metagenome
########################################
metaG_annotations<- read.csv("data/metaG_full_annotations.txt",
                             row.names = 1)

#filter out proteins that had no taxonomic or functional annotation
metaG_annotations <- metaG_annotations %>%
  select(taxon_id, COG20_FUNCTION_accession, GO_accession, 
         Pfam_accession,InterPro_accession, Hamap_accession) %>% 
  filter(if_any(everything(), ~ !is.na(.)))

########################################
#import metaproteome data
########################################
#generate sample list of metaproteomes
metaP_samples <- data.frame(Sample_name = c("C1_MP","C2_MP","C3_MP","J1_MP","J2_MP","J3_MP","T0_MP"),
                            Sample_PD_IDs= c("F2","F4","F6","F8","F10","F13","F15"),
                            Type = c(rep("Control",3),rep("Jelly",3),"T0"),
                            Replicate = c(1:3, 1:3, 1), 
                            Fraction = "MetaP",
                            row.names = c("C1_MP","C2_MP","C3_MP","J1_MP","J2_MP","J3_MP","T0_MP"))
#import metaP data
metaP_raw <- read.csv(paste(wd,"metaP_analysis/metaP-fractions/DECOMB-frac-metaP-concensus_Proteins.txt", sep=""),sep="\t", h= T)

#rename sample columns and filter out low confidence proteins and replace NA with 0
metaP_filt <- metaP_raw %>% 
  dplyr::rename(gene_caller_id = Accession) %>% 
  select_at(vars(!contains("Found"))) %>% 
  rename_with(~gsub("Abundance\\.|\\.Sample","",.), everything()) %>% 
  rename_at(all_of(metaP_samples$Sample_PD_IDs), ~ metaP_samples$Sample_name) %>% 
  filter(Number.of.PSMs >=2 , Number.of.Unique.Peptides>=1)%>% 
  mutate_if(is.numeric, funs(replace_na(., 0))) %>% 
  mutate_if(is.numeric,as.integer)


########################################
#for convenient analysis the metaproteome data is imported into phyloseq
########################################
#produce protein counts table
metaP_prot_counts<- metaP_filt %>% select(c("gene_caller_id", metaP_samples$Sample_name))
metaP_prot_counts<- otu_table(data.frame(metaP_prot_counts[, all_of(metaP_samples$Sample_name)], row.names = metaP_prot_counts$gene_caller_id), taxa_are_rows = TRUE)

#gene ids as taxonomy table
metaP_annotation<- tax_table(as.matrix(data.frame(metaG_annotations, row.names = metaG_annotations$gene_caller_id)))

#merge into phyloseq object
metaP_obj0<- phyloseq(metaP_prot_counts, metaP_annotation, sample_data(metaP_samples))

#remove proteins that were not observed 
metaP_obj0<- prune_taxa(taxa_sums(metaP_obj0)>0,metaP_obj0)

#save metaproteome phyloseq
saveRDS(metaP_obj0, "data/metaP_ps_raw.rds")

########################################
#import exoproteome data
########################################
#generate sample list of exoproteomes
exoP_sample <- data.frame(Sample_name = c("C1_exoP","C2_exoP","C3_exoP","J1_exoP","J2_exoP","J3_exoP","T0_exoP"),
                          Sample_PD_IDs = c("F1","F2","F3","F4","F5","F6","F7"),
                          Type = c(rep("Control",3),rep("Jelly",3),"T0"),
                          Replicate = c(1:3, 1:3, 1), 
                          Fraction = "ExoP",
                          row.names = c("C1_exoP","C2_exoP","C3_exoP","J1_exoP","J2_exoP","J3_exoP","T0_exoP"))
#import ExoP data
exoP_raw <- read.csv(paste(wd,"metaP_analysis/exoP-fractions/DECOMB-frac-exoP-concensus_Proteins.txt", sep=""),sep="\t", h= T)

#rename sample columns and filter out low confidence proteins and replace NA with 0
exoP_filt <- exoP_raw %>% 
  dplyr::rename(gene_caller_id = Accession) %>% 
  select_at(vars(!contains("Found"))) %>% 
  rename_with(~gsub("Abundance\\.|\\.Sample","",.), everything()) %>% 
  rename_at(all_of(exoP_sample$Sample_PD_IDs), ~ exoP_sample$Sample_name) %>% 
  filter(Number.of.PSMs >=2 , Number.of.Unique.Peptides>=1)%>% 
  mutate_if(is.numeric, funs(replace_na(., 0))) %>% 
  mutate_if(is.numeric,as.integer)


########################################
#merge the dataset into a phyloseq object for convenient data management
########################################
#protein counts as otu_table
exoP_counts<- exoP_filt %>% select(c("gene_caller_id", exoP_sample$Sample_name))
exoP_counts<- otu_table(data.frame(exoP_counts[, all_of(exoP_sample$Sample_name)], row.names = exoP_counts$gene_caller_id), taxa_are_rows = TRUE)

#gene ids as taxonomy table
annotation<- tax_table(as.matrix(data.frame(metaG_annotations, row.names = metaG_annotations$gene_caller_id)))

#merge into phyloseq object
exoP_obj0<- phyloseq(exoP_counts, annotation, sample_data(exoP_sample))

#remove proteins that were not observed 
exoP_obj0<- prune_taxa(taxa_sums(exoP_obj0)>0,exoP_obj0)

#save metaproteome phyloseq
saveRDS(exoP_obj0, "data/exoP_ps_raw.rds")

#print session info and clean the workspace
sessionInfo()
rm(list = ls())
gc()
