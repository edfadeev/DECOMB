# https://uclouvain-cbio.github.io/BSS2019/figs/cancer_9x9.html#data-exploration

require(dplyr)
require(tidyr)
require(MSqRob)
require(MSnbase)
require(tidyverse)
require(limma)
require(stringr)

#set working directory
wd <- "/Users/eduardfadeev/Google Drive (dr.eduard.fadeev@gmail.com)/DECOMB/"

#import taxonomy and annotation of each gene in the reference metagenome
#taxonomy table
tax_table <- read.csv(paste(wd,"metaG_analysis/metaG_anvio/05_ANVIO/spades_merged/spades-tax-names.txt",sep=""),
                     sep="\t", h= T)
#taxonomic classification
gene_tax_table <- read.csv(paste(wd,"metaG_analysis/metaG_anvio/05_ANVIO/spades_merged/spades-genes-taxonomy.txt",sep=""),
                           sep="\t", h= T)
#annotation
gene_fun_table <- read.csv(paste(wd,"metaG_analysis/metaG_anvio/05_ANVIO/spades_merged/spades-genes-fun-merged.txt",sep=""),
                          sep=" ", h= T)
#merge all information together
genes_meta <- merge(gene_tax_table,tax_table, by ="taxon_id", all.x = TRUE) %>% 
                merge(gene_fun_table, by ="gene_callers_id", all = TRUE) %>% 
                  select(c("gene_callers_id", "contig", "start", "stop", "partial", "accession","function.","e_value", "taxon_id", "t_genus", "t_species","aa_sequence"))
                

#import sample list of metaproteomes
prot_sample <- read.csv(paste(wd,"metaP_analysis/quant-no-groups/DECOMB-consensus-no-groups_InputFiles.txt", sep=""),sep="\t", h= T) %>% 
  mutate(Replicate = case_when(grepl("20200622_TinkaraTinta2",File.Name) == TRUE ~ 2, TRUE ~ 1),
          File.Name =   tools::file_path_sans_ext(gsub("Z:\\\\EFadeev\\\\DECOMB_raw_data\\\\20200622_TinkaraTinta2\\\\|Z:\\\\EFadeev\\\\DECOMB_raw_data\\\\20200615_TinkaraTinta\\\\",
                                                      "", File.Name)),
          Fraction = case_when(grepl("MP",File.Name) == TRUE ~ "metaP", TRUE ~ "ExoP"),
           Type = case_when(grepl("C",File.Name) ==TRUE ~"Control", 
                            grepl("J",File.Name) ==TRUE ~"Jelly",
                            grepl("T0",File.Name) ==TRUE ~ "T0"),
         Sample.ID = paste(File.ID, File.Name, Type, Replicate, sep="_")) # %>% 
  # separate(File.Name, into = c("Sample","Fraction"), sep ="_", remove = FALSE)

#import proteins data
protein_raw <- read.csv(paste(wd,"metaP_analysis/quant-no-groups/DECOMB-consensus-no-groups_Proteins.txt", sep=""),sep="\t", h= T)

#rename sample columns and filter out low confidence proteins
protein_filt <- protein_raw %>% 
  dplyr::rename(gene_callers_id = Accession) %>% 
  select_at(vars(!contains("Found"))) %>% 
  rename_with(~gsub("Abundance\\.|\\.Sample","",.), everything()) %>% 
  rename_at(vars(prot_sample$File.ID), ~ prot_sample$Sample.ID) %>% 
  left_join(genes_meta, by ="gene_callers_id", all = TRUE) %>% 
  filter(Number.of.PSMs >=2 , Number.of.Unique.Peptides>=1)

#sum the two runs together

# NAAF transformation
protein_trans <- protein_filt %>% dplyr::mutate_at(c(prot_sample$Sample.ID), funs(replace_na(., 0))) %>% 
  dplyr::mutate_at(c(prot_sample$Sample.ID), funs(./Number.of.AAs)) %>% 
  dplyr::mutate_at(c(prot_sample$Sample.ID), funs(NAAF = ./ sum(.)))

#NMDS
prot.pca <- vegan::metaMDS(t(protein_trans[74:115]))


protein_nmds <- as.data.frame(prot.pca$points) %>% tibble::rownames_to_column() %>% 
                   separate(rowname, into = c("ID","Sample","Fraction","Type","Replicate", "Method"), sep ="_", remove = TRUE) %>% 
                    mutate(Type = case_when(Sample %in% c("C1","C2","C3") ~"Control", 
                                            Sample %in% c("J1","J2","J3") ~"Jelly",
                                                      Sample == "T0" ~ "T0"))

ggplot(protein_nmds, aes(x=MDS1,y=MDS2, colour = Type, shape = as.factor(Type), label = Replicate))+
  geom_point(size = 5)+
  geom_text(size = 5, nudge_y = -0.2)



# Quantile normalisation : the aim is to give different distributions the
# same statistical properties
quantile_normalisation <- function(df){
  
  # Find rank of values in each column
  df_rank <- map_df(df,rank,ties.method="average")
  # Sort observations in each column from lowest to highest 
  df_sorted <- map_df(df,sort)
  # Find row mean on sorted columns
  df_mean <- rowMeans(df_sorted)
  
  # Function for substiting mean values according to rank 
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  # Replace value in each column with mean according to rank 
  df_final <- map_df(df_rank,index_to_mean, my_mean=df_mean)
  
  return(df_final)
}


dat_norm <- protein_no_na[20:61] %>% quantile_normalisation()
  


protein_no_na <- protein_filt %>% mutate_if(is.numeric, funs(replace_na(., 0))) #%>% mutate_if(is.numeric, funs(log(., 2)))


prot_mds <- metaMDS(dat_norm)


















peptide_raw <- read.csv(paste(wd,"metaP_analysis/quant-no-groups/DECOMB-consensus-no-groups_PeptideGroups.txt", sep=""),sep="\t", h= T)



peptide_no_mod <- peptide_raw %>%filter(Modifications == "")

pepData<-readMSnSet2(peptide_no_mod,ecol=17:58,fnames="Annotated.Sequence",sep="\t")


sampleNames(pepData) <- pepData %>% sampleNames %>% str_replace(., pattern="Abundance\\.", replacement="") %>% 
                            str_replace(., pattern="\\.Sample", replacement="")    

pepData<-selectFeatureData(pepData,fcol=c("Checked","Confidence","Annotated.Sequence",
                                          "Modifications.in.Master.Proteins","Contaminant",
                                          "Qvality.PEP","Qvality.q.value","Number.of.Protein.Groups",
                                          "Number.of.Proteins","Number.of.PSMs","Master.Protein.Accessions"))



plotNA(pepData)

plotDensities(exprs(pepData))
nrow(pepData)

pepData_log <- log(pepData, base = 2)
plotDensities(exprs(pepData_log))

pepData_quant <- normalise(pepData_log, "quantiles")

plotDensities(exprs(pepData_quant))

pepData<-selectFeatureData(pepData)



plotMDS(exprs(pepData_quant),col=as.double(pData(pepData)$condition))



                
                
plotDensities(protein_filt[,c(19:61)])





test <- metaMDS(protein_filt[,c(19:61)],distance = "bray", k = 2, plot = TRUE)


plot(x, display = c("sites", "species"), choices = c(1, 2),
     type = "p", shrink = FALSE,  ...)




