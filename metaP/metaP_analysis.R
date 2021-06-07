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
  mutate(File.Name =   tools::file_path_sans_ext(gsub("Z:\\\\EFadeev\\\\DECOMB_raw_data\\\\20200622_TinkaraTinta2\\\\|Z:\\\\EFadeev\\\\DECOMB_raw_data\\\\20200615_TinkaraTinta\\\\",
                                                      "", File.Name)),
         Sample.ID = paste(File.ID,File.Name, sep="_")) %>% 
  separate(File.Name, into = c("Sample","Fraction"), sep ="_", remove = FALSE)

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




