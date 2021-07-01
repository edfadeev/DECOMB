# https://uclouvain-cbio.github.io/BSS2019/figs/cancer_9x9.html#data-exploration

require(dplyr)
require(tidyr)
require(MSqRob)
require(MSnbase)
require(tidyverse)
require(limma)
require(stringr)

#set working directory
#macOS
wd <- "/Users/eduardfadeev/Google Drive (dr.eduard.fadeev@gmail.com)/DECOMB/"

#Windows
wd <- "D:/Postdoc-Vienna/DECOMB/"

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


sample_nmds <- as.data.frame(prot.pca$points) %>% tibble::rownames_to_column() %>% 
                   separate(rowname, into = c("ID","Sample","Fraction","Type","Replicate", "Method"), sep ="_", remove = TRUE) %>% 
                    mutate(Type = case_when(Sample %in% c("C1","C2","C3") ~"Control", 
                                            Sample %in% c("J1","J2","J3") ~"Jelly",
                                                      Sample == "T0" ~ "T0"))

protein_nmds <- as.data.frame(prot.pca$species) %>% tibble::rownames_to_column() %>% 
  separate(rowname, into = c("ID","Sample","Fraction","Type","Replicate", "Method"), sep ="_", remove = TRUE) %>% 
  mutate(Type = case_when(Sample %in% c("C1","C2","C3") ~"Control", 
                          Sample %in% c("J1","J2","J3") ~"Jelly",
                          Sample == "T0" ~ "T0"))

ggplot()+
  geom_point(data=protein_nmds, aes(x=MDS1,y=MDS2), size = 1, alpha = 0.3, colour = "gray50")+
  geom_point(data=sample_nmds, aes(x=MDS1,y=MDS2, colour = Type, shape = as.factor(Type), label = Replicate),
             size = 5)+
  geom_text(data=sample_nmds, aes(x=MDS1,y=MDS2, colour = Type, shape = as.factor(Type), label = Replicate),
            size = 5, nudge_y = -0.2)+
  theme_bw()


#enrichment test
require(DESeq2)


protein_no_na <- protein_filt %>% mutate_if(is.numeric, funs(replace_na(., 0))) %>% 
  mutate_if(is.numeric,as.integer)

dds <- DESeqDataSetFromMatrix(countData=protein_no_na[,c(4,20:61)], 
                              colData=prot_sample, 
                              design=~Type, tidy = TRUE)

prot.DEseq <- DESeq(dds)
prot.DEseq.res <- results(prot.DEseq)



require("msmsTests")

samples_frac <- prot_sample %>% filter(Type != "T0", Fraction == "metaP") %>% select(c(9:11))
row.names(samples_frac)<- samples_frac$Sample.ID

protein_frac <- protein_no_na %>% select(c("gene_callers_id", samples_frac$Sample.ID))

prot_MSnSet <- readMSnSet2(protein_frac, ecol = samples_frac$Sample.ID, fnames = 1)

phenoData(prot_MSnSet) <- AnnotatedDataFrame(samples_frac)
                    

e <- pp.msms.data(prot_MSnSet)

null.f <- "y~1"
alt.f <- "y~Type"
div <- apply(exprs(e),2,sum)
edgeR_res <- msms.edgeR(e, form0= null.f, form1= alt.f, div=div, fnm="Type")


test<- test.results(edgeR_res, e, gpf = pData(e)$Type, gp1="Control",gp2="Jelly",
                    alpha = 0.05, minLFC=1, div= div,
                    method="BH")

res.volcanoplot(test$tres,max.pval=0.1, maxy=2, maxx = 30)



results<- data.frame(test$tres, gene_callers_id = as.integer(row.names(test$tres))) %>% 
            left_join(genes_meta, by = "gene_callers_id") %>% 
  filter(DEP =="TRUE")

