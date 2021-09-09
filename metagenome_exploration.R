########################################
#Plot KEGG pathways of each bin
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
require(ggplot2)
require(ggpubr)
require(DESeq2)
require(vegan)
require(pathview)


tol21rainbow<- c("#771155", "#AA4488","#CC99BB","#114477", 
                 "#4477AA","#117744","#117777","#88CCAA", 
                 "#77CCCC","#00ffff","#44AA77","#44AAAA", 
                 "#777711","#AAAA44","#DDDD77","#774411", 
                 "#AA7744","#DDAA77","#771122","#AA4455", 
                 "#DD7788")

###################
#Explore taxonomic composition of the metagenome
###################
#taxonomy table
tax_table <- read.csv(paste(wd,"metaG_analysis/metaG_anvio/05_ANVIO/spades-tax-names.txt",sep=""),
                      sep="\t", h= T)
#taxonomic classification
gene_tax <- read.csv(paste(wd,"metaG_analysis/metaG_anvio/05_ANVIO/spades-genes-taxonomy.txt",sep=""),
                           sep="\t", h= T)


#generate full taxonomy table
gene_tax_table <- merge(gene_tax,tax_table, by ="taxon_id", all.x = TRUE) %>% 
                    rename("t_order"="Class","t_class"="Order") #switch columns due to some taxonomy parsing bug in anvio


#summarize number of genes per taxon
gene_tax_table_summary <-gene_tax_table %>% group_by(Class) %>% 
                            summarize(Total = n())


#remove taxa that had less than 1000 genes or an unknown taxonomy
gene_tax_table_summary_fixed <- gene_tax_table_summary %>% 
                                  mutate(Class= case_when(Total < 1000 ~ "Other taxa",
                                                          grepl('Unknown', Class) ~ "Other taxa",
                                                          TRUE~ Class))
#fix the order of the classes for the figure
taxa_classes <- unique(gene_tax_table_summary_fixed$Class)
taxa_classes<- taxa_classes[!taxa_classes %in% c("Other taxa")]

gene_tax_table_summary_fixed$Class <- factor(gene_tax_table_summary_fixed$Class,
                                 levels=c(taxa_classes,"Other taxa"))


#plot taxonomy composition

#summarize once more for smoother visualization with black lines
gene_tax_table_summary_fixed <-gene_tax_table_summary_fixed %>% group_by(Class) %>% 
  summarize(Total.genes = sum(Total))

tax_comp.p <- ggplot(gene_tax_table_summary_fixed, aes(x= 1, y= Total.genes, fill = Class))+
  geom_bar(stat = "identity", color='black')+
  #scale_y_continuous(breaks = cumsum(gene_tax_table_summary_fixed$Total.genes) - gene_tax_table_summary_fixed$Total.genes/2,#produce coordinates for the annotation 
   #                  labels = gene_tax_table_summary_fixed$Class)+ # the labels
  coord_polar(theta ="y", direction = 1)+
  scale_fill_manual(values = tol21rainbow)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks=element_blank(),  # the axis ticks
        axis.title=element_blank(),  # the axis labels
        axis.text.y=element_blank(),
        text=element_text(size=14),legend.position = "bottom")
  

#save the plot
ggsave(paste0(wd,"/R_figures/metaG_tax_comp.pdf"), 
       plot = tax_comp.p,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)



###################
#Import metabolic estimates for each gene
###################
metaG_kofam_hits <- read.csv(paste(wd,"metaG_analysis/metaG_anvio/05_ANVIO/spades_kofam_hits.txt",sep=""), sep="\t", h= T)

metaG_KEGG_modules <- read.csv(paste(wd,"metaG_analysis/metaG_anvio/05_ANVIO/spades_modules.txt",sep=""), sep="\t", h= T)


metaG_KEGG_modules_by_gene <- metaG_KEGG_modules %>%
                                 tidyr::separate_rows(gene_caller_ids_in_module, sep = ',')
                                

#summarize KEGG modules
metaG_KEGG_modules_summary<- metaG_KEGG_modules_by_gene %>% 
                                filter(module_class =="Pathway modules") %>% 
                                  group_by(module_class, module_subcategory) %>% 
                                    summarize(Total.genes = n())



#removeplace modules with less than 15000 genes
metaG_KEGG_modules_summary_fixed <- metaG_KEGG_modules_summary %>% 
  mutate(module_subcategory= case_when(Total.genes < 1500 ~ "Other (<1500 genes)",
                          TRUE~ module_subcategory))

#fix the order of the classes for the figure
KEGG_modules <- unique(metaG_KEGG_modules_summary_fixed$module_subcategory)
KEGG_modules<- KEGG_modules[!KEGG_modules %in% c("Other (<1500 genes)")]

metaG_KEGG_modules_summary_fixed$module_subcategory <- factor(metaG_KEGG_modules_summary_fixed$module_subcategory,
                                             levels=c(KEGG_modules,"Other (<1500 genes)"))

#plot only pathway modules
metaG_KEGG_modules_summary_pathway<- metaG_KEGG_modules_summary_fixed %>% filter(module_class =="Pathway modules")

#summarize once more for smoother visualization with black lines
metaG_KEGG_modules_summary_pathway <-metaG_KEGG_modules_summary_pathway %>% group_by(module_subcategory) %>% 
  summarize(Total.genes = sum(Total.genes))

KEGG_modules_comp.p <- ggplot(metaG_KEGG_modules_summary_pathway, aes(x= 1, y= Total.genes, fill = module_subcategory))+
  geom_bar(stat = "identity", color='black')+
  #scale_y_continuous(breaks = cumsum(gene_tax_table_summary_fixed$Total.genes) - gene_tax_table_summary_fixed$Total.genes/2,#produce coordinates for the annotation 
  #                  labels = gene_tax_table_summary_fixed$Class)+ # the labels
  coord_polar(theta ="y", direction = 1)+
  scale_fill_manual(values = tol21rainbow)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks=element_blank(),  # the axis ticks
        axis.title=element_blank(),  # the axis labels
        axis.text=element_blank(),
        text=element_text(size=14),legend.position = "bottom")



#save the plot
ggsave(paste0(wd,"/R_figures/metaG_KEGG_comp.pdf"), 
       plot = KEGG_modules_comp.p,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)



################################################################################
#merge the taxonomy with KEGG modules and other annotations for metaP analyses
################################################################################
#gene annotation list by sources
gene_ids <- read.csv(paste(wd,"metaG_analysis/metaG_anvio/05_ANVIO/spades-gene-calls.txt",sep=""),
                     sep="\t", h= T)

#select only relevant columns and calculate the length of each protein
gene_annotations_df<- gene_ids %>% select("gene_callers_id","contig", "aa_sequence") %>% 
                        group_by(aa_sequence) %>% 
                         mutate(prot_length = nchar(aa_sequence))

for (src in c("COG20_FUNCTION","KeggGhostKoala","GO","Pfam","InterPro", "Hamap")){
  gene_annotations <- read.csv(paste(wd,"metaG_analysis/metaG_anvio/05_ANVIO/spades-",src,"-functions.txt",sep=""),
                               sep="\t", h= T)%>% select(gene_callers_id, accession, function.) %>% 
    group_by(gene_callers_id)%>%
    summarise_each(funs(paste(unique(.), collapse='|')),matches('^\\D+$')) %>% 
    plyr::rename(replace= c("accession"=paste(src, "accession", sep ="_"), "function."=paste(src, "function", sep ="_")))
  
  gene_annotations_df<- merge(gene_annotations_df,gene_annotations, by ="gene_callers_id", all = TRUE )
}


#merge taxonomy and gene ids together and calculate protein length
metaG_annotations <- gene_tax_table %>% 
  merge(gene_annotations_df, by ="gene_callers_id", all = TRUE) %>% 
  dplyr::rename("gene_caller_id"="gene_callers_id")

write.csv(metaG_annotations, "data/metaG_full_annotations.txt")  


#Explore un-annotated genes 
metaG_annotations_unknowns<-metaG_annotations %>%
  filter(is.na(COG20_FUNCTION_accession) & is.na(t_phylum))



######################################
#Investigate metabolic pathways in the bins
######################################
#Import metabolism estimates for the bins
modules_info <- read.csv(paste(wd,"metaG_analysis/metaG_anvio/07_METABOLISM/modules_info.txt",sep=""), sep="\t", h= T)

Bins_kofam_hits <- read.csv(paste(wd,"metaG_analysis/metaG_anvio/07_METABOLISM/Bins_kofam_hits.txt",sep=""), sep="\t", h= T)

Bins_modules <- read.csv(paste(wd,"metaG_analysis/metaG_anvio/07_METABOLISM/Bins_modules.txt",sep=""), sep="\t", h= T)


######################################
#plot KEGG modules number per bin
######################################
#subset complete modules only in categories of interest and summarize
Bins_modules_sum<- Bins_modules %>% 
  group_by(genome_name, module_category, module_subcategory) %>% 
  filter(module_completeness >0.7 & 
           module_category %in% c("Amino acid metabolism",
                                  "Carbohydrate metabolism",
                                  "Energy metabolism",
                                  "Lipid metabolism")) %>%
  summarize(Total = n()) %>% 
  mutate(genome_name = factor(genome_name, levels = c("Bin_84_1","Bin_76_1","Bin_5_2",
                                                      "Bin_2_1","Bin_5_3","Bin_38_1","Bin_2_2",
                                                      "Bin_102_1","Bin_134_1","Bin_134_2",
                                                      "Bin_179_1","Bin_115_2","Bin_115_1",
                                                      "Bin_12_1")))

#plot
Bins_modules_plot<- ggplot(Bins_modules_sum, aes(x=genome_name, y = Total, fill = module_subcategory))+
  geom_col()+
  facet_wrap(module_category~., scale = "free_y")+
  scale_fill_manual(values = tol21rainbow)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(size=14),legend.position = "bottom", 
        axis.title.x = element_blank(), axis.text.x = element_text(angle=90))


#save the plot
ggsave(paste0(wd,"/R_figures/KEGG_modules_per_bin.pdf"), 
       plot = Bins_modules_plot,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)


######################################
#Explore KEGG pathways per bin
######################################
Bins_kofam_hits_parsed <- Bins_kofam_hits %>% 
  separate_rows(modules_with_ko, sep = ',') %>% 
  dplyr::rename("kegg_module" = "modules_with_ko")


#merge the identified modules and the genes calls
Bins_gene_calls_KEGG_modules<-Bins_modules %>% 
  separate_rows(gene_caller_ids_in_module, sep = ',') %>% 
  dplyr::rename("gene_caller_id" = "gene_caller_ids_in_module") %>% 
  mutate(gene_caller_id=as.integer(gene_caller_id),
         genome_name = factor(genome_name, levels = c("Bin_84_1","Bin_76_1","Bin_5_2",
                                                         "Bin_2_1","Bin_5_3","Bin_38_1","Bin_2_2",
                                                         "Bin_179_1","Bin_115_2","Bin_115_1","Bin_102_1",
                                                         "Bin_12_1","Bin_150_1_1"))) %>% 
  left_join(Bins_kofam_hits_parsed[,c("genome_name","gene_caller_id","contig","kegg_module","ko","ko_definition")],
            by = c("genome_name","kegg_module","gene_caller_id"))


#export for KEGG website only the more complete modules
Bins_gene_calls_KEGG_modules %>%  
  filter(module_completeness >0.7) %>% 
  mutate(colour = "yellow") %>% 
  select("ko","colour") %>% 
  write.table("data/KEGG/Bins_ko.txt",
              row.names = FALSE,
              col.names = FALSE,quote = FALSE)


############################################################################
#plot the pathways of interest for each bin
############################################################################
#pathways of interest
pathways<- c("00010", #Glycolysis / Gluconeogenesis
             "00020", #TCA cycle
             "01200", #carbon metabolism
             "00190", #Oxidative phosphorylation
             "00071", #Fatty acid degradation
             "00250", #Alanine, aspartate and glutamate metabolism
             "00260", #Glycine, serine and threonine metabolism
             "00270", #Cysteine and methionine metabolism
             "00280", #Valine, leucine and isoleucine degradation
             "00310", #Lysine degradation
             "00330", #Arginine and proline metabolism
             "00340", #Histidine metabolism
             "00350", #Tyrosine metabolism
             "00360", #Phenylalanine metabolism
             "00380", #Tryptophan metabolism
             "00400" #Phenylalanine, tyrosine and tryptophan biosynthesis
)



#plot in a loop all the pathways of interest for each bin
for (bin in unique(Bins_gene_calls_KEGG_modules$genome_name)) {
Bins_gene_calls_KEGG_modules_complete<- Bins_gene_calls_KEGG_modules %>%  
    filter(genome_name == bin & module_completeness >0.7) 

  pathview(gene.data = Bins_gene_calls_KEGG_modules_complete$ko, 
           pathway.id =pathways,
           species = "ko", 
           keys.align = "y", 
           kegg.native = T, both.dirs = TRUE,
           discrete	=list(gene=TRUE),
           res = 300, cex = 0.25,
           out.suffix = bin)
}

############################################################################
#explore each bin manualy
############################################################################

#Bin 84_1 - Pseudoalteromonas phenolica - length 3.92Mbp (C93/R0)

Bin_84_1_KEGG_modules <- Bins_gene_calls_KEGG_modules %>% 
  filter(genome_name == "Bin_84_1" & module_completeness >0.70) %>% 
  select("genome_name","kegg_module", "module_category", "module_subcategory", 
         "module_completeness","contig","gene_caller_id","ko","ko_definition")

#Bin_115_1 - Bermanella sp002683575 (Pseudomonadales) - length 4.86Mbp (C98.6/R0)
Bin_115_1_KEGG_modules <- Bins_gene_calls_KEGG_modules %>% 
  filter(genome_name == "Bin_115_1" & module_completeness >0.70) %>% 
  select("genome_name","kegg_module", "module_category", "module_subcategory", 
         "module_completeness","contig","gene_caller_id","ko","ko_definition")

#Bin_115_2 - Glaciecola sp000155775 (Enterobacterales) - length 2.22Mbp (C100/R0)
Bin_115_2_KEGG_modules <- Bins_gene_calls_KEGG_modules %>% 
  filter(genome_name == "Bin_115_2" & module_completeness >0.75) %>% 
  select("genome_name","kegg_module", "module_category", "module_subcategory", 
         "module_completeness","contig","gene_caller_id","ko","ko_definition")

#Bin_76_1 - family Nitrincolaceae - length 3.78Mbp (C100/R0)
Bin_76_1_KEGG_modules <- Bins_gene_calls_KEGG_modules %>% 
  filter(genome_name == "Bin_76_1" & module_completeness >0.75) %>% 
  select("genome_name","kegg_module", "module_category", "module_subcategory", 
         "module_completeness","contig","gene_caller_id","ko","ko_definition")

#Bin_179_1 - family Cellvibrionaceae - length 3.18Mbp (C93/R2.8)
Bin_179_1_KEGG_modules <- Bins_gene_calls_KEGG_modules %>% 
  filter(genome_name == "Bin_179_1" & module_completeness >0.75) %>% 
  select("genome_name","kegg_module", "module_category", "module_subcategory", 
         "module_completeness","contig","gene_caller_id","ko","ko_definition")

#Bin_38_1 - Vibrio - length 1.7Mbp (C83/R0)
Bin_38_1_KEGG_modules <- Bins_gene_calls_KEGG_modules %>% 
  filter(genome_name == "Bin_38_1" & module_completeness >0.75) %>% 
  select("genome_name","kegg_module", "module_category", "module_subcategory", 
         "module_completeness","contig","gene_caller_id","ko","ko_definition")

#Bin_102_1 - Kordiimonas lacus - length 2.18Mbp (C73/R0)
Bin_102_1_KEGG_modules <- Bins_gene_calls_KEGG_modules %>% 
  filter(genome_name == "Bin_102_1" & module_completeness >0.75) %>% 
  select("genome_name","kegg_module", "module_category", "module_subcategory", 
         "module_completeness","contig","gene_caller_id","ko","ko_definition")

#Bin_12_1 - Saccharospirillum sp003054965 - length 3.75Mbp (C95.8/R1.4)
Bin_12_1_KEGG_modules <- Bins_gene_calls_KEGG_modules %>% 
  filter(genome_name == "Bin_12_1" & module_completeness >0.75) %>% 
  select("genome_name","kegg_module", "module_category", "module_subcategory", 
         "module_completeness","contig","gene_caller_id","ko","ko_definition")

#Bin_5_2 - family Nitrincolaceae - length 5.19Mbp (C100/R0)
Bin_5_2_KEGG_modules <- Bins_gene_calls_KEGG_modules %>% 
  filter(genome_name == "Bin_5_2" & module_completeness >0.75) %>% 
  select("genome_name","kegg_module", "module_category", "module_subcategory", 
         "module_completeness","contig","gene_caller_id","ko","ko_definition")

#Bin_5_3 - Reinekea blandensis (Pseudomonadales) - length 4.08Mbp (C100/R2.8)
Bin_5_3_KEGG_modules <- Bins_gene_calls_KEGG_modules %>% 
  filter(genome_name == "Bin_5_3" & module_completeness >0.75) %>% 
  select("genome_name","kegg_module", "module_category", "module_subcategory", 
         "module_completeness","contig","gene_caller_id","ko","ko_definition")

#Bin_2_1 - Alteromonas - length 3.7Mbp (C78.9/R0)
Bin_2_1_KEGG_modules <- Bins_gene_calls_KEGG_modules %>% 
  filter(genome_name == "Bin_2_1" & module_completeness >0.75) %>% 
  select("genome_name","kegg_module", "module_category", "module_subcategory", 
         "module_completeness","contig","gene_caller_id","ko","ko_definition")

#Bin_2_2 - Alteromonas - length 5.19Mbp (C80.3/R0)
Bin_2_2_KEGG_modules <- Bins_gene_calls_KEGG_modules %>% 
  filter(genome_name == "Bin_2_2" & module_completeness >0.75) %>% 
  select("genome_name","kegg_module", "module_category", "module_subcategory", 
         "module_completeness","contig","gene_caller_id","ko","ko_definition")

#Bin_150_1_1 - Rhodobacteraceae - length 2.6Mbp (C97.2/R1.4)
Bin_150_1_1_KEGG_modules <- Bins_gene_calls_KEGG_modules %>% 
  filter(genome_name == "Bin_150_1_1" & module_completeness >0.75) %>% 
  select("genome_name","kegg_module", "module_category", "module_subcategory", 
         "module_completeness","contig","gene_caller_id","ko","ko_definition")


