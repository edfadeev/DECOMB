#load libraries
require(dplyr)
require(tidyr)
require(ggplot2)


require(phyloseq)

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
#generate full taxonomy table
gene_tax_table <- merge(read.csv("data/metagenome/spades-tax-names.txt",sep="\t", h= T),#taxonomy table
                        read.csv("data/metagenome/spades-genes-taxonomy.txt",sep="\t", h= T),#taxonomic classification
                        by ="taxon_id", all.x = TRUE) %>% 
                  dplyr::rename(Class = t_order ,Order = t_class) #switch columns due to some taxonomy parsing bug in anvio

#summarize number of genes per taxon and remove taxa that had less than 1000 genes or an unknown class
gene_tax_table_summary <-gene_tax_table %>% 
                            group_by(Class) %>% 
                            summarize(Total = n()) %>% 
                            mutate(Class= case_when(Total < 1000 ~ "Other taxa (<1000 genes)",
                                                    grepl('Unknown', Class) ~ "Unknown class",
                                                    TRUE~ Class))
#fix the order of the classes for the figure
taxa_classes <- unique(gene_tax_table_summary$Class)
taxa_classes<- taxa_classes[!taxa_classes %in% c("Other taxa (<1000 genes)")]
gene_tax_table_summary$Class <- factor(gene_tax_table_summary$Class,
                                 levels=c(taxa_classes,"Other taxa (<1000 genes)"))

#summarize once more for smoother visualization and plot
gene_tax_table_summary %>% 
  group_by(Class) %>% 
  summarize(Total.genes = sum(Total)) %>%
  ggplot(aes(x= 1, y= Total.genes, fill = Class))+
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
ggsave("/Figures/metaG_tax_comp.pdf",
       plot = last_plot(),
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)



###################
#Import metabolic estimates for each gene
###################

metaG_KEGG_modules_by_gene <- read.csv("data/metagenome/spades-Kofam_modules.txt", sep="\t", h= T)%>%
  tidyr::separate_rows(gene_caller_ids_in_module, sep = ',')

#summarize KEGG modules and remove modules with less than 1000 genes
metaG_KEGG_modules_summary<- metaG_KEGG_modules_by_gene %>% 
  filter(module_class =="Pathway modules") %>% 
  group_by(module_class, module_subcategory) %>% 
  summarize(Total.genes = n()) %>% 
  mutate(module_subcategory= case_when(Total.genes < 1000 ~ "Other modules (<1000 genes)",
                                       TRUE~ module_subcategory))

#fix the order of the classes for the figure
KEGG_modules <- unique(metaG_KEGG_modules_summary$module_subcategory)
KEGG_modules<- KEGG_modules[!KEGG_modules %in% c("Other modules (<1000 genes)")]
metaG_KEGG_modules_summary$module_subcategory <- factor(metaG_KEGG_modules_summary$module_subcategory,
                                       levels=c(KEGG_modules,"Other modules (<1000 genes)"))

#summarize once more for smoother visualization with black lines
metaG_KEGG_modules_summary %>% 
  group_by(module_subcategory) %>% 
  summarize(Total.genes = sum(Total.genes)) %>% 
  ggplot(aes(x= 1, y= Total.genes, fill = module_subcategory))+
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
ggsave("Figures/metaG_taxa_KEGG_combined.pdf", 
       plot = last_plot(),
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)







################################################################################
#merge the taxonomy with KEGG modules and other annotations for metaP analyses
################################################################################

metaG_kofam_hits <- read.csv("data/metagenome/spades-Kofam_hits.txt", sep="\t", h= T)




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



