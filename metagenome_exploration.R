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
                                  mutate(Class= case_when(Total < 1000 ~ "Other taxa (<1000 genes)",
                                                          grepl('Unknown', Class) ~ "Other taxa (<1000 genes)",
                                                          TRUE~ Class))
#fix the order of the classes for the figure
taxa_classes <- unique(gene_tax_table_summary_fixed$Class)
taxa_classes<- taxa_classes[!taxa_classes %in% c("Other taxa (<1000 genes)")]

gene_tax_table_summary_fixed$Class <- factor(gene_tax_table_summary_fixed$Class,
                                 levels=c(taxa_classes,"Other taxa (<1000 genes)"))


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
metaG_kofam_hits <- read.csv(paste(wd,"metaG_analysis/metaG_anvio/05_ANVIO/spades-Kofam_hits.txt",sep=""), sep="\t", h= T)

metaG_KEGG_modules <- read.csv(paste(wd,"metaG_analysis/metaG_anvio/05_ANVIO/spades-Kofam_modules.txt",sep=""), sep="\t", h= T)


metaG_KEGG_modules_by_gene <- metaG_KEGG_modules %>%
                                 tidyr::separate_rows(gene_caller_ids_in_module, sep = ',')
                                

#summarize KEGG modules
metaG_KEGG_modules_summary<- metaG_KEGG_modules_by_gene %>% 
                                filter(module_class =="Pathway modules") %>% 
                                  group_by(module_class, module_subcategory) %>% 
                                    summarize(Total.genes = n())



#remove modules with less than 1000 genes
metaG_KEGG_modules_summary_fixed <- metaG_KEGG_modules_summary %>% 
  mutate(module_subcategory= case_when(Total.genes < 1000 ~ "Other (<1000 genes)",
                          TRUE~ module_subcategory))

#fix the order of the classes for the figure
KEGG_modules <- unique(metaG_KEGG_modules_summary_fixed$module_subcategory)
KEGG_modules<- KEGG_modules[!KEGG_modules %in% c("Other (<1000 genes)")]

metaG_KEGG_modules_summary_fixed$module_subcategory <- factor(metaG_KEGG_modules_summary_fixed$module_subcategory,
                                             levels=c(KEGG_modules,"Other (<1000 genes)"))

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


#save combined figure
ggarrange(tax_comp.p, KEGG_modules_comp.p, #heights = c(2,1.2),
          ncol = 2, nrow = 1, align = "v", legend = "bottom",
          legend.grob = do.call(rbind, c(list(get_legend(tax_comp.p),get_legend(KEGG_modules_comp.p)), size="first")))

#save the plot
ggsave(paste0(wd,"/R_figures/metaG_taxa_KEGG_combined.pdf"), 
       plot = last_plot(),
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



