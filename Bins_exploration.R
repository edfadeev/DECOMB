#load libraries
require(dplyr)
require(tidyr)
require(ggplot2)
require(pathview)

source("scripts/extra_functions.R")

######################################
#Investigate metabolic pathways in the bins
######################################
Bins_modules <- read.csv("data/metagenome/Bins_modules.txt", sep="\t", h= T)


Bins_total<- Bins_modules %>% 
  group_by(genome_name) %>% 
  filter(module_class == "Pathway modules" &
           module_is_complete == "True" &
           #module_completeness >0.9 & 
           module_category %in% c("Amino acid metabolism",
                                  "Carbohydrate metabolism",
                                  "Energy metabolism",
                                  "Lipid metabolism",
                                  "Metabolism of cofactors and vitamins")) %>%
  summarize(Total = n())

#list the bins with more modules than the mean (26)
top_bins<- Bins_total %>% filter(Total > 26) %>% select(genome_name) %>% tibble::deframe()

#subset complete modules only in categories of interest and summarize
Bins_modules_sum<- Bins_modules %>% 
  group_by(genome_name, module_category, module_subcategory) %>% 
  filter(genome_name %in% top_bins &
           module_class == "Pathway modules" &
           module_is_complete == "True" &
           #module_completeness >0.9 & 
           module_category %in% c("Amino acid metabolism",
                                  "Carbohydrate metabolism",
                                  "Energy metabolism",
                                  "Lipid metabolism",
                                  "Metabolism of cofactors and vitamins")) %>%
  summarize(Total = n())

#plot
Bins_modules.p <- ggplot(Bins_modules_sum, aes(x=genome_name, y = Total, fill = module_subcategory))+
  geom_col()+
  #facet_wrap(~genome_name, scale = "free_y")+
  scale_fill_manual(values = tol21rainbow)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(size=20),legend.position = "right", 
        axis.title.x = element_blank(), axis.text.x = element_text(angle=90))


#save the plot
ggsave("Figures/KEGG_modules_per_bin.png", 
       plot = Bins_modules.p,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)

######################################
#Explore KEGG pathways per bin
######################################
Bins_kofam_hits_parsed <- read.csv("data/metagenome/Bins_kofam_hits.txt",
                                   sep="\t", h= T,quote="")%>% 
                          separate_rows(modules_with_ko, sep = ',') %>% 
                          dplyr::rename("kegg_module" = "modules_with_ko")


#merge the identified modules and the genes calls of only the complete modules
Bins_gene_calls_KEGG_modules<-Bins_modules %>% 
  separate_rows(gene_caller_ids_in_module, sep = ',') %>% 
  dplyr::rename("gene_caller_id" = "gene_caller_ids_in_module") %>% 
  filter(module_is_complete == "True") %>% 
  mutate(gene_caller_id=as.integer(gene_caller_id)) %>% 
  left_join(Bins_kofam_hits_parsed[,c("genome_name","gene_caller_id","contig","kegg_module","ko","ko_definition")],
            by = c("genome_name","kegg_module","gene_caller_id"))

#export for KEGG website 
#Bins_gene_calls_KEGG_modules %>%  
#  mutate(colour = "yellow") %>% 
#  select("ko","colour") %>% 
#  unique() %>% 
#  write.table("data/KEGG/Bin_ko.txt",
#              row.names = FALSE,
#              col.names = FALSE,quote = FALSE)

############################################################################
#explore each bin manualy
############################################################################

#Bin 84 - Pseudoalteromonas phenolica - length 4.27Mbp (C85/R0)
Bin_84_KEGG_modules <- Bins_gene_calls_KEGG_modules %>% 
  filter(genome_name == "Bin_84" & module_completeness >0.75) %>% 
  select("genome_name","kegg_module", "module_category", "module_subcategory", 
         "module_completeness","contig","gene_caller_id","ko","ko_definition")

#Bin_76 - Marinobacterium jannaschii - length 5.32Mbp (C86/R0)
Bin_76_KEGG_modules <- Bins_gene_calls_KEGG_modules %>% 
  filter(genome_name == "Bin_76" & module_completeness >0.75) %>% 
  select("genome_name","kegg_module", "module_category", "module_subcategory", 
         "module_completeness","contig","gene_caller_id","ko","ko_definition")


#Bin_102 - Kordiimonas lacus - length 3.24Mbp (C100/R2)
Bin_102_KEGG_modules <- Bins_gene_calls_KEGG_modules %>% 
  filter(genome_name == "Bin_102" & module_completeness >0.75) %>% 
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
    filter(genome_name == bin & module_is_complete == "True") 
  
  pathview(gene.data = Bins_gene_calls_KEGG_modules_complete$ko, 
           pathway.id =pathways,
           species = "ko", 
           keys.align = "y", 
           kegg.native = T, both.dirs = TRUE,
           discrete	=list(gene=TRUE),
           res = 300, cex = 0.25,
           out.suffix = bin,
           kegg.dir="./Figures/KEGG/")
}



#to directly view the pathway map
see_pathview <- function(..., save_image = FALSE)
{
  msg <- capture.output(pathview::pathview(...), type = "message")
  msg <- grep("image file", msg, value = T)
  filename <- sapply(strsplit(msg, " "), function(x) x[length(x)])
  img <- png::readPNG(filename)
  grid::grid.raster(img)
  if(!save_image) invisible(file.remove(filename))
}

see_pathview(gene.data = gse16873.d[, 1], 
             pathway.id = demo.paths$sel.paths[1],
             species = "hsa", 
             out.suffix = "gse16873", 
             kegg.native = T)
