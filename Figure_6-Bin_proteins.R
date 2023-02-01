require(pathview)
require(phyloseq)
require(DESeq2)
require(dplyr)
require(ggplot2)

source("scripts/extra_functions.R")

#########################################################
#Explore all proteins in binned genes                 ###
#########################################################
#load metaproteome phyloseq object
#conduct NSAF transformation
#subset to the genes which were observed also in proteins
metaP_in_bins<- readRDS("data/metaproteome/metaP_runB_merged.rds") %>% 
                add_nsaf(., "prot_length") %>% 
                subset_taxa(., !is.na(Bin)) %>% 
                prune_taxa(taxa_sums(.)>0, .)

#melt phyloseq into a dataframe for ploting
metaP_in_bins.long <- psmelt(metaP_in_bins) %>% 
  mutate(Type = factor(Type, levels =c("Cellular","Exocellular")))

#sum up each bin
prot_nsaf.Bin.agg <- metaP_in_bins.long %>% 
  group_by(Sample_name, Type, Treatment, Bin) %>% 
  filter(Abundance> 0) %>% 
  summarize(Tot.abundance = sum(Abundance),
            Total_p = length(Abundance)) %>% 
  group_by(Type, Treatment, Bin) %>% 
  summarize(Tot.abundance = mean(Tot.abundance))

#remove below 1% ra
Bins <- sort(as.character(unique(prot_nsaf.Bin.agg$Bin[!prot_nsaf.Bin.agg$Tot.abundance<0.01])))

prot_nsaf.Bin.agg$Bin[prot_nsaf.Bin.agg$Tot.abundance<0.01] <- "Other bins"
prot_nsaf.Bin.agg$Bin[is.na(prot_nsaf.Bin.agg$Bin)] <- "Other bins"

prot_nsaf.Bin.agg$Bin <- factor(prot_nsaf.Bin.agg$Bin,
                                levels=c(Bins,"Other bins"))

#plot
prot_Bin.p<- ggplot(prot_nsaf.Bin.agg, 
                         aes(x = Treatment, y = Tot.abundance, 
                             fill = Bin)) + 
  facet_wrap(.~Type, scales = "free_x") +
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values = tol21rainbow)+ 
  #guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Protein proportions (>1%) \n")+
  theme_EF+
  theme(legend.position = "bottom")


#save the plot
ggsave("./Figures/Figure_7-Proteins_by_bin.pdf", 
       plot = prot_Bin.p,
       units = "mm",
       width = 90, height = 90, 
       scale = 3,
       dpi = 300)

#abundance
Prot_abund_per_bin <- prot_nsaf.Bin.agg %>% 
  mutate(Type = factor(Type, levels =c("Cellular","Exocellular"))) %>% 
  group_by(Type, Treatment, Sample_name, Bin) %>% 
  summarize(Total_prot = sum(Tot.abundance)) %>% 
  group_by(Type, Treatment, Bin) %>% 
  summarize(Min_abund = min(Total_prot),
            Max_abund = max(Total_prot)) 

#how many proteins there are from each bin?
#sum all proteins per sample
metaP_totals<- readRDS("data/metaproteome/metaP_runB_merged.rds") %>% 
  add_nsaf(., "prot_length") %>% 
  psmelt() %>% 
  mutate(Type = factor(Type, levels =c("Cellular","Exocellular"))) %>% 
  group_by(Sample_name, Type, Treatment) %>% 
  filter(Abundance> 0) %>% 
  summarize(Proteins = length(Abundance))

#calculate proportion of all proteins per bin per sample
total_prot_per_bin <- prot_nsaf.Bin.agg %>% 
  mutate(Type = factor(Type, levels =c("Cellular","Exocellular"))) %>% 
  left_join(metaP_totals, by = c("Type", "Treatment", "Sample_name")) %>% 
  group_by(Type, Treatment, Sample_name, Bin) %>% 
  summarize(Prop = Total_p/Proteins,
            Total_p = Total_p,
            Proteins = Proteins)

#########################################################
#Explore sig. enriched proteins in binned genes       ###
#########################################################
#import prot. enrichment results 
DESeq_res <- read.csv("data/DESEq_res.csv", row.names = 1) %>% 
  mutate(Fraction = factor(Fraction, levels =c("Cellular","Exocellular")))

#calculate how many enriched proteins there are per bin
DESeq_res_bins_totals <- DESeq_res %>% 
  filter(!is.na(Bin)) %>% 
  group_by(Bin, Type) %>%  summarize (n())

#explore manualy Bin_84 (Pseudoalteromonas)
DESeq_res_bin_84_COG20_fun <- DESeq_res %>%  
  filter(Bin =="Bin_84", Type %in% c("Cteno-OM","Control")) %>% 
  select(Type, Fraction, COG20_CATEGORY_accession, COG20_CATEGORY_function, COG20_FUNCTION_function, log2FoldChange)

DESeq_res_bin_84_KEGG <- DESeq_res %>%  
  filter(Bin =="Bin_84", Type %in% c("Cteno-OM","Control")) %>% 
  select(Type, Fraction, COG20_FUNCTION_function, contains("KEGG"),log2FoldChange)


#########################################################
#Explore KEGG pathways of Bin_84                      ###
#########################################################
#import annotations by gene
Bin_Kofam_hits <- read.csv("./data/metagenome/Bins_kofam_hits.txt", 
                           sep="\t", h= T, quote = "")%>% 
  dplyr::rename("gene_callers_id" = "gene_caller_id",
                "Bin" = "db_name") %>% 
  filter(Bin =="Bin_84")

Bin_84_Kofam_merged <- Bin_Kofam_hits %>% 
  left_join(DESeq_res[,c("gene_callers_id","Bin","log2FoldChange")],
            by = c("gene_callers_id","Bin")) %>% 
  filter(log2FoldChange>0 | is.na(log2FoldChange))

Bins_gene_calls_KEGG_modules<-Bin_84_Kofam_merged %>% 
  filter(!is.na(log2FoldChange)) %>% 
  tidyr::separate_rows(modules_with_ko, sep = ',')%>% 
  dplyr::rename("kegg_module" = "modules_with_ko") %>% 
  select(kegg_module, ko) %>% 
  unique() %>% 
  group_by(kegg_module) %>% 
  summarize(Protein_KO = paste0(ko, collapse = ",")) 

#Investigate metabolic pathways in the bin and produce a summary table
Bin_84_KEGG_summary_table <- read.csv("data/metagenome/Bins_modules.txt", sep="\t", h= T) %>% 
  filter(#module_completeness >0.7 &
         genome_name == "Bin_84") %>% 
         select(module_category, module_subcategory, module_name,
                kegg_module, module_completeness, kofam_hits_in_module)  %>% 
  left_join(Bins_gene_calls_KEGG_modules, by = "kegg_module") %>% 
  mutate(module_completeness= paste0(round(module_completeness,2)*100,"%")) %>% 
  filter(!is.na(Protein_KO))

#export the summary table
write.table(Bin_84_KEGG_summary_table, "data/Table_2-KEGG_modules.txt", sep="\t", row.names = FALSE, quote = FALSE)


#produce list of KOs for mapping
DESeq_res_bin_84_KO_list <- Bin_84_Kofam_merged %>%  
  mutate(log2FoldChange = case_when(is.na(log2FoldChange)==TRUE ~ -1, 
                                    log2FoldChange>0 ~ 1)) %>% 
  select(ko,log2FoldChange)%>% 
  tibble::deframe()

#export KOs for exploration on KEGG website
KOs_for_website <- data.frame(KOs = names(DESeq_res_bin_84_KO_list), colours = DESeq_res_bin_84_KO_list) %>% 
  mutate(colours = case_when(colours ==1 ~ "red,black",
                             TRUE ~ "yellow,black")) %>% 
  write.table("KOs_Bin_84.txt", sep="\t", row.names = FALSE, quote = FALSE)

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
             "00400", #Phenylalanine, tyrosine and tryptophan biosynthesis
             "02010",  #ABC transporters
             "04974",  #Protein digestion and absorption
             "00910",  #Nitrogen metabolism
             "00920",  #Sulfur metabolism
             "00071", #Fatty acid degradation
)


#plot in a loop all the pathways of interest for each bin
pathview(gene.data = DESeq_res_bin_84_KO_list, 
         #pathway.id =pathways,
         pathway.id = "01212",
         species = "ko", 
         keys.align = "y", 
         kegg.native = T,
         low = "yellow",
         high = "red",
         res = 300, cex = 0.25,
         out.suffix = "Bin_84",
         kegg.dir="./Figures/KEGG/")



#print session info and clean the workspace
sessionInfo()
rm(list = ls())
gc()