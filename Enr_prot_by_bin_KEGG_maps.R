require(pheatmap)
require(phyloseq)
require(DESeq2)
require(dplyr)
require(ggplot2)

#import prot. enrichment results 
DESeq_res <- read.csv("data/DESEq_res.csv", row.names = 1)

#import KEGG annoptations of the metagenome and merge with  genes list per bin
metaG_KEGG_modules_by_gene <- read.csv("data/metagenome/spades-KOfam-functions.txt", sep="\t", h= T) %>% 
                                  left_join(read.csv("./data/metagenome/Refined-bins-gene-calls.txt",
                                                     sep="\t", h= T), by = c("gene_callers_id"))

#########################################################
#Merge enriched proteins data with metagenome KEGG IDs###
#########################################################
DESeq_res_bins <- metaG_KEGG_modules_by_gene %>% left_join(DESeq_res[c("gene_callers_id", "Type","Fraction","log2FoldChange")], by = c("gene_callers_id"))


#calculate how many enriched proteins there are per bin
DESeq_res_bins_totals <- DESeq_res_bins %>% 
                          filter(Type %in% c("Jelly","Control") & !is.na(Bin)) %>% 
                          group_by(Bin, Type, Fraction) %>%  summarize (n())



#########################################################
#Explore enr. proteins in Bin_84 (Pseudoalteromonas)  ###
#########################################################
DESeq_res_bin_84 <- DESeq_res_bins %>%  
  filter(Bin =="Bin_84") %>% 
  mutate(log2FoldChange = case_when(is.na(log2FoldChange)==TRUE ~ -1, 
                                    TRUE ~ log2FoldChange)) %>% 
  filter(#Type =="Jelly", 
         !is.na(accession)) %>% 
  select(accession,log2FoldChange)%>% 
  tibble::deframe()

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
pathview(gene.data = DESeq_res_bin_84, 
           pathway.id =pathways,
           species = "ko", 
           keys.align = "y", 
           kegg.native = T,
            #map.null=FALSE,
            #both.dirs = TRUE,
           #discrete	=list(gene=TRUE),
           res = 300, cex = 0.25,
           out.suffix = "Bin_84",
           kegg.dir="./Figures/KEGG/")

