require(pheatmap)
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
#subset to the genes observed also in proteins
metaP_in_bins<- readRDS("data/metaproteome/metaP_ps_runB.rds") %>% 
                add_nsaf(., "prot_length") %>% 
                subset_taxa(., !is.na(Bin)) %>% 
                prune_taxa(taxa_sums(.)>0, .)

#melt phyloseq into a dataframe for ploting
metaP_in_bins.long <- psmelt(metaP_in_bins) 

#sum up each taxa
prot_nsaf.Bin.agg <- metaP_in_bins.long %>% group_by(SampleID, Fraction, Group, Bin) %>% 
  summarize(Tot.abundance = sum(Abundance)) %>% 
  mutate(Fraction = factor(Fraction, levels =c("MP","EH","EL")))

#remove below 1% ra
Bins <- sort(as.character(unique(prot_nsaf.Bin.agg$Bin[!prot_nsaf.Bin.agg$Tot.abundance<0.005])))

prot_nsaf.Bin.agg$Bin[prot_nsaf.Bin.agg$Tot.abundance<0.005] <- "Other bins"
prot_nsaf.Bin.agg$Bin[is.na(prot_nsaf.Bin.agg$Bin)] <- "Other bins"

prot_nsaf.Bin.agg$Bin <- factor(prot_nsaf.Bin.agg$Bin,
                                    levels=c(Bins,"Other bins"))

prot_Bin.p<- ggplot(prot_nsaf.Bin.agg, 
                         aes(x = SampleID, y = Tot.abundance, fill = Bin)) + 
  facet_wrap(.~Fraction, scales = "free_x") +
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values = tol21rainbow)+ 
  #guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Protein proportions (>0.5%) \n")+
  geom_hline(aes(yintercept=-Inf)) + 
  geom_vline(aes(xintercept=-Inf)) +
  geom_vline(aes(xintercept=Inf))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),axis.text.x = element_text(angle = 90),
        text=element_text(size=14),legend.position = "bottom", 
        axis.title.x = element_blank())


#save the plot
ggsave("./Figures/Proteins_by_bin.png", 
       plot = prot_Bin.p,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)


#########################################################
#Explore sig. enriched proteins in binned genes       ###
#########################################################
#import prot. enrichment results 
DESeq_res <- read.csv("data/DESEq_res.csv", row.names = 1)

Kofam_hits <- read.csv("./data/metagenome/spades-Kofam_hits.txt", sep="\t", h= T) %>% 
  dplyr::rename("gene_callers_id" = "gene_caller_id", "contig" = "contig_name")

#import KEGG annoptations of the metagenome and merge with  genes list per bin
metaG_KEGG_modules_by_gene <- read.csv("data/metagenome/spades-KOfam-functions.txt", sep="\t", h= T) %>% 
  left_join(read.csv("./data/metagenome/Refined-bins-gene-calls.txt",
                     sep="\t", h= T), by = c("gene_callers_id")) %>% 
  left_join(Kofam_hits, by = c("gene_callers_id","contig"))


DESeq_res_bins <- metaG_KEGG_modules_by_gene %>% 
  left_join(DESeq_res[c("gene_callers_id", "Type","Fraction","log2FoldChange")], by = c("gene_callers_id")) %>% 
  filter(!is.na(Bin))


#calculate how many enriched proteins there are per bin
DESeq_res_bins_totals <- DESeq_res_bins %>% 
                          filter(Type %in% c("Jelly","Control") & !is.na(Bin)) %>% 
                          group_by(Bin, Type, Fraction) %>%  summarize (n())



#########################################################
#Explore enr. proteins in Bin_84 (Pseudoalteromonas)  ###
#########################################################
#explore manualy
DESeq_res_bin_84_KEGG_modules <- DESeq_res_bins %>%  
  filter(Bin =="Bin_84") %>% 
  filter(!is.na(log2FoldChange))

#produce list of KOs for mapping
DESeq_res_bin_84_KO_list <- DESeq_res_bins %>%  
  filter(Bin =="Bin_84") %>% 
  mutate(log2FoldChange = case_when(is.na(log2FoldChange)==TRUE ~ -1, 
                                    TRUE ~ log2FoldChange)) %>% 
  filter(!is.na(accession))
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
           pathway.id ="00190",
           species = "ko", 
           keys.align = "y", 
           kegg.native = T,
            #map.null=FALSE,
            #both.dirs = TRUE,
           #discrete	=list(gene=TRUE),
           res = 300, cex = 0.25,
           out.suffix = "Bin_84",
           kegg.dir="./Figures/KEGG/")

