########################################
#Metaproteome data analysis
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

#conduct NSAF transformation
#https://github.com/moldach/proteomics-spectralCount-normalization/blob/master/nsaf.R
#https://rdrr.io/github/DanielSprockett/reltools/man/add_nsaf.html
add_nsaf=function(ps, prot_length){
  if(ps@otu_table@taxa_are_rows == TRUE){
    mat <- (otu_table(ps))
  }else{
    mat <- t((otu_table(ps)))
  }
  prot_len <- unlist(as.numeric(tax_table(ps)[,prot_length])) # Unlist your protein lengths before you sweep
  mat_prop <- sweep(mat,1,prot_len,"/") # Divide spectral counts (SpC) for a protein by its length (L)
  mat_sum <- as.data.frame(colSums(mat_prop)) # Get the column sums for each cell-line/treatment
  mat_sum <- mat_sum[,1]
  mat_nsaf <- sweep(mat_prop,2,mat_sum,"/") # Normalize by dividing by the sum of all SpC/L for all proteins identified 
  otu_table(ps) <- otu_table(mat_nsaf, taxa_are_rows = TRUE)
  return(ps)
}

#load metaproteome phyloseq object
metaP_obj0<- readRDS("metaP/metaP_ps_raw.rds")
exoP_obj0<- readRDS("metaP/exoP_ps_raw.rds")

###################
#Plot number of proteins per sample
###################
#merge the metaP and exoP into a single phyloseq object
all_prot_obj0 <- merge_phyloseq(metaP_obj0, exoP_obj0)

# number of proteins per sample
prot_per_sample <- estimate_richness(all_prot_obj0, split = TRUE, measures = "Observed") %>% 
  mutate(Sample_name = row.names(.)) %>% 
  left_join(sample_data(all_prot_obj0), by = "Sample_name") %>% 
  separate(Sample_name, into = c("Sample_name","Fraction"), sep ="_")

prot_counts_bar.p<- list()
for (frac in c("MP","exoP")){
  sub<- prot_per_sample %>% dplyr::filter(Fraction ==frac)
  prot_counts_bar.p[[frac]] <- ggplot(sub, aes(x = Sample_name, y = Observed, fill = Type)) + 
    facet_wrap(Fraction~.) +
    geom_col()+
    ylab("# of proteins \n")+
    #scale_y_log10()+
    theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          text=element_text(size=14),legend.position = "bottom", 
          axis.title.x = element_blank(), axis.text.x = element_text(angle=90))
}

ggarrange(prot_counts_bar.p[["MP"]], prot_counts_bar.p[["exoP"]],
          ncol = 2, nrow = 1, align = "hv")


#save the plot
ggsave(paste0(wd,"/R_figures/total_prot.png"), 
       plot = last_plot(),
       units = "cm",
       width = 30, height = 15, 
       #scale = 1,
       dpi = 300)

###################
#Barplots
###################
#conduct NSAF transformation
prot_nsaf<- add_nsaf(all_prot_obj0, "prot_length")

prot_nsaf.long <- psmelt(prot_nsaf) %>% 
  separate(Sample_name, sep ="_", into = c("Sample","Fraction"))

#remove below 3% ra
taxa_classes <- unique(prot_nsaf.long$t_genus[!prot_nsaf.long$Abundance<0.001])

prot_nsaf.long$t_genus[prot_nsaf.long$Abundance<0.001] <- "Other taxa"

prot_nsaf.long$t_genus <- factor(prot_nsaf.long$t_genus,
                                 levels=c(taxa_classes,"Other taxa"))

ggplot(prot_nsaf.long, 
       aes(x = Sample, y = Abundance,
           fill = t_genus)) + 
  facet_grid(.~Fraction, space= "fixed") +
  geom_bar(position="stack", stat="identity")+
  #scale_fill_manual(values = phyla.col )+ 
  #guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Sequence proportions (%) \n")+
  geom_hline(aes(yintercept=-Inf)) + 
  geom_vline(aes(xintercept=-Inf)) +
  geom_vline(aes(xintercept=Inf))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),axis.text.x = element_text(angle = 90),
        text=element_text(size=14),legend.position = "bottom", 
        axis.title.x = element_blank())

prot_nsaf.long_sub<- prot_nsaf.long %>%  filter(t_genus %in% c("Pseudoalteromonas", "Alteromonas", "Vibrio", "Synechococcus"))

ggplot(prot_nsaf.long_sub, 
       aes(x = Sample, y = Abundance,
           fill = t_genus)) + 
  facet_grid(.~Fraction, space= "fixed") +
  geom_bar(position="stack", stat="identity")+
  #scale_fill_manual(values = phyla.col )+ 
  #guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Sequence proportions (%) \n")+
  geom_hline(aes(yintercept=-Inf)) + 
  geom_vline(aes(xintercept=-Inf)) +
  geom_vline(aes(xintercept=Inf))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),axis.text.x = element_text(angle = 90),
        text=element_text(size=14),legend.position = "bottom", 
        axis.title.x = element_blank())



###################
#Generate PCA plot
###################
#stabilize the dataset using gemetric mean 
# calculate geometric means prior to estimate size factors
prot_frac_agg.dds <- phyloseq_to_deseq2(metaP_obj0, ~1)
prot_frac_agg.dds = estimateSizeFactors(prot_frac_agg.dds)
prot_frac_agg.dds <- estimateDispersions(prot_frac_agg.dds, fitType='local')
prot.vst <- getVarianceStabilizedData(prot_frac_agg.dds)

metaP_obj.vst<-metaP_obj0
otu_table(metaP_obj.vst)<- otu_table(prot.vst, taxa_are_rows = TRUE)
rm(prot.vst, prot_frac_agg.dds)

#produce ordination
metaP_dist <- phyloseq::distance(metaP_obj.vst, method = "euclidean")
metaP_pca <- ordinate(metaP_obj.vst, method = "RDA",metaP_dist )
metaP_pca.df <- plot_ordination(metaP_obj.vst, metaP_pca, axes = c(1,2,3),justDF = TRUE) 

#extract explained variance
metaP_pca.evals <- 100 * summary(metaP_pca)$cont$importance[2, c("PC1","PC2")]

metaP_ordination_plot<- ggplot(data = metaP_pca.df, aes(x = PC1, y = PC2))+
  geom_point(fill = "black", size = 6,alpha = 0.8) +
  geom_point(aes(colour = Type), size = 4,alpha = 0.8) +
  geom_text(aes(x = PC1, y = PC2,label = Sample_name), 
            nudge_y= -1,size=4)+
  labs(x = sprintf("PC1 [%s%%]", round(metaP_pca.evals[1], 2)), 
       y = sprintf("PC2 [%s%%]", round(metaP_pca.evals[2], 2)))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(size=14),legend.position = "bottom")


#test significance of clustering
df <- as(sample_data(metaP_obj.vst), "data.frame")
d <- phyloseq::distance(metaP_obj.vst, "euclidean")
adonis_all <- adonis2(d ~ Type  , df)
adonis_all

#posthoc to check which ponds are different
groups <- df[["Type"]]
mod <- betadisper(d, groups)
permutest(mod)

#dispersion is different between groups
plot(mod)
boxplot(mod)
mod.HSD <- TukeyHSD(mod)
mod.HSD
plot(mod.HSD)


#stabilize the dataset using gemetric mean 
# calculate geometric means prior to estimate size factors
prot_frac_agg.dds <- phyloseq_to_deseq2(exoP_obj0, ~1)
prot_frac_agg.dds = estimateSizeFactors(prot_frac_agg.dds)
prot_frac_agg.dds <- estimateDispersions(prot_frac_agg.dds, fitType='local')
prot.vst <- getVarianceStabilizedData(prot_frac_agg.dds)

exoP_obj.vst<-exoP_obj0
otu_table(exoP_obj.vst)<- otu_table(prot.vst, taxa_are_rows = TRUE)
rm(prot.vst, prot_frac_agg.dds)

#produce ordination
exoP_pca <- ordinate(exoP_obj.vst, method = "RDA", distance = "euclidean")
exoP_pca.df <- plot_ordination(exoP_obj.vst, exoP_pca, axes = c(1,2,3),justDF = TRUE) 

#extract explained variance
exoP_pca.evals <- 100 * summary(exoP_pca)$cont$importance[2, c("PC1","PC2")]

exoP_ordination_plot<- ggplot(data = exoP_pca.df, aes(x = PC1, y = PC2))+
  geom_point(fill = "black", size = 6,alpha = 0.8) +
  geom_point(aes(colour = Type), size = 4,alpha = 0.8) +
  geom_text(aes(x = PC1, y = PC2,label = Sample_name), 
            nudge_y= -1,size=4)+
  labs(x = sprintf("PC1 [%s%%]", round(exoP_pca.evals[1], 2)), 
       y = sprintf("PC2 [%s%%]", round(exoP_pca.evals[2], 2)))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(size=14),legend.position = "bottom")


#test significance of clustering
df <- as(sample_data(exoP_obj.vst), "data.frame")
d <- phyloseq::distance(exoP_obj.vst, "euclidean")
adonis_all <- adonis2(d ~ Type  , df)
adonis_all

#posthoc to check which ponds are different
groups <- df[["Type"]]
mod <- betadisper(d, groups)
permutest(mod)

#dispersion is different between groups
plot(mod)
boxplot(mod)
mod.HSD <- TukeyHSD(mod)
mod.HSD
plot(mod.HSD)


#plot both PCAs
ggarrange(metaP_ordination_plot, exoP_ordination_plot, rows = 1,
          cols = 2, align = "hv", labels = c("MetaP","ExoP"))

#save the plot
ggsave(paste0(wd,"/R_figures/PCAs.pdf"), 
       plot = last_plot(),
       units = "cm",
       width = 30, height = 15, 
       #scale = 1,
       dpi = 300)

###################
#Identify enriched proteins in MetaP
###################
metaP_obj0_no_T0<- subset_samples(metaP_obj0, Type != "T0")
#remove proteins that were not observed 
metaP_obj0_no_T0<- prune_taxa(taxa_sums(metaP_obj0_no_T0)>0,metaP_obj0_no_T0)

metaP_obj0.ddsMat <- phyloseq_to_deseq2(metaP_obj0_no_T0, ~Type)
metaP_obj0.ddsMat = estimateSizeFactors(metaP_obj0.ddsMat)
metaP_obj0.ddsMat <- estimateDispersions(metaP_obj0.ddsMat)
metaP.DEseq <- DESeq(metaP_obj0.ddsMat, fitType="local")
metaP.DEseq.res <- results(metaP.DEseq)

#extract only significant proteins
metaP.DEseq.res.sig <- as(metaP.DEseq.res, "data.frame") %>%  filter(padj < 0.1)

metaP.DEseq.res.sig <- cbind(as(metaP.DEseq.res.sig, "data.frame"),
                                           as(tax_table(metaP_obj0)[rownames(metaP.DEseq.res.sig), ], "matrix"))

###################
#Import metabolism estimates for the bins
###################
modules_info <- read.csv(paste(wd,"metaG_analysis/metaG_anvio/07_METABOLISM/modules_info.txt",sep=""), sep="\t", h= T)

Bins_kofam_hits <- read.csv(paste(wd,"metaG_analysis/metaG_anvio/07_METABOLISM/Bins_kofam_hits.txt",sep=""), sep="\t", h= T)

Bins_modules <- read.csv(paste(wd,"metaG_analysis/metaG_anvio/07_METABOLISM/Bins_modules.txt",sep=""), sep="\t", h= T)


#merge the identified modules and the genes calls
Bins_gene_calls_KEGG_modules<-Bins_modules %>% 
  tidyr::separate_rows(gene_caller_ids_in_module, sep = ',') %>% 
  dplyr::rename("gene_caller_id" = "gene_caller_ids_in_module") %>% 
  mutate(gene_caller_id=as.integer(gene_caller_id),
         genome_name = factor(genome_name, levels = c("Bin_84_1","Bin_76_1","Bin_5_2",
                                                      "Bin_2_1","Bin_5_3","Bin_38_1","Bin_2_2",
                                                      "Bin_179_1","Bin_115_2","Bin_115_1","Bin_102_1",
                                                      "Bin_12_1","Bin_150_1_1"))) %>% 
  left_join(Bins_kofam_hits_parsed[,c("genome_name","gene_caller_id","contig","kegg_module","ko","ko_definition")],
            by = c("genome_name","kegg_module","gene_caller_id"))


#merge the enriched proteins with the bins
metaP.DEseq.res.sig <- metaP.DEseq.res.sig  %>%
  mutate(gene_caller_id = as.integer(gene_caller_id)) %>% 
  left_join(., Bins_gene_calls_KEGG_modules, by = "gene_caller_id") %>% 
  select(-c("start","stop", "direction", "call_type", "source","unique_id","version","contig.y"))


write.table(metaP.DEseq.res.sig, "metaP/metaP_DEseq_res_sig.txt")

#see how many enriched proteins were associated with bins
metaP.DEseq.res.sig %>% group_by(genome_name) %>% 
  dplyr::summarize(num_of_prot= n()) %>% 
  dplyr::arrange(desc(num_of_prot))


#subset each bins proteins and summarize according to KO modules 
#Bin 84_1 - Pseudoalteromonas phenolica - length 3.92Mbp (C93/R0)
metaP.DEseq.res.sig_Bin_P_phenolica<-metaP.DEseq.res.sig %>% filter(genome_name == "Bin_84_1") 


Bins_gene_calls_KEGG_modules_Bin_84_1<- Bins_gene_calls_KEGG_modules %>% 
  mutate(gene_caller_id = as.character(gene_caller_id)) %>% 
  left_join(metaP.DEseq.res.sig, by = "gene_caller_id") %>% 
  filter(genome_name == "Bin_84_1") 

#export for KEGG website
Bins_gene_calls_KEGG_modules_Bin_84_1_for_KEGG <- Bins_gene_calls_KEGG_modules_Bin_84_1 %>%  
  mutate(colour = case_when(log2FoldChange>0 ~"red", 
                            log2FoldChange<0 ~"blue",
                            is.na(log2FoldChange)~"orange")) %>% 
  select("ko","colour")
  
write.table(Bins_gene_calls_KEGG_modules_Bin_84_1_for_KEGG, "Bin_84_KEGG.txt",
            row.names = FALSE,
            col.names = FALSE,quote = FALSE)











metaP.DEseq.res.sig_Bin_84_1<-metaP.DEseq.res.sig %>% filter(genome_name == "Bin_84_1")

Bin_84_1_logFC <- metaP.DEseq.res.sig_Bin_84_1$log2FoldChange
names(Bin_84_1_logFC)<- metaP.DEseq.res.sig_Bin_84_1$ko


pv.out <- pathview(gene.data = Bin_84_1_logFC, 
                   pathway.id ="01230",
                   species = "hsa",
                   kegg.native = T,
                   keys.align = "y", 
                   gene.idtype="KEGG",
                   low = "blue", mid = "gray", high = "#ffa500", bin = 30,
                   out.suffix = "Bin_84_1")




#explore metabolic capacities
test <- Bins_modules %>% filter(module_completeness > 0.5) %>% group_by(genome_name, module_subcategory) %>% 
  summarize(total= n())


ggplot(test, aes(x = genome_name, y = total))+
  geom_col()+
  facet_wrap(module_subcategory~., scales = "free_y")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(size=14),legend.position = "bottom", 
        axis.title.x = element_blank(), axis.text.x = element_text(angle=90))





#subset each bins proteins and summarize according to KO modules 
#Bin 84_1 - Pseudoalteromonas phenolica - length 3.92Mbp (C93/R0)
metaP.DEseq.res.sig_Bin_P_phenolica<-metaP.DEseq.res.sig %>% filter(genome_name == "Bin_84_1") 

#check how many proteins in each KO module
metaP.DEseq.res.sig_Bin_P_phenolica %>% 
  separate_rows(modules_with_ko, sep = ',') %>% 
  group_by(modules_with_ko) %>% summarize(num_of_prot= n()) %>% 
  arrange(desc(num_of_prot))

#check for potentially duplicated rows as a result of different potential KO annotations (will be resolved at a later stage)
metaP.DEseq.res.sig_Bin_P_phenolica[duplicated(metaP.DEseq.res.sig_Bin_P_phenolica$start),]




metaP.DEseq.res.sig_Bin_Vibrio<-metaP.DEseq.res.sig %>% filter(genome_name == "Bin_38_1") 

metaP.DEseq.res.sig_Bin_Bermanella<-metaP.DEseq.res.sig %>% filter(genome_name == "Bin_115_1") 

metaP.DEseq.res.sig_Bin_Glaciecola<-metaP.DEseq.res.sig %>% filter(genome_name == "Bin_115_2") 

metaP.DEseq.res.sig_Bin_Alteromonas<-metaP.DEseq.res.sig %>% filter(genome_name %in% c("Bin_2_1","Bin_2_2"))

#not assigned to any bin
metaP.DEseq.res.sig_no_bin<-metaP.DEseq.res.sig %>% filter(is.na(genome_name) == TRUE) 


#split posible annotations into separated rows
test <- prot_frac_agg_metaP.DEseq.res.sig_no_bin 
