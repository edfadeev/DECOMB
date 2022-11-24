########################################
#Metaproteome QC analysis
########################################
#load libraries
require(dplyr)
require(tidyr)
require(phyloseq)
require(ggplot2)
require(ggpubr)
require(rstatix)
require(DESeq2)
require(vegan)

#load metaproteome phyloseq object
metaP_obj0<- readRDS("data/metaproteome/metaP_ps_raw.rds")

###################
#Plot number of proteins per sample
###################
# number of proteins per sample
prot_per_sample <- estimate_richness(metaP_obj0, split = TRUE, measures = "Observed") %>% 
  mutate(Sample_name = row.names(.)) %>% 
  left_join(sample_data(metaP_obj0), by = "Sample_name")

prot_counts_bar.p<- list()

for (frac in c("MP","EH","EL")){
  sub_df<- prot_per_sample %>% dplyr::filter(Fraction == frac)
  
  prot_counts_bar.p[[frac]] <- ggplot(sub_df, aes(x = Sample_name, y = Observed,
                                                  fill = Treatment)) + 
    facet_wrap(Fraction~.) +
    scale_fill_manual(values =c("T0"="gray50",
                                "Jelly"="red",
                                "Control"="blue"))+
    geom_col()+
    ylab("# of proteins \n")+
    #scale_y_log10()+
    theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          text=element_text(size=14),legend.position = "bottom", 
          axis.title.x = element_blank(), axis.text.x = element_text(angle=90))
}

ggarrange(prot_counts_bar.p[["MP"]], prot_counts_bar.p[["EH"]],prot_counts_bar.p[["EL"]],
          ncol = 3, nrow = 1, align = "hv")


#save the plot
ggsave("./Figures/total_prot.pdf", 
       plot = last_plot(),
       units = "cm",
       width = 30, height = 15, 
       #scale = 1,
       dpi = 300)

#test differences between runs
prot_per_sample_test <- prot_per_sample   %>%
  group_by(Fraction) %>% 
  t_test(Observed ~ Run, paired = TRUE, p.adjust.method = "BH") %>%
  add_significance()


###################
#Protein overlaps between replicates
###################
y<- list()
overlaps_table<- data.frame(Sample=character(), Run1=numeric(),Run2=numeric(), Overlap=numeric())

for (i in sample_data(metaP_obj0)$SampleID){
  sub_group <- subset_samples(metaP_obj0, SampleID == i) %>% 
    prune_taxa(taxa_sums(.)>0,.)
  for (n in c("A","B")){
    sub_sample <- subset_samples(sub_group, Run == n)
    sub_sample <- prune_taxa(taxa_sums(sub_sample)>0,sub_sample)
    y[[paste(i,n, sep = "_")]] <- as.character(row.names(otu_table(sub_sample)))
  }
  overlap<- VennDiagram::calculate.overlap(y)
  overlap.df<- data.frame(Sample= paste(i), Run1 = length(overlap$a1),Run2 = length(overlap$a2), Overlap= length(overlap$a3))
  overlaps_table<- rbind(overlaps_table, overlap.df)
  y<- list()
}

overlaps_table<- overlaps_table %>% unique() %>% 
  mutate(Run1_prop = signif(Overlap/Run1, digits = 2),
         Run2_prop = signif(Overlap/Run2, digits = 2),
         Type = case_when(grepl("C",Sample) ==TRUE ~"Control", 
                          grepl("J",Sample) ==TRUE ~"Jelly",
                          grepl("T0",Sample) ==TRUE ~ "T0"))


###################
#Generate PCA plot 
###################
#conduct variance stabiliozation of the metaP dataset
metaP_obj0.ddsMat <- phyloseq_to_deseq2(metaP_obj0, ~Run)
metaP.DEseq.vsd <- vst(metaP_obj0.ddsMat,  fitType="local", blind=FALSE)

#plot PCA
metaP_pca.df <- plotPCA(metaP.DEseq.vsd, 
                        intgroup=c("SampleID","Treatment","Run","Fraction"),
                        returnData=TRUE)

#extract explained variance
percentVar <- round(100 * attr(metaP_pca.df, "percentVar"))

#plot
metaP_ordination_plot<- ggplot(data = metaP_pca.df,
                               aes(x = PC1, y = PC2, colour = Run, shape = Fraction))+
  geom_point(fill = "black", size = 6) +
  geom_point(size = 4,alpha = 0.8) +
  geom_text(aes(x = PC1, y = PC2,label = SampleID), 
            nudge_y= -8,size=4)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_colour_manual(values =c(#"MP"="gray50",
    "A"="red",
    "B"="blue"))+
  coord_fixed()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(size=14),legend.position = "bottom")

#save the plot
ggsave("./Figures/metaP_runs_ord.pdf", 
       plot = metaP_ordination_plot,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)

#test whether the differences between the runs are significant
df <- colData(metaP.DEseq.vsd)[,c("SampleID","Treatment","Run","Fraction")]
sample_distance <- dist(t(assay(metaP.DEseq.vsd)))
adonis_all <- adonis2(sample_distance ~ Treatment+Fraction+Run, df)
adonis_all

#posthoc to check which ponds are different
groups <- df[["Run"]]
mod <- betadisper(sample_distance, groups)
permutest(mod)

#dispersion is different between groups
plot(mod)
boxplot(mod)
mod.HSD <- TukeyHSD(mod)
mod.HSD
plot(mod.HSD)

#based on the fact that the there are potentially sig. differences
#between the runs, I decided to focus only on Run 2, which also consisted
#of more protein observations. Separated tests showed that the differences
#between the treatments are similar in terms of sig. and var. in both replicates

#conduct variance stabiliozation of the metaP dataset
metaP_runB <- subset_samples(metaP_obj0, Run =="B") %>% 
  prune_taxa(taxa_sums(.)>0,.)

#generate metadata
data_runB <- as(sample_data(metaP_runB),"data.frame") %>% 
  mutate(#Fraction = case_when(Fraction %in% c("EL","EH") ~ "exoP",
    #        TRUE ~ "MP"),
    Group = paste(Treatment,Replicate, Fraction, sep ="_"))

sample_data(metaP_runB) <- sample_data(data_runB)

#save the new phyloseq
saveRDS(metaP_runB, "data/metaproteome/metaP_ps_runB.rds")
