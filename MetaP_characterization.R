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

source("scripts/extra_functions.R")

#load metaproteome phyloseq object
metaP_obj0<- readRDS("data/metaproteome/metaP_ps_raw.rds")

###################
#Plot number of proteins per sample
###################
# number of proteins per sample
prot_per_sample <- estimate_richness(metaP_obj0, split = TRUE, measures = "Observed") %>% 
  mutate(Sample_name = row.names(.)) %>% 
  left_join(sample_data(metaP_obj0), by = "Sample_name")

total_prot.p <- ggplot(prot_per_sample, aes(x = Sample_name, y = Observed,
                            fill = Treatment)) + 
                        facet_wrap(Fraction~., scales = "free_x") +
                        scale_fill_manual(values =c("T0"="gray50",
                                                    "Jelly"="red",
                                                      "Control"="blue"))+
  geom_col()+
  ylab("# of proteins \n")+
  #scale_y_log10()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(size=20),legend.position = "bottom", 
        axis.title.x = element_blank(), axis.text.x = element_text(angle=90))

#save the plot
ggsave("./Figures/total_prot.pdf", 
       plot = total_prot.p,
       units = "cm",
       width = 30, height = 30, 
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
overlaps_table<- data.frame()

for (i in sample_data(metaP_obj0)$SampleID){
  sub_group <- subset_samples(metaP_obj0, SampleID == i) %>% 
    prune_taxa(taxa_sums(.)>0,.)
  for (n in c("A","B")){
    sub_sample <- subset_samples(sub_group, Run == n)
    sub_sample <- prune_taxa(taxa_sums(sub_sample)>0,sub_sample)
    y[[n]] <- as.character(row.names(otu_table(sub_sample)))
  }
  overlap.df <- pres_abs_matrix(y) %>% 
    mutate(Both = case_when(A==1 & 
                              B==1 ~1, 
                            TRUE~ 0),
           A_u = case_when(A==1 & 
                             B==0 ~1,
                           TRUE~ 0),
           B_u = case_when(A==0 & 
                             B == 1 ~ 1,
                           TRUE~ 0)) %>% 
    summarize_all(sum) %>% 
    mutate(Sample = i)
  overlaps_table<- rbind(overlaps_table, overlap.df)
  y<- list()
}

overlaps_table<- overlaps_table %>% unique() %>% 
  mutate(A_prop = signif(Both/A, digits = 2),
         B_prop = signif(Both/B, digits = 2),
         Type = case_when(grepl("C",Sample) ==TRUE ~"Control", 
                          grepl("J",Sample) ==TRUE ~"Jelly",
                          grepl("T0",Sample) ==TRUE ~ "T0"))

overlaps_table_mean <- overlaps_table %>% 
                        group_by(Type) %>% 
                        summarize_all(c(mean = mean, se= se))

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
  geom_point(colour = "black", size = 5) +
  geom_point(size = 4) +
  geom_text(aes(x = PC1, y = PC2,label = SampleID), 
            nudge_y= -8, size=4, colour = "gray50")+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_colour_manual(values =c(#"MP"="gray50",
    "A"="orange",
    "B"="lightblue"))+
  coord_fixed()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(size=20),legend.position = "bottom")

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
    Sample_name = gsub("_.*","", SampleID),
    Group = paste(Treatment,Replicate, Fraction, sep ="_"))

sample_data(metaP_runB) <- sample_data(data_runB)

#save the new phyloseq
saveRDS(metaP_runB, "data/metaproteome/metaP_ps_runB.rds")

###################
#Protein overlaps between fraction
###################
y<- list()
frac_overlaps<- data.frame()

for (i in sample_data(metaP_runB)$Sample_name){
  sub_group <- subset_samples(metaP_runB, Sample_name == i) %>% 
    prune_taxa(taxa_sums(.)>0,.)
  for (n in c("MP","EH","EL")){
    sub_sample <- subset_samples(sub_group, Fraction == n) %>% 
      prune_taxa(taxa_sums(.)>0,.)
    y[[n]] <- as.character(row.names(otu_table(sub_sample)))
  }
  overlap.df <- pres_abs_matrix(y) %>% 
    mutate(Shared_all = case_when(MP==1 & 
                                    EH==1 & 
                                    EL == 1 ~ 1,
                                  TRUE~ 0),
           Shared_exoP = case_when(MP==0 & 
                                     EH==1 &
                                     EL == 1 ~ 1,
                                   TRUE~ 0),
           MP_u = case_when(MP==1 & 
                              EH==0 & 
                              EL == 0 ~ 1,
                            TRUE~ 0),
           EH_u = case_when(MP==0 & 
                              EH==1 & 
                              EL == 0 ~ 1,
                            TRUE~ 0),
           EL_u = case_when(MP==0 & 
                              EH==0 & 
                              EL == 1 ~ 1,
                            TRUE~ 0)) %>% 
    summarize_all(sum) %>% 
    mutate(Sample = i)
  frac_overlaps<- rbind(frac_overlaps, overlap.df)
  y<- list()
}

frac_overlaps<- frac_overlaps %>% unique() %>% 
  mutate(MP_prop = signif(Shared_all/MP, digits = 2),
         EH_prop = signif(Shared_all/EH, digits = 2),
         EL_prop = signif(Shared_all/EL, digits = 2),
         Type = case_when(grepl("C",Sample) ==TRUE ~"Control", 
                          grepl("J",Sample) ==TRUE ~"Jelly",
                          grepl("T0",Sample) ==TRUE ~ "T0"))
