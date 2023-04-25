########################################
#Metaproteome composition analysis
########################################
#load libraries
require(dplyr)
require(tidyr)
require(phyloseq)
require(ggplot2)
require(ggpubr)

source("scripts/extra_functions.R")

#load metaproteome phyloseq object
metaP_runB_merged<- readRDS("data/metaproteome/metaP_runB_merged.rds")

###################
#Plot number of proteins per sample
###################
# number of proteins per sample
prot_per_sample <- estimate_richness(metaP_runB_merged, split = TRUE, measures = "Observed") %>% 
  mutate(Sample_name = row.names(.)) %>% 
  left_join(sample_data(metaP_runB_merged), by = "Sample_name") %>% 
  mutate(Sample= gsub("_.*","", Sample_name))

#plot
total_prot.p <- prot_per_sample %>% 
  ggplot(aes(x = Treatment, y = Observed, fill = Treatment, group = Sample_name)) + 
  facet_wrap(Type~., scales = "free_x") +
  scale_fill_manual(values =c("Inoculum"="gray50",
                              "Cteno-OM"="red", 
                              "Control"="blue"))+
  geom_bar(position = position_dodge(width=1), 
           stat="identity")+
  geom_text(aes(x = Treatment, y = Observed+100, label = Sample),
            position = position_dodge(width=1))+
  ylab("Number of proteins \n")+
  #scale_y_log10()+
  theme_EF+
  theme(legend.position = "bottom")

#save the plot
ggsave("./Figures/Figure_S1-Total_prot.pdf", 
       plot = total_prot.p,
       units = "mm",
       width = 90, height = 90, 
       scale = 3,
       dpi = 300)

#calculate mean and SE for each fraction
prot_per_sample %>% group_by(Treatment, Type) %>% 
  summarise(Min= min(Observed), Max = max(Observed),
            Mean = mean(Observed), SE=se(Observed))

###################
#Protein overlaps between fraction
###################
sample_data(metaP_runB_merged)$Sample= gsub("_.*","", sample_data(metaP_runB_merged)$Sample_name)

y<- list()
frac_overlaps<- data.frame()

for (i in unique(sample_data(metaP_runB_merged)$Sample)){
  sub_group <- subset_samples(metaP_runB_merged, Sample == i) %>% 
    prune_taxa(taxa_sums(.)>0,.)
  for (n in c("Cellular","Exocellular")){
    sub_sample <- subset_samples(sub_group, Type == n) %>% 
      prune_taxa(taxa_sums(.)>0,.)
    y[[paste(i,n, sep ="_")]] <- as.character(row.names(t(otu_table(sub_sample,taxa_are_rows =FALSE))))
  }
}  
overlap.df <- pres_abs_matrix(y) 

#overlap between fractions in each sample
#build an empty datafram
col.names<- c("Sample", "Shared", "Total")
summary <- data.frame(matrix(nrow =0, ncol = 3))
names(summary)<- col.names

#summarize per culture
for (i in c("C1","C2","C3","J1","J2","J3","T0")){
shared<- overlap.df%>% select(contains(i)) %>% 
    mutate(Shared = case_when(rowSums(.) == 2 ~ 1,TRUE~ 0),
           Total = case_when(rowSums(.) > 0 ~ 1,TRUE~ 0)) %>% 
    select(Shared, Total) %>% 
    summarize_all(sum)

summary<- rbind(summary, data.frame(Sample = i, 
                                    Shared = shared$Shared, 
                                    Total = shared$Total))
}

#calculate proportion
summary %>% 
  mutate(Prop = Shared/Total) %>% 
  summarize_all(se)

#summarize shared proteins between all samples
overlap.df%>% mutate(Shared_all = case_when(rowSums(.) == 14 ~ 1,TRUE~ 0),
                     Shared_exo = case_when(rowSums(.[grepl("_Exocellular",names(.))])== 7 ~ 1,
                                                 TRUE~ 0),
                     Shared_Cel = case_when(rowSums(.[grepl("_Cellular",names(.))])== 7 ~ 1,
                                                 TRUE~ 0)) %>% 
                  select(contains("Shared")) %>% 
                  summarize_all(sum)


###################
#Ordination plot
###################
#conduct variance stabiliozation of the metaP dataset
metaP_obj0.ddsMat <- phyloseq_to_deseq2(metaP_runB_merged, ~Treatment)
metaP.DEseq.vsd <- vst(metaP_obj0.ddsMat,  fitType="local", blind=FALSE)

#plot PCA
metaP_pca.df <- plotPCA(metaP.DEseq.vsd, 
                        intgroup=c("Treatment","Type"),
                        returnData=TRUE)%>% 
  mutate(Sample= gsub("_.*","", name))

#extract explained variance
percentVar <- round(100 * attr(metaP_pca.df, "percentVar"))

#plot
metaP_ordination_plot<- ggplot(data = metaP_pca.df,
                               aes(x = PC1, y = PC2, 
                                   colour = Treatment, shape = Type))+
  geom_point(colour = "black", size = 7) +
  geom_point(size = 5) +
  scale_colour_manual(values =c("Inoculum"="gray50",
                              "Cteno-OM"="red", 
                              "Control"="blue"))+
  geom_text(aes(x = PC1, y = PC2,label = Sample),
           nudge_y= -3, size=5, colour = "gray50")+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  #coord_fixed()+
  theme_EF+
  theme(legend.position = "bottom")

#save the plot
ggsave("./Figures/Figure_3-metaP_ordination.pdf", 
       plot = metaP_ordination_plot,
       units = "mm",
       width = 90, height = 90, 
       scale = 3,
       dpi = 300)

#test whether the differences between the runs are significant
df <- colData(metaP.DEseq.vsd)[,c("Treatment","Type")]
sample_distance <- dist(t(assay(metaP.DEseq.vsd)))
adonis_all <- adonis2(sample_distance ~ Treatment*Type, df)
adonis_all

#posthoc to check which ponds are different
groups <- df[["Treatment"]]
mod <- betadisper(sample_distance, groups)
permutest(mod)

#dispersion is different between groups
plot(mod)
boxplot(mod)
mod.HSD <- TukeyHSD(mod)
mod.HSD



###################
#Compare number of proteins between the ML powder and the microcosms
###################
#import sample list
samples_df <- read.csv("./data/metaproteome/Mleidyi_transcriptome_MASCOT_InputFiles.txt",
                       sep="\t", h= T) %>% 
  #mutate(File.Name = case_when(File.Name =="T0_EL_200616145603.raw" ~ "T0_EL (2).raw",
  #                            TRUE ~ File.Name)) %>% 
  mutate(Sample_name= gsub("\\.raw","",basename(File.Name)))


#select only replicate B (in consistency with the rest of the analysis)
samples_df_selected <- samples_df %>% filter(grepl('_2|(2)', Sample_name))%>% 
  pull(Sample_name)

#import PD output
PD_df <- read.csv("data/metaproteome/Mleidyi_transcriptome_MASCOT_Proteins.txt",
                  sep="\t", h= T)%>% 
  
  #select only Master Proteins with unique peptides
  filter(Master=="IsMasterProtein") %>%  
  filter(Number.of.PSMs >=2 , Number.of.Unique.Peptides>=1)%>% 
  
  #replace filenames with sample names
  rename_with(~gsub("Abundance\\.|\\.Sample","",.), everything()) %>% 
  rename_at(all_of(samples_df$File.ID), ~ samples_df$Sample_name) %>% 
  
  #subset the samples of interest
  select(c(samples_df_selected)) %>% 
  
  #summarize the EL and EH fractions for the microcosms and remove the separated fractions
  mutate_all(funs(replace_na(., 0))) %>% 
  mutate(J1 = rowSums(across(starts_with("J1"))),
         J2 = rowSums(across(starts_with("J2"))),
         J3 = rowSums(across(starts_with("J3")))) %>% 
  select(-contains(c("EH","EL"))) %>% 
  
  #summarize the fractions of the powder
  mutate(ML_powder=rowSums(across(contains("_")))) %>% 
  select(ML_powder,starts_with("J")) %>% 
  
  #remove proteins that were not observed in the powder
  filter(ML_powder>0)


#summarize observed proteins per sample and plot
data.frame(Sample = names(PD_df), Obs_prot = colSums(PD_df != 0)) %>% 
  ggplot(aes(x=Sample, y=Obs_prot, label = Obs_prot))+
  geom_col()+ geom_text(nudge_y = 100)+
  ylab("Number of proteins \n")+
  #scale_y_log10()+
  theme_EF

#save the plot
ggsave("./Figures/Figure_S1-Cteno_prot.pdf", 
       plot = last_plot(),
       units = "mm",
       width = 90, height = 90, 
       scale = 3,
       dpi = 300)

###################
#print session info and clean the workspace
###################
sessionInfo()
rm(list = ls())
gc()