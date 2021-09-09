<<<<<<< HEAD
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

tol21rainbow<- c("#771155", "#AA4488","#CC99BB","#114477", 
                 "#4477AA","#117744","#117777","#88CCAA", 
                 "#77CCCC","#00ffff","#44AA77","#44AAAA", 
                 "#777711","#AAAA44","#DDDD77","#774411", 
                 "#AA7744","#DDAA77","#771122","#AA4455", 
                 "#DD7788")

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
  prot_counts_bar.p[[frac]] <- ggplot(sub, aes(x = Sample_name, y = Observed,
                                               fill = Type)) + 
    facet_wrap(Fraction~.) +
    scale_fill_manual(values =c("T0"="gray50",
                                "Jelly"="#019360",
                                "Control"="#A5081A"))+
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
ggsave(paste0(wd,"/R_figures/total_prot.pdf"), 
       plot = last_plot(),
       units = "cm",
       width = 30, height = 15, 
       #scale = 1,
       dpi = 300)

###################
#Barplots
###################
#conduct NSAF transformation
metaP_obj0_nsaf<- add_nsaf(metaP_obj0, "prot_length")
exoP_obj0_nsaf<- add_nsaf(exoP_obj0, "prot_length")

all_prot_obj0_nsaf<- merge_phyloseq(metaP_obj0_nsaf,exoP_obj0_nsaf)

#melt phyloseq into a dataframe for ploting
prot_nsaf.long <- psmelt(all_prot_obj0_nsaf) %>% 
  separate(Sample_name, sep ="_", into = c("Sample","Fraction"))

#sum up each taxa
prot_nsaf.class.agg <- prot_nsaf.long %>% group_by(Fraction,Sample, Class, Order) %>% 
  summarize(Tot.abundance = sum(Abundance)) %>% 
  mutate(Fraction = factor(Fraction, levels =c("MP","exoP")))

#remove below 1% ra
taxa_classes <- sort(as.character(unique(prot_nsaf.class.agg$Order[!prot_nsaf.class.agg$Tot.abundance<0.01])))

prot_nsaf.class.agg$Order[prot_nsaf.class.agg$Tot.abundance<0.01] <- "Other taxa"
prot_nsaf.class.agg$Order[is.na(prot_nsaf.class.agg$Order)] <- "Other taxa"

prot_nsaf.class.agg$Order <- factor(prot_nsaf.class.agg$Order,
                                      levels=c(taxa_classes,"Other taxa"))

prot_tax_comp.p<- ggplot(prot_nsaf.class.agg, 
       aes(x = Sample, y = Tot.abundance,
           fill = Order)) + 
  facet_grid(.~Fraction, space= "fixed") +
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values = tol21rainbow)+ 
  #guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Protein proportions (>1%) \n")+
  geom_hline(aes(yintercept=-Inf)) + 
  geom_vline(aes(xintercept=-Inf)) +
  geom_vline(aes(xintercept=Inf))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),axis.text.x = element_text(angle = 90),
        text=element_text(size=14),legend.position = "bottom", 
        axis.title.x = element_blank())


#save the plot
ggsave(paste0(wd,"/R_figures/metaP_tax_order_comp.pdf"), 
       plot = prot_tax_comp.p,
       units = "cm",
       width = 30, height = 15, 
       #scale = 1,
       dpi = 300)


###################
#Generate PCA plots
###################
#metaproteome
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
metaP_pca.df <- plot_ordination(metaP_obj.vst, metaP_pca, axes = c(1,2,3),justDF = TRUE) %>% 
                  separate(Sample_name, into = c("Sample_name","Fraction"), sep ="_")
#extract explained variance
metaP_pca.evals <- 100 * summary(metaP_pca)$cont$importance[2, c("PC1","PC2")]

metaP_ordination_plot<- ggplot(data = metaP_pca.df, aes(x = PC1, y = PC2))+
  geom_point(fill = "black", size = 6,alpha = 0.8) +
  geom_point(aes(colour = Type), size = 4,alpha = 0.8) +
  geom_text(aes(x = PC1, y = PC2,label = Sample_name), 
            nudge_y= -1.5,size=4)+
  labs(x = sprintf("PC1 [%s%%]", round(metaP_pca.evals[1], 2)), 
       y = sprintf("PC2 [%s%%]", round(metaP_pca.evals[2], 2)))+
  scale_colour_manual(values =c("T0"="gray50",
                              "Jelly"="#019360",
                              "Control"="#A5081A"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(size=14),legend.position = "bottom")

#exoproteome
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
exoP_pca.df <- plot_ordination(exoP_obj.vst, exoP_pca, axes = c(1,2,3),justDF = TRUE) %>% 
  separate(Sample_name, into = c("Sample_name","Fraction"), sep ="_")

#extract explained variance
exoP_pca.evals <- 100 * summary(exoP_pca)$cont$importance[2, c("PC1","PC2")]

exoP_ordination_plot<- ggplot(data = exoP_pca.df, aes(x = PC1, y = PC2))+
  geom_point(fill = "black", size = 6,alpha = 0.8) +
  geom_point(aes(colour = Type), size = 4,alpha = 0.8) +
  geom_text(aes(x = PC1, y = PC2,label = Sample_name), 
            nudge_y= -1,size=4)+
  labs(x = sprintf("PC1 [%s%%]", round(exoP_pca.evals[1], 2)), 
       y = sprintf("PC2 [%s%%]", round(exoP_pca.evals[2], 2)))+
  scale_colour_manual(values =c("T0"="gray50",
                              "Jelly"="#019360",
                              "Control"="#A5081A"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(size=14),legend.position = "bottom")




#plot both PCAs
ggarrange(metaP_ordination_plot, exoP_ordination_plot,align = "hv", 
          labels = c("Metaproteome","Exoproteome"))

#save the plot
ggsave(paste0(wd,"/R_figures/PCAs.pdf"), 
       plot = last_plot(),
       units = "cm",
       width = 30, height = 15, 
       #scale = 1,
       dpi = 300)


###################
#Test significance of separation
###################
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


=======
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
require(pheatmap)

tol21rainbow<- c("#771155", "#AA4488","#CC99BB","#114477", 
                 "#4477AA","#117744","#117777","#88CCAA", 
                 "#77CCCC","#00ffff","#44AA77","#44AAAA", 
                 "#777711","#AAAA44","#DDDD77","#774411", 
                 "#AA7744","#DDAA77","#771122","#AA4455", 
                 "#DD7788")

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
  prot_counts_bar.p[[frac]] <- ggplot(sub, aes(x = Sample_name, y = Observed,
                                               fill = Type)) + 
    facet_wrap(Fraction~.) +
    scale_fill_manual(values =c("T0"="gray50",
                                "Jelly"="#019360",
                                "Control"="#A5081A"))+
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
ggsave(paste0(wd,"/R_figures/total_prot.pdf"), 
       plot = last_plot(),
       units = "cm",
       width = 30, height = 15, 
       #scale = 1,
       dpi = 300)

###################
#Taxonomic compositions
###################
#conduct NSAF transformation
metaP_obj0_nsaf<- add_nsaf(metaP_obj0, "prot_length")
exoP_obj0_nsaf<- add_nsaf(exoP_obj0, "prot_length")

all_prot_obj0_nsaf<- merge_phyloseq(metaP_obj0_nsaf,exoP_obj0_nsaf)

#melt phyloseq into a dataframe for ploting
prot_nsaf.long <- psmelt(all_prot_obj0_nsaf) %>% 
  separate(Sample_name, sep ="_", into = c("Sample","Fraction"))

#sum up each taxa
prot_nsaf.class.agg <- prot_nsaf.long %>% group_by(Fraction,Sample, Class, Order) %>% 
  summarize(Tot.abundance = sum(Abundance)) %>% 
  mutate(Fraction = factor(Fraction, levels =c("MP","exoP")))

#remove below 1% ra
taxa_classes <- sort(as.character(unique(prot_nsaf.class.agg$Order[!prot_nsaf.class.agg$Tot.abundance<0.01])))

prot_nsaf.class.agg$Order[prot_nsaf.class.agg$Tot.abundance<0.01] <- "Other taxa"
prot_nsaf.class.agg$Order[is.na(prot_nsaf.class.agg$Order)] <- "Other taxa"

prot_nsaf.class.agg$Order <- factor(prot_nsaf.class.agg$Order,
                                      levels=c(taxa_classes,"Other taxa"))

prot_tax_comp.p<- ggplot(prot_nsaf.class.agg, 
       aes(x = Sample, y = Tot.abundance,
           fill = Order)) + 
  facet_grid(.~Fraction, space= "fixed") +
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values = tol21rainbow)+ 
  #guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Protein proportions (>1%) \n")+
  geom_hline(aes(yintercept=-Inf)) + 
  geom_vline(aes(xintercept=-Inf)) +
  geom_vline(aes(xintercept=Inf))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),axis.text.x = element_text(angle = 90),
        text=element_text(size=14),legend.position = "bottom", 
        axis.title.x = element_blank())


#save the plot
ggsave(paste0(wd,"/R_figures/metaP_tax_order_comp.pdf"), 
       plot = prot_tax_comp.p,
       units = "cm",
       width = 30, height = 15, 
       #scale = 1,
       dpi = 300)


#print session info and clean the workspace
sessionInfo()
rm(list = ls())
gc()
>>>>>>> 36ce85cc62f48d8a14ac55cd84210dcf25ab8082
