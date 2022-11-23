########################################
#Metaproteome data analysis
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
metaP_obj0<- readRDS("data/metaproteome/metaP_ps_raw.rds")

sample_data(metaP_obj0)$Fraction<- gsub(" ","", sample_data(metaP_obj0)$Fraction)



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

ggarrange(prot_counts_bar.p[["MP"]], prot_counts_bar.p[["EH"]],prot_counts_bar.p[["EL"]],
          ncol = 3, nrow = 1, align = "hv")


#save the plot
ggsave(paste0(wd,"./Figures/total_prot.pdf"), 
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
  mutate(Run1_prop = signif(100*Overlap/Run1, digits = 0),
         Run2_prop = signif(100*Overlap/Run2, digits = 0),
         Type = case_when(grepl("C",Sample) ==TRUE ~"Control", 
                          grepl("J",Sample) ==TRUE ~"Jelly",
                          grepl("T0",Sample) ==TRUE ~ "T0"))

######################################
#Sum the two technical replicates and merge the EL and EH fractions
######################################
#merge technical replicates
metaP_tech_agg <- merge_samples(metaP_obj0, "SampleID", fun = sum)

#generate aggregated metadata
meta_agg <- as(sample_data(metaP_obj0),"data.frame") %>% 
            select(SampleID, Replicate, Fraction, Treatment) %>% 
            unique() %>% 
            mutate(Fraction = case_when(Fraction %in% c("EL","EH") ~ "exoP",
                                        TRUE ~ "MP"),
                   Sample = paste(Treatment,Replicate, Fraction, sep ="_"))

sample_data(metaP_tech_agg)<- sample_data(meta_agg)

#merge exoP fractions
metaP_agg <- merge_samples(metaP_tech_agg, "Sample", fun = sum)

#generate aggregated metadata
meta_agg <- as(sample_data(metaP_tech_agg),"data.frame") %>% 
  select(Sample, Fraction, Replicate, Treatment) %>% 
  unique() %>% 
  mutate(Group = paste(Treatment, Replicate, sep ="_"))

rownames(meta_agg) <- meta_agg$Sample 

sample_data(metaP_agg)<- sample_data(meta_agg)

###################
#Summary of merged samples
###################
prot_per_sample_agg <- estimate_richness(metaP_agg, split = TRUE, measures = "Observed") %>% 
  mutate(Sample = row.names(.)) %>% 
  left_join(sample_data(metaP_agg), by = "Sample")

###################
#Taxonomic compositions
###################
#conduct NSAF transformation
metaP_agg_nsaf<- add_nsaf(metaP_agg, "prot_length")

#melt phyloseq into a dataframe for ploting
prot_nsaf.long <- psmelt(metaP_agg_nsaf) 

#sum up each taxa
prot_nsaf.class.agg <- prot_nsaf.long %>% group_by(Fraction, Group, Class, Order) %>% 
  summarize(Tot.abundance = sum(Abundance)) %>% 
  mutate(Fraction = factor(Fraction, levels =c("MP","exoP")))

#remove below 1% ra
taxa_classes <- sort(as.character(unique(prot_nsaf.class.agg$Class[!prot_nsaf.class.agg$Tot.abundance<0.01])))

prot_nsaf.class.agg$Class[prot_nsaf.class.agg$Tot.abundance<0.01] <- "Other taxa"
prot_nsaf.class.agg$Class[is.na(prot_nsaf.class.agg$Class)] <- "Other taxa"

prot_nsaf.class.agg$Class <- factor(prot_nsaf.class.agg$Class,
                                      levels=c(taxa_classes,"Other taxa"))

prot_tax_comp.p<- ggplot(prot_nsaf.class.agg, 
       aes(x = Group, y = Tot.abundance,
           fill = Class)) + 
  facet_grid(.~Fraction, space= "fixed") +
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values = class_col)+ 
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
ggsave("./Figures/metaP_tax_order_comp.pdf", 
       plot = prot_tax_comp.p,
       units = "cm",
       width = 30, height = 15, 
       #scale = 1,
       dpi = 300)

###################
#Functional compositions
###################
#sum up each COG
prot_nsaf.COG.agg <- prot_nsaf.long %>% group_by(Fraction, Group, COG20_CATEGORY_function) %>% 
  mutate(COG20_CATEGORY_function = case_when(is.na(COG20_CATEGORY_function) ~ "Unk",
                                             TRUE ~ COG20_CATEGORY_function)) %>% 
  summarize(Tot.abundance = sum(Abundance)) %>% 
  mutate(Fraction = factor(Fraction, levels =c("MP","exoP")))

#remove below 1% ra
COG_categories <- sort(as.character(unique(prot_nsaf.COG.agg$COG20_CATEGORY_function[!prot_nsaf.COG.agg$Tot.abundance<0.02])))

prot_nsaf.COG.agg$COG20_CATEGORY_function[prot_nsaf.COG.agg$Tot.abundance<0.02] <- "Other"
prot_nsaf.COG.agg$COG20_CATEGORY_function[is.na(prot_nsaf.COG.agg$COG20_CATEGORY_function)] <- "Other"

prot_nsaf.COG.agg$COG20_CATEGORY_function <- factor(prot_nsaf.COG.agg$COG20_CATEGORY_function,
                                    levels=c(COG_categories,"Other"))

#large colours range
require(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


prot_COG.p<- ggplot(prot_nsaf.COG.agg, 
                         aes(x = Group, y = Tot.abundance,
                             fill = COG20_CATEGORY_function)) + 
  facet_grid(.~Fraction, space= "fixed") +
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values = col_vector)+ 
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
ggsave("./Figures/metaP_COG20_comp.pdf", 
       plot = prot_COG.p,
       units = "cm",
       width = 30, height = 15, 
       #scale = 1,
       dpi = 300)





#print session info and clean the workspace
sessionInfo()
rm(list = ls())
gc()

