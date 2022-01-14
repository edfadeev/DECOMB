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
metaP_obj0<- readRDS("data/metaP_ps_raw.rds")
exoP_obj0<- readRDS("data/exoP_ps_raw.rds")

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

