########################################
#import metaproteome data
########################################
#import sample list of metaproteomes
prot_sample <- data.frame(Sample_name = c("C1_MP","C2_MP","C3_MP","J1_MP","J2_MP","J3_MP","T0_MP"),
                          Type = c(rep("Control",3),rep("Jelly",3),"T0"),
                          Replicate = c(1:3, 1:3, 1), 
                          row.names = c("C1_MP","C2_MP","C3_MP","J1_MP","J2_MP","J3_MP","T0_MP"))
#import metaP data
protein_raw <- read.csv(paste(wd,"metaP_analysis/metaP-fractions/DECOMB-frac-metaP-concensus_Proteins.txt", sep=""),sep="\t", h= T)

sampleIDS<- c("F2","F4","F6","F8","F10","F13","F15")
samplenames<- c("C1_MP","C2_MP","C3_MP","J1_MP","J2_MP","J3_MP","T0_MP")

#rename sample columns and filter out low confidence proteins and replace NA with 0
protein_filt <- protein_raw %>% 
  dplyr::rename(gene_caller_id = Accession) %>% 
  select_at(vars(!contains("Found"))) %>% 
  rename_with(~gsub("Abundance\\.|\\.Sample","",.), everything()) %>% 
  rename_at(all_of(sampleIDS), ~ samplenames) %>% 
  filter(Number.of.PSMs >=2 , Number.of.Unique.Peptides>=1)%>% 
  mutate_if(is.numeric, funs(replace_na(., 0))) %>% 
  mutate_if(is.numeric,as.integer)

########################################
#merge the dataset into a phyloseq object for convenient data management
########################################
#protein counts as otu_table
protein_counts<- protein_filt %>% select(c("gene_caller_id", samplenames))
protein_counts<- otu_table(data.frame(protein_counts[, all_of(samplenames)], row.names = protein_counts$gene_caller_id), taxa_are_rows = TRUE)

#gene ids as taxonomy table
annotation<- tax_table(as.matrix(data.frame(genes_meta, row.names = genes_meta$gene_caller_id)))

#metadata
sample_meta <- sample_data(prot_sample)

#merge into phyloseq object
prot_obj0<- phyloseq(protein_counts, annotation, sample_meta)

#remove proteins that were not observed 
prot_obj0<- prune_taxa(taxa_sums(prot_obj0)>0,prot_obj0)

########################################
#import exoproteome data
########################################
#import sample list of metaproteomes
exoP_sample <- data.frame(Sample_name = c("C1_exoP","C2_exoP","C3_exoP","J1_exoP","J2_exoP","J3_exoP","T0_exoP"),
                          Type = c(rep("Control",3),rep("Jelly",3),"T0"),
                          Replicate = c(1:3, 1:3, 1), 
                          row.names = c("C1_exoP","C2_exoP","C3_exoP","J1_exoP","J2_exoP","J3_exoP","T0_exoP"))
#import ExoP data
exoP_raw <- read.csv(paste(wd,"metaP_analysis/exoP-fractions/DECOMB-frac-exoP-concensus_Proteins.txt", sep=""),sep="\t", h= T)

exoP_sampleIDS<- c("F1","F2","F3","F4","F5","F6","F7")

#rename sample columns and filter out low confidence proteins and replace NA with 0
exoP_filt <- exoP_raw %>% 
  dplyr::rename(gene_caller_id = Accession) %>% 
  select_at(vars(!contains("Found"))) %>% 
  rename_with(~gsub("Abundance\\.|\\.Sample","",.), everything()) %>% 
  rename_at(all_of(exoP_sampleIDS), ~ exoP_sample$Sample_name) %>% 
  filter(Number.of.PSMs >=2 , Number.of.Unique.Peptides>=1)%>% 
  mutate_if(is.numeric, funs(replace_na(., 0))) %>% 
  mutate_if(is.numeric,as.integer)


########################################
#merge the dataset into a phyloseq object for convenient data management
########################################
#protein counts as otu_table
exoP_counts<- exoP_filt %>% select(c("gene_caller_id", exoP_sample$Sample_name))
exoP_counts<- otu_table(data.frame(exoP_counts[, all_of(exoP_sample$Sample_name)], row.names = exoP_counts$gene_caller_id), taxa_are_rows = TRUE)

#gene ids as taxonomy table
annotation<- tax_table(as.matrix(data.frame(genes_meta, row.names = genes_meta$gene_caller_id)))

#metadata
exoP_meta <- sample_data(exoP_sample)

#merge into phyloseq object
exoP_obj0<- phyloseq(exoP_counts, annotation, exoP_meta)

#remove proteins that were not observed 
exoP_obj0<- prune_taxa(taxa_sums(exoP_obj0)>0,exoP_obj0)

###################
#Plot number of proteins per sample
###################
# number of proteins per sample
prot_per_sample <- estimate_richness(prot_obj0, split = TRUE, measures = "Observed") %>% 
                    mutate(Sample_name = row.names(.)) %>% 
                      left_join(prot_sample, by = "Sample_name")

ggplot(prot_per_sample, aes(x = Sample_name, y = Observed, fill = Type)) + 
  #facet_wrap(Fraction~.) +
  geom_col()+
  ylab("# of proteins \n")+
  #scale_y_log10()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(size=14),legend.position = "bottom", 
        axis.title.x = element_blank(), axis.text.x = element_text(angle=90))

###################
#Barplots
###################
#conduct NSAF transformation
prot_nsaf<- add_nsaf(prot_obj0, "prot_length")

prot_nsaf.long <- psmelt(prot_nsaf)

#remove below 3% ra
taxa_classes <- unique(prot_nsaf.long$t_genus[!prot_nsaf.long$Abundance<0.005])

prot_nsaf.long$t_genus[prot_nsaf.long$Abundance<0.005] <- "Other taxa"

prot_nsaf.long$t_genus <- factor(prot_nsaf.long$t_genus,
                                 levels=c(sort(taxa_classes),"Other taxa"))

prot_nsaf.long_sub<- prot_nsaf.long %>%  filter(t_genus %in% c("Pseudoalteromonas", "Alteromonas", "Vibrio", "Synechococcus"))

ggplot(prot_nsaf.long_sub, 
       aes(x = Sample_name, y = Abundance,
           fill = t_genus)) + 
  #facet_grid(Type~Fraction, space= "fixed") +
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
prot_frac_agg.dds <- phyloseq_to_deseq2(prot_obj0, ~1)
prot_frac_agg.dds = estimateSizeFactors(prot_frac_agg.dds)
prot_frac_agg.dds <- estimateDispersions(prot_frac_agg.dds)
prot.vst <- getVarianceStabilizedData(prot_frac_agg.dds)

#make sure that the dimensions of the protein table and the DEseq object are matching
dim(prot.vst)
dim(otu_table(prot_obj0))


prot_frac_agg.vst<-prot_obj0
otu_table(prot_frac_agg.vst)<- otu_table(prot.vst, taxa_are_rows = TRUE)

#produce ordination
ps_obj_pca <- ordinate(prot_frac_agg.vst, method = "RDA", distance = "euclidean")
ps_obj_pca.df <- plot_ordination(prot_frac_agg.vst, ps_obj_pca, axes = c(1,2,3),justDF = TRUE) 

ps_obj_pca.df$Sample_name<- row.names(ps_obj_pca.df)

#extract explained variance
ps_obj_pca.evals <- 100 * summary(ps_obj_pca)$cont$importance[2, c("PC1","PC2")]

ggplot(data = ps_obj_pca.df, aes(x = PC1, y = PC2))+
  geom_point(fill = "black", size = 6,alpha = 0.8) +
  geom_point(aes(colour = Type), size = 4,alpha = 0.8) +
  geom_text(aes(x = PC1, y = PC2,label = Sample_name), 
            nudge_y= -1,size=4)+
  labs(x = sprintf("PC1 [%s%%]", round(ps_obj_pca.evals[1], 2)), 
       y = sprintf("PC2 [%s%%]", round(ps_obj_pca.evals[2], 2)))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(size=14),legend.position = "bottom")

#save the plot
ggsave(paste0(wd,"/R_figures/PCA_plot.png"), 
       plot = ordination_plot,
       units = "cm",
       width = 15, height = 15, 
       #scale = 1,
       dpi = 300)

#test significance of clustering
df <- as(sample_data(prot_frac_agg.vst), "data.frame")
d <- phyloseq::distance(prot_frac_agg.vst, "euclidean")
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

###################
#Identify enriched proteins in MetaP
###################
prot_frac_agg_metaP<- subset_samples(prot_obj0, Type != "T0")
prot_frac_agg_metaP.ddsMat <- phyloseq_to_deseq2(prot_frac_agg_metaP, ~Type)
prot_frac_agg_metaP.ddsMat = estimateSizeFactors(prot_frac_agg_metaP.ddsMat)
prot_frac_agg_metaP.ddsMat <- estimateDispersions(prot_frac_agg_metaP.ddsMat)
prot_frac_agg_metaP.DEseq <- DESeq(prot_frac_agg_metaP.ddsMat, fitType="local")
prot_frac_agg_metaP.DEseq.res <- results(prot_frac_agg_metaP.DEseq)

#extract only significant proteins
prot_frac_agg_metaP.DEseq.res.sig <- as(prot_frac_agg_metaP.DEseq.res, "data.frame") %>%  filter(padj < 0.1)

prot_frac_agg_metaP.DEseq.res.sig <- cbind(as(prot_frac_agg_metaP.DEseq.res.sig, "data.frame"),
                                           as(tax_table(prot_obj0)[rownames(prot_frac_agg_metaP.DEseq.res.sig), ], "matrix"))



###################
#Import metabolism estimates for the bins
###################
modules_info <- read.csv(paste(wd,"metaG_analysis/metaG_anvio/08_METABOLISM/modules_info.txt",sep=""), sep="\t", h= T)

Bins_kofam_hits <- read.csv(paste(wd,"metaG_analysis/metaG_anvio/08_METABOLISM/Bins_kofam_hits.txt",sep=""), sep="\t", h= T)

Bins_modules <- read.csv(paste(wd,"metaG_analysis/metaG_anvio/08_METABOLISM/Bins_modules.txt",sep=""), sep="\t", h= T)


#merge the enriched proteins with the bins
prot_frac_agg_metaP.DEseq.res.sig_Bin<- prot_frac_agg_metaP.DEseq.res.sig %>% plyr::rename(replace = c("gene_callers_id"="gene_caller_id")) %>%
  mutate(gene_caller_id = as.integer(gene_caller_id)) %>% 
  left_join(., Bins_kofam_hits, by = "gene_caller_id")

#see how many enriched proteins were associated with bins
prot_frac_agg_metaP.DEseq.res.sig_Bin %>% group_by(genome_name) %>% summarize(n())


#subset each bins proteins
prot_frac_agg_metaP.DEseq.res.sig_Bin_P_phenolica<-prot_frac_agg_metaP.DEseq.res.sig_Bin %>% filter(genome_name == "Bin_84_1") 

prot_frac_agg_metaP.DEseq.res.sig_Bin_Vibrio<-prot_frac_agg_metaP.DEseq.res.sig_Bin %>% filter(genome_name == "Bin_38_1") 

prot_frac_agg_metaP.DEseq.res.sig_Bin_Bermanella<-prot_frac_agg_metaP.DEseq.res.sig_Bin %>% filter(genome_name == "Bin_115_1") 

prot_frac_agg_metaP.DEseq.res.sig_Bin_Glaciecola<-prot_frac_agg_metaP.DEseq.res.sig_Bin %>% filter(genome_name == "Bin_115_2") 

prot_frac_agg_metaP.DEseq.res.sig_Bin_Alteromonas<-prot_frac_agg_metaP.DEseq.res.sig_Bin %>% filter(genome_name %in% c("Bin_2_1","Bin_2_2"))

#not assigned to any bin
prot_frac_agg_metaP.DEseq.res.sig_no_bin<-prot_frac_agg_metaP.DEseq.res.sig_Bin %>% filter(is.na(genome_name) == TRUE) 


#split posible annotations into separated rows
test <- prot_frac_agg_metaP.DEseq.res.sig_no_bin %>% separate_rows(GO_accession, sep = '\\|')

