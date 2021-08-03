# https://uclouvain-cbio.github.io/BSS2019/figs/cancer_9x9.html#data-exploration

require(dplyr)
require(tidyr)
require(stringr)
require(phyloseq)
require(ggplot2)
require(vegan)
require(microbiome)
require(pheatmap)
require(limma)
require(ggpubr)

#require(MSqRob)
#require(pathview)
#require(tidyverse)
#require(gplots)

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

#define function #https://github.com/hms-dbmi/UpSetR/issues/85
pres_abs_matrix <- function (input) {
  # Same as original fromList()...
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  # ... Except now it conserves your original value names!
  row.names(data) <- elements
  return(data)
}

#set working directory
#macOS
wd <- "/Users/eduardfadeev/Google Drive (dr.eduard.fadeev@gmail.com)/DECOMB/"

#Linux 
wd <- "~/Data/Postdoc-Vienna/DECOMB/"

#Windows
wd <- "D:/Postdoc-Vienna/DECOMB/"

########################################
#import taxonomy and annotation of each gene in the reference metagenome
########################################
#taxonomy table
tax_table <- read.csv(paste(wd,"metaG_analysis/metaG_anvio/05_ANVIO/spades_merged/spades-tax-names.txt",sep=""),
                     sep="\t", h= T)
#taxonomic classification
gene_tax_table <- read.csv(paste(wd,"metaG_analysis/metaG_anvio/05_ANVIO/spades_merged/spades-genes-taxonomy.txt",sep=""),
                           sep="\t", h= T)

#gene annotation list by sources
gene_annotations<- list()
gene_annotations[["gene_ids"]] <- read.csv(paste(wd,"metaG_analysis/metaG_anvio/05_ANVIO/spades_merged/spades-gene-calls.txt",sep=""),
                           sep="\t", h= T)
                           
for (src in c("COG20_FUNCTION","KeggGhostKoala","GO","Pfam","InterPro", "Hamap")){
gene_annotations[[src]] <- read.csv(paste(wd,"metaG_analysis/metaG_anvio/05_ANVIO/spades_merged/spades-",src,"-functions.txt",sep=""),
                           sep="\t", h= T)
}


#merge taxonomy and gene ids together and calculate protein length
genes_meta <- merge(gene_tax_table,tax_table, by ="taxon_id", all.x = TRUE) %>% 
                merge(gene_annotations[["gene_ids"]], by ="gene_callers_id", all = TRUE) %>% 
                  select(c("gene_callers_id", "contig", "start", "stop","aa_sequence",
                           "taxon_id", "t_genus", "t_species")) %>% 
                    group_by(aa_sequence) %>% 
                    mutate(prot_length = nchar(aa_sequence))
                
########################################
#import metaproteome data
########################################
#import sample list of metaproteomes
prot_sample <- read.csv(paste(wd,"metaP_analysis/quant-no-groups/DECOMB-consensus-no-groups_InputFiles.txt", sep=""),sep="\t", h= T) %>% 
  mutate(Run = case_when(grepl("20200622_TinkaraTinta2",File.Name) == TRUE ~ 2, TRUE ~ 1),
          File.Name =   tools::file_path_sans_ext(gsub("Z:\\\\EFadeev\\\\DECOMB_raw_data\\\\20200622_TinkaraTinta2\\\\|Z:\\\\EFadeev\\\\DECOMB_raw_data\\\\20200615_TinkaraTinta\\\\",
                                                      "", File.Name)),
          Fraction = case_when(grepl("MP",File.Name) == TRUE ~ "metaP", TRUE ~ "ExoP"),
           Type = case_when(grepl("C",File.Name) ==TRUE ~"Control", 
                            grepl("J",File.Name) ==TRUE ~"Jelly",
                            grepl("T0",File.Name) ==TRUE ~ "T0"),
         Replicate = gsub("_.*","",File.Name),
         Sample.ID = paste(File.ID, File.Name, Type, Run, sep="_"),
         Group = as.factor(paste(Replicate,Fraction, Run, sep="_"))) # %>% 
  # separate(File.Name, into = c("Sample","Fraction"), sep ="_", remove = FALSE)

#import proteins data
protein_raw <- read.csv(paste(wd,"metaP_analysis/quant-no-groups/DECOMB-consensus-no-groups_Proteins.txt", sep=""),sep="\t", h= T)

#rename sample columns and filter out low confidence proteins and replace NA with 0
protein_filt <- protein_raw %>% 
  dplyr::rename(gene_callers_id = Accession) %>% 
  select_at(vars(!contains("Found"))) %>% 
  rename_with(~gsub("Abundance\\.|\\.Sample","",.), everything()) %>% 
  rename_at(vars(prot_sample$File.ID), ~ prot_sample$Sample.ID) %>% 
  filter(Number.of.PSMs >=2 , Number.of.Unique.Peptides>=1)%>% 
  mutate_if(is.numeric, funs(replace_na(., 0))) %>% 
  mutate_if(is.numeric,as.integer)

 # left_join(gene_annotations[["gene_ids"]], by ="gene_callers_id", all = TRUE) %>% 

########################################
#merge the dataset into a phyloseq object for convenient data management
########################################
#protein counts as otu_table
protein_counts<- protein_filt %>% select(c("gene_callers_id", prot_sample$Sample.ID))
protein_counts<- otu_table(data.frame(protein_counts[, prot_sample$Sample.ID], row.names = protein_counts$gene_callers_id), taxa_are_rows = TRUE)

#gene ids as taxonomy table
annotation<- tax_table(as.matrix(data.frame(genes_meta, row.names = genes_meta$gene_callers_id)))

#metadata
sample_meta <- sample_data(data.frame(prot_sample, row.names =prot_sample$Sample.ID))

#merge into phyloseq object
prot_obj0<- phyloseq(protein_counts, annotation, sample_meta)

#remove proteins that were not observed 
prot_obj0<- prune_taxa(taxa_sums(prot_obj0)>0,prot_obj0)

###################
#Plot number of proteins per sample
###################
# number of proteins per sample
prot_per_sample <- estimate_richness(prot_obj0, split = TRUE, measures = "Observed") %>% 
                   mutate(File.Name = row.names(.)) %>% 
                     mutate(Fraction = case_when(grepl("MP",File.Name) == TRUE ~ "metaP", TRUE ~ "ExoP"),
                            Type = case_when(grepl("C",File.Name) ==TRUE ~"Control", 
                              grepl("J",File.Name) ==TRUE ~"Jelly",
                              grepl("T0",File.Name) ==TRUE ~ "T0")) %>% 
                      mutate(Type = factor(Type, levels = c("T0","Control","Jelly")),
                             Sample.ID = gsub("F._|F.._","",File.Name),
                             Tech.Rep = gsub(".*_","",File.Name))

for (frac in c("metaP","ExoP")){
  sub<- prot_per_sample %>% filter(Fraction ==frac)
prot_counts_bar.p[[frac]] <- ggplot(sub, aes(x = Sample.ID, y = Observed, fill = Type)) + 
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

ggarrange(prot_counts_bar.p[["metaP"]], prot_counts_bar.p[["ExoP"]],
          ncol = 1, nrow = 2, align = "v")


#test differences between runs
prot_per_sample_test <- prot_per_sample   %>%
  group_by(Fraction) %>% 
  rstatix::t_test(Observed ~ Replicate, paired = TRUE, p.adjust.method = "BH") %>%
  add_significance()

###################
#Combine EL and EH exoP fractions
###################
#aggregate the fractions
prot_frac_agg <- merge_samples(prot_obj0, "Group", fun = sum)
otu_table(prot_frac_agg)<- t(otu_table(prot_frac_agg))
#generate aggregated metadata
meta <- sample_data(prot_obj0)
meta <- as(meta, "data.frame") %>% select(Group, Fraction, Type, Run) %>% 
  unique() %>% mutate(Fraction = as.factor(Fraction),
                      Type = as.factor(Type))
row.names(meta)<- meta$Group

sample_data(prot_frac_agg)<- sample_data(meta)

###################
#Upset protein overlaps between replicates
###################
y<- list()

for (i in sample_names(prot_frac_agg)){
  sub_sample <- prune_samples(i, prot_frac_agg)
  sub_sample <- prune_taxa(taxa_sums(sub_sample)>0,sub_sample)
  y[[i]] <- as.character(row.names(otu_table(sub_sample)))
}

#generate overlap matrix
protein_overlaps <- pres_abs_matrix(y)    
protein_overlaps$gene_callers_id <- rownames(protein_overlaps)

#plot
upset(protein_overlaps, number.angles = 30,
      sets = as.vector(sample_names(prot_frac_agg)),
      keep.order = TRUE, 
      mainbar.y.label = "No. of overlaping proteins",
      order.by = "freq" ) #,empty.intersections = "on")

###################
#Generate NMDS plot
###################
#conduct NSAF transformation
prot_nsaf<- add_nsaf(prot_frac_agg, "prot_length")

# Do CLR transformation to protein abundances
prot_nsaf_clr <- microbiome::transform(prot_nsaf, "clr")

ps_obj_nmds <- ordinate(prot_nsaf_clr, method = "RDA", distance = "euclidean")
ps_obj_nmds.df <- plot_ordination(prot_nsaf_clr, ps_obj_nmds, axes = c(1,2,3),justDF = TRUE)

ps_obj_nmds.df$Sample<- row.names(ps_obj_nmds.df)

#extract explained variance
ps_obj_nmds.evals <- 100 * summary(ps_obj_nmds)$cont$importance[2, c("PC1","PC2")]

ordination_plot<- ggplot(data = ps_obj_nmds.df, aes(x = PC1, y = PC2, shape = Fraction))+
  geom_point(fill = "black", size = 6,alpha = 0.8) +
  geom_point(aes(colour = Type), size = 4,alpha = 0.8) +
  geom_text(aes(x = PC1, y = PC2,label = Sample), 
            nudge_y= -1,size=5)+
  labs(x = sprintf("PC1 [%s%%]", round(ps_obj_nmds.evals[1], 2)), 
       y = sprintf("PC2 [%s%%]", round(ps_obj_nmds.evals[2], 2)))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(size=14),legend.position = "bottom")

#save the plot
ggsave(paste0(wd,"/R_figures/PCoA_plot.pdf"), 
       plot = NMDS_plot,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)

#test significance of clustering
df <- as(sample_data(prot_nsaf_clr), "data.frame")
d <- phyloseq::distance(prot_nsaf_clr, "euclidean")
adonis_all <- adonis2(d ~ Fraction + Type + Run  , df)
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
ps_obj_nsaf_metaP<- subset_samples(prot_nsaf_clr, Type != "T0" & Fraction =="metaP")
ps_obj_nsaf_metaP_abund<- abundances(ps_obj_nsaf_metaP)

# Prepare the design matrix which states the groups for each sample
design <- cbind(intercept = 1, Grp2vs1 = sample_data(ps_obj_nsaf_metaP)[["Type"]])
rownames(design) <- row.names(sample_data(ps_obj_nsaf_metaP))
design <- as.data.frame(design[sample_names(ps_obj_nsaf_metaP), ])

# Fit the limma model
fit <- lmFit(ps_obj_nsaf_metaP_abund, design, adjust="BH")
fit <- eBayes(fit)

#explore results
results <- decideTests(fit[,"Grp2vs1"])
summary(results)

# NOTE: results and p-values are given for all groupings in the design matrix
coef.index <- 2

# Limma P-values
pvalues.limma = fit$p.value[, 2]

# Limma effect sizes
efs.limma <-  fit$coefficients[, "Grp2vs1"]

# QQ plot
qqt(fit$t[, coef.index], df = fit$df.residual + fit$df.prior); abline(0,1)

# Volcano
volcanoplot(fit, coef = coef.index, highlight = coef.index)

# Summarise
knitr::kable(topTable(fit, coef = coef.index, p.value=0.1), digits = 2)

#plot heat map of the first 100 enriched proteins
tab_100 <- topTable(fit, n=100, sort.by = "logFC", coef=coef.index) %>% filter(adj.P.Val<0.1)
prot_abund_100<- ps_obj_nsaf_metaP_abund[row.names(tab_100),]




pheatmap(prot_abund_100)










###################
#Plot proportions of proteins per sample
###################


prot_NSAF.long <- psmelt(ps_obj_nsaf)

prot_NSAF.agg<- prot_NSAF.long %>% 
  group_by(Sample) %>% 
  summarize(Abundance, "sum")


prot_NSAF_bar.p <- ggplot(prot_NSAF.long, aes(x = Sample, y = Abundance)) + 
  facet_grid(Fraction~., space= "fixed") +
  geom_col()+
  ylab("# of proteins \n")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(size=14),legend.position = "bottom", 
        axis.title.x = element_blank(), axis.text.x = element_text(angle=90))






#merge runs of the analysis
ps_obj_agg <- merge_samples(ps_obj0, "Group", fun = sum)

#generate aggregated metadata
meta <- sample_data(ps_obj0)
meta <- as(meta, "data.frame") %>% select(Group, Fraction, Type) %>% 
   unique() %>% mutate(Fraction = as.factor(Fraction),
                                               Type = as.factor(Type))
row.names(meta)<- meta$Group

sample_data(ps_obj_agg)<- sample_data(meta)


















###################
#plot a heatmap
###################
metaP_agg<- subset_samples(ps_obj_nsaf, Fraction =="metaP")

# Do log10 transformation to protein abundances
ps_clr_abund <- microbiome::abundances(ps_obj_nsaf)

#select the most abundant proteins
top_prot<- sort(rowSums(ps_clr_abund), decreasing = T)
ps_clr_abund_top<- ps_clr_abund[names(top_prot[1:50]),]

#plot heatmap
pheatmap(ps_clr_abund_top)
















#plot heatmap
pheatmap(otu_table(ps_obj_nsaf)[1:100,])



#metaP control vs. jelly
ps_obj0_no_T0<- subset_samples(ps_obj0, Type != "T0" & Fraction == "metaP")
ps_obj0_no_T0<- prune_taxa(taxa_sums(ps_obj0_no_T0)>0,ps_obj0_no_T0)


ps_obj_agg <- merge_samples(ps_obj0_no_T0, "Group", fun = mean)
#generate metadata for limma
meta <- sample_data(ps_obj0_no_T0)

meta <- as(meta, "data.frame") %>% select(Group, Fraction, Type) %>% 
  filter(Type != "T0") %>% unique() %>% mutate(Fraction = as.factor(Fraction),
                                               Type = as.factor(Type))
row.names(meta)<- meta$Group

sample_data(ps_obj_agg)<- sample_data(meta)


# Get protein abundances and sample metadata
prot_abund_log10 <- microbiome::abundances(microbiome::transform(ps_obj_agg, "clr"))

ps_clr_sub<- subset_samples(ps_clr, Type != "T0" & Fraction =="metaP")

ps_clr_sub <- microbiome::abundances(ps_clr_sub)

# Prepare the design matrix which states the groups for each sample
design <- cbind(intercept = 1, Grp2vs1 = meta[["Type"]])
rownames(design) <- row.names(meta)
design <- as.data.frame(design[colnames(ps_clr_sub), ])

# NOTE: results and p-values are given for all groupings in the design matrix
coef.index <- 1

# Fit the limma model
fit <- lmFit(ps_clr_sub, design, adjust="BH")
fit <- eBayes(fit)

# Limma P-values
pvalues.limma = fit$p.value[, 2]

# Limma effect sizes
efs.limma <-  fit$coefficients[, "Grp2vs1"]

# Summarise
knitr::kable(topTable(fit, coef = coef.index, p.value=0.1), digits = 2)

# QQ plot
qqt(fit$t[, coef.index], df = fit$df.residual + fit$df.prior); abline(0,1)

# Volcano
volcanoplot(fit, coef = coef.index, highlight = coef.index)


#plot heat map of the first 100 enriched proteins
tab_100 <- topTable(fit, n=Inf, sort.by = "logFC", coef=coef.index) %>% filter(adj.P.Val<0.1)
prot_abund_100<- ps_clr_sub[row.names(tab_100),]

heatmap.2(prot_abund_100, col= bluered(10), scale = "none", margins = c(10, 5))


#map KOs of enriched proteins on KEGG maps
enriched_prot <- topTable(fit, n=Inf, coef=coef.index) %>% filter(adj.P.Val<0.1)
enriched_prot$gene_callers_id<- as.integer(row.names(enriched_prot))

#merge protein data with annotation and exlude those without KO
enriched_prot <- left_join(enriched_prot, genes_meta, by ="gene_callers_id") %>% filter(KEGG !="")


pv.out <- pathview(gene.data = enriched_prot$KEGG, 
                   pathway.id ="00190 ",
                   species = "ko", 
                   keys.align = "y", 
                   kegg.native = T, both.dirs = TRUE, 
                   low = "blue", mid = "gray", high = "red", bin = 20)





#Exop control vs. jelly
ExoP_no_T0<- subset_samples(ps_obj0, Type != "T0" & Fraction == "ExoP")
ExoP_no_T0<- prune_taxa(taxa_sums(ExoP_no_T0)>0,ExoP_no_T0)

ExoP_agg <- merge_samples(ExoP_no_T0, "Group", fun = mean)

#generate metadata for limma
meta <- sample_data(ExoP_agg)

meta <- as(meta, "data.frame") %>% select(Group, Fraction, Type) %>% 
  filter(Type != "T0") %>% unique() %>% mutate(Fraction = as.factor(Fraction),
                                               Type = as.factor(Type))
row.names(meta)<- meta$Group

sample_data(ExoP_agg)<- sample_data(meta)


# Get protein abundances and sample metadata
ExoP_agg_log10 <- microbiome::abundances(microbiome::transform(ExoP_agg, "clr"))

# Prepare the design matrix which states the groups for each sample
design <- cbind(intercept = 1, Grp2vs1 = meta[["Type"]])
rownames(design) <- row.names(meta)
design <- as.data.frame(design[colnames(ExoP_agg_log10), ])

# NOTE: results and p-values are given for all groupings in the design matrix
coef.index <- 1

# Fit the limma model
ExoP_fit <- lmFit(ExoP_agg_log10, design, adjust="BH")
ExoP_fit <- eBayes(ExoP_fit)

# Limma P-values
pvalues.limma = ExoP_fit$p.value[, 2]

# Limma effect sizes
efs.limma <-  ExoP_fit$coefficients[, "Grp2vs1"]

# Summarise
knitr::kable(topTable(ExoP_fit, coef = coef.index, p.value=0.1), digits = 2)

# QQ plot
qqt(ExoP_fit$t[, coef.index], df = ExoP_fit$df.residual + ExoP_fit$df.prior); abline(0,1)

# Volcano
volcanoplot(ExoP_fit, coef = coef.index, highlight = coef.index)


#plot heat map of the first 100 enriched proteins
ExoP_top_100 <- topTable(ExoP_fit, n=100, sort.by = "logFC", coef=coef.index) %>% filter(adj.P.Val<0.1)
ExoP_top_100<- ExoP_agg_log10[row.names(ExoP_top_100),]

gplots::heatmap.2(ExoP_top_100, col= bluered(10), scale = "none", margins = c(10, 5))


#map KOs of enriched proteins on KEGG maps
ExoP_enriched <- topTable(ExoP_fit, n=Inf, coef=coef.index) %>% filter(adj.P.Val<0.1)
ExoP_enriched$gene_callers_id<- as.integer(row.names(ExoP_enriched))

#merge protein data with annotation and exlude those without KO
ExoP_enriched <- left_join(ExoP_enriched, genes_meta, by ="gene_callers_id") %>% filter(KEGG !="")


pv.out <- pathview(gene.data = ExoP_enriched$KEGG, 
                   pathway.id ="00190 ",
                   species = "ko", 
                   keys.align = "y", 
                   kegg.native = T, both.dirs = TRUE, 
                   low = "blue", mid = "gray", high = "red", bin = 20)

















ggplot(test_sub, aes(x= KEGG, y = logFC, colour = t_genus))+
  geom_point(size = 5)+
  coord_flip()+
  geom_vline(xintercept=0)+
  theme(legend.position = "bottom")


                   








topTable(fit, n=100, coef=coef.index) %>% filter(adj.P.Val<0.1 & abs(logFC)>2)








results <- decideTests(fit)

vennDiagram(results)




                   


otu_table(ps_obj_agg)<- otu_table(prot_abund_log10, taxa_are_rows = T)



#NMDS plot
ps_obj_nsaf.ord <- ordinate(ps_obj_agg, method = "MDS", distance = "euclidean")
ps_obj_nsaf.ord.df <- plot_ordination(ps_obj_agg, ps_obj_nsaf.ord, axes = c(1,2,3),justDF = TRUE)

ps_obj_nsaf.ord.df$Sample<- row.names(ps_obj_nsaf.ord.df)

ggplot(data = ps_obj_nsaf.ord.df)+
  geom_point(aes(x = Axis.1, y = Axis.2), 
             fill = "black", size = 5,alpha = 0.8) +
  geom_text(aes(x = Axis.1, y = Axis.2,label = Sample), 
            nudge_y= -0.01,size=5)


+
  #geom_point(data = Dor.ord.df, aes(x = NMDS1, y = NMDS2, colour = Mic.Season, shape = Year), 
  #           size = 7) +
  geom_point(aes(x = NMDS1, y = NMDS2, colour = Temp_degC, shape = location), 
             size = 3,alpha = 0.8) +
  geom_text(aes(x = NMDS1, y = NMDS2,label = Comment), 
            nudge_y= -8,size=5)+
  scale_colour_gradient(low = "blue", high = "yellow")+
  annotate(geom="text", x=-100, y=100, label= paste0("Stress = ", round(Dor_ps.gm_mean.ord$stress,2)),
           color="red", size = 5)+
  coord_fixed()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(size=14),legend.position = "bottom")





#sum the two runs together

# NAAF transformation
protein_trans <- protein_filt %>% dplyr::mutate_at(c(prot_sample$Sample.ID), funs(replace_na(., 0))) %>% 
  dplyr::mutate_at(c(prot_sample$Sample.ID), funs(./Number.of.AAs)) %>% 
  dplyr::mutate_at(c(prot_sample$Sample.ID), funs(NAAF = ./ sum(.)))

#NMDS
prot.pca <- vegan::metaMDS(t(protein_trans[74:115]))


sample_nmds <- as.data.frame(prot.pca$points) %>% tibble::rownames_to_column() %>% 
                   separate(rowname, into = c("ID","Sample","Fraction","Type","Replicate", "Method"), sep ="_", remove = TRUE) %>% 
                    mutate(Type = case_when(Sample %in% c("C1","C2","C3") ~"Control", 
                                            Sample %in% c("J1","J2","J3") ~"Jelly",
                                                      Sample == "T0" ~ "T0"))

protein_nmds <- as.data.frame(prot.pca$species) %>% tibble::rownames_to_column() %>% 
  separate(rowname, into = c("ID","Sample","Fraction","Type","Replicate", "Method"), sep ="_", remove = TRUE) %>% 
  mutate(Type = case_when(Sample %in% c("C1","C2","C3") ~"Control", 
                          Sample %in% c("J1","J2","J3") ~"Jelly",
                          Sample == "T0" ~ "T0"))

ggplot()+
  geom_point(data=protein_nmds, aes(x=MDS1,y=MDS2), size = 1, alpha = 0.3, colour = "gray50")+
  geom_point(data=sample_nmds, aes(x=MDS1,y=MDS2, colour = Type, shape = as.factor(Type), label = Replicate),
             size = 5)+
  geom_text(data=sample_nmds, aes(x=MDS1,y=MDS2, colour = Type, shape = as.factor(Type), label = Replicate),
            size = 5, nudge_y = -0.2)+
  theme_bw()


#enrichment test
require(DESeq2)


protein_no_na <- protein_filt %>% mutate_if(is.numeric, funs(replace_na(., 0))) %>% 
  mutate_if(is.numeric,as.integer)

dds <- DESeqDataSetFromMatrix(countData=protein_no_na[,c(4,20:61)], 
                              colData=prot_sample, 
                              design=~Type, tidy = TRUE)

prot.DEseq <- DESeq(dds)
prot.DEseq.res <- results(prot.DEseq)



require("msmsTests")

samples_frac <- prot_sample %>% filter(Type != "T0", Fraction == "metaP") %>% select(c(9:11))
row.names(samples_frac)<- samples_frac$Sample.ID

protein_frac <- protein_no_na %>% select(c("gene_callers_id", samples_frac$Sample.ID))

prot_MSnSet <- readMSnSet2(protein_frac, ecol = samples_frac$Sample.ID, fnames = 1)

phenoData(prot_MSnSet) <- AnnotatedDataFrame(samples_frac)
                    

e <- pp.msms.data(prot_MSnSet)

null.f <- "y~1"
alt.f <- "y~Type"
div <- apply(exprs(e),2,sum)
edgeR_res <- msms.edgeR(e, form0= null.f, form1= alt.f, div=div, fnm="Type")


test<- test.results(edgeR_res, e, gpf = pData(e)$Type, gp1="Control",gp2="Jelly",
                    alpha = 0.05, minLFC=1, div= div,
                    method="BH")

res.volcanoplot(test$tres,max.pval=0.1, maxy=2, maxx = 30)



results<- data.frame(test$tres, gene_callers_id = as.integer(row.names(test$tres))) %>% 
            left_join(genes_meta, by = "gene_callers_id") %>% 
  filter(DEP =="TRUE")




heatmap(as.matrix(protein_frac[2:13]))


MAplot(qnt[, c(4, 2)], cex = .9, col = cls, pch = pch, show.statistics = FALSE)






#add reference sequence and replace variants with ASVs

prot_AAs<- as.vector(genes_meta$aa_sequence)

prot_AAs_length <- length(genes_meta$aa_sequence)

ps2 <- Biostrings::AAStringSet(prot_AAs)
names(ps2) <- taxa_names(ps)




