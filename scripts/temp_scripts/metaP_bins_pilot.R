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
require(rstatix)

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
gene_ids <- read.csv(paste(wd,"metaG_analysis/metaG_anvio/05_ANVIO/spades_merged/spades-gene-calls.txt",sep=""),
                     sep="\t", h= T)
gene_annotations_df<- gene_ids

for (src in c("COG20_FUNCTION","KeggGhostKoala","GO","Pfam","InterPro", "Hamap")){
  gene_annotations <- read.csv(paste(wd,"metaG_analysis/metaG_anvio/05_ANVIO/spades_merged/spades-",src,"-functions.txt",sep=""),
                               sep="\t", h= T)%>% select(gene_callers_id, accession, function.) %>% 
    group_by(gene_callers_id)%>%
    summarise_each(funs(paste(unique(.), collapse='|')),matches('^\\D+$')) %>% 
    plyr::rename(replace= c("accession"=paste(src, "accession", sep ="_"), "function."=paste(src, "function", sep ="_")))
  
  gene_annotations_df<- merge(gene_annotations_df,gene_annotations, by ="gene_callers_id", all = TRUE )
}

#merge taxonomy and gene ids together and calculate protein length
genes_meta <- merge(gene_tax_table,tax_table, by ="taxon_id", all.x = TRUE) %>% 
  merge(gene_annotations_df, by ="gene_callers_id", all = TRUE) %>% 
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
         Group = as.factor(paste(Replicate,Fraction, Run, sep="_")))  %>% 
  mutate(Sample.ID = gsub("_200616145603","", Sample.ID),
         File.Name = gsub("_200616145603","", File.Name)) #correct sample name
# separate(File.Name, into = c("Sample","Fraction"), sep ="_", remove = FALSE)

#import proteins data
protein_raw <- read.csv(paste(wd,"metaP_analysis/Bins_prot/Bin_115_1/Bin-115_1-consensus_Proteins.txt", sep=""),sep="\t", h= T)

#rename sample columns and filter out low confidence proteins and replace NA with 0
protein_filt <- protein_raw %>% 
  dplyr::rename(gene_callers_id = Accession) %>% 
  select_at(vars(!contains("Found"))) %>% 
  rename_with(~gsub("Abundance\\.|\\.Sample","",.), everything()) %>% 
  rename_at(vars(prot_sample$File.ID), ~ prot_sample$Sample.ID) %>% 
  filter(Number.of.PSMs >=2 , Number.of.Unique.Peptides>=1)%>% 
  mutate_if(is.numeric, funs(replace_na(., 0))) %>% 
  mutate_if(is.numeric,as.integer)



protein_filt_ann<- protein_filt %>% left_join(genes_meta, by ="gene_callers_id")

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
         Sample.ID = gsub("Control_|Jelly_|T0_","", gsub("F._|F.._","",File.Name)),
         Tech.Rep = gsub(".*_","",File.Name))

prot_counts_bar.p<- list()
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
          ncol = 2, nrow = 1, align = "hv")

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
                      Type = as.factor(Type),
                      Merge = gsub("_1|_2","",Group))
row.names(meta)<- meta$Group

sample_data(prot_frac_agg)<- sample_data(meta)

###################
#Generate PCA plot
###################
#conduct NSAF transformation
prot_nsaf<- add_nsaf(prot_frac_agg, "prot_length")

# Do CLR transformation to protein abundances
prot_nsaf_clr <- microbiome::transform(prot_nsaf, "clr")

ps_obj_nmds <- ordinate(prot_nsaf_clr, method = "RDA", distance = "euclidean")
ps_obj_nmds.df <- plot_ordination(prot_nsaf_clr, ps_obj_nmds, axes = c(1,2,3),justDF = TRUE) 

ps_obj_nmds.df$Sample<- gsub("metaP_|ExoP_","",row.names(ps_obj_nmds.df))

#extract explained variance
ps_obj_nmds.evals <- 100 * summary(ps_obj_nmds)$cont$importance[2, c("PC1","PC2")]

ordination_plot<- ggplot(data = ps_obj_nmds.df, aes(x = PC1, y = PC2, shape = Fraction))+
  geom_point(fill = "black", size = 6,alpha = 0.8) +
  geom_point(aes(colour = Type), size = 4,alpha = 0.8) +
  geom_text(aes(x = PC1, y = PC2,label = Sample), 
            nudge_y= -0.2,size=4)+
  labs(x = sprintf("PC1 [%s%%]", round(ps_obj_nmds.evals[1], 2)), 
       y = sprintf("PC2 [%s%%]", round(ps_obj_nmds.evals[2], 2)))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(size=14),legend.position = "bottom")

ordination_plot

###################
#Identify enriched proteins in MetaP
###################
ps_obj_nsaf_metaP<- subset_samples(prot_nsaf_clr, Type != "T0" & Fraction =="metaP" & Run =="1")
ps_obj_nsaf_metaP_abund<- abundances(ps_obj_nsaf_metaP)

# Prepare the design matrix which states the groups for each sample
design <- cbind(intercept = 1, Grp2vs1 = as.factor(sample_data(ps_obj_nsaf_metaP)[["Type"]]))
rownames(design) <- row.names(sample_data(ps_obj_nsaf_metaP))
design <- as.data.frame(design[sample_names(ps_obj_nsaf_metaP), ])

# Fit the limma model
fit <- lmFit(ps_obj_nsaf_metaP_abund, design, adjust="fdr")
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
tab_100 <- topTable(fit, n=200, sort.by = "logFC", coef=coef.index) %>% filter(adj.P.Val<0.1)
prot_abund_100<- ps_obj_nsaf_metaP_abund[row.names(tab_100),]
pheatmap(prot_abund_100)


protein_no_na <- protein_filt %>% mutate_if(is.numeric, funs(replace_na(., 0))) %>% 
  mutate_if(is.numeric,as.integer)

prot_sample_sub<- prot_sample %>% filter(Group %in% colnames(ps_obj_nsaf_metaP_abund))

row.names(prot_sample_sub)<- prot_sample_sub$Group

dds <- DESeqDataSetFromMatrix(countData=ps_obj_nsaf_metaP_abund, 
                              colData= prot_sample_sub, 
                              design=~Type, tidy = TRUE)

prot.DEseq <- DESeq(dds)
prot.DEseq.res <- results(prot.DEseq)

