
########################################
#Metaproteome data analysis
########################################

#set working directory
#macOS
wd <- "/Users/eduardfadeev/Google Drive (dr.eduard.fadeev@gmail.com)/DECOMB/"

#Linux 
wd <- "~/Data/Postdoc-Vienna/DECOMB/"

#Windows
wd <- "F:/My Drive/DECOMB/"

#load libraries
require(dplyr)
require(tidyr)
require(phyloseq)
require(ggplot2)
require(ggpubr)
require(DESeq2)
require(vegan)
require(pheatmap)

#load metaproteome phyloseq object
metaP_obj0<- readRDS("data/metaP_ps_raw.rds")

###################
#Generate PCA plot 
###################
#conduct variance stabiliozation of the metaP dataset
metaP_obj0.ddsMat <- phyloseq_to_deseq2(metaP_obj0, ~Type)
metaP.DEseq.vsd <- vst(metaP_obj0.ddsMat,  fitType="local", blind=FALSE)

#plot PCA
metaP_pca.df <- plotPCA(metaP.DEseq.vsd, intgroup=c("Sample_name","Replicate","Type"),  returnData=TRUE)

metaP_pca.df<-  metaP_pca.df %>% mutate(SampleID = gsub("_MP","",Sample_name))

#extract explained variance
percentVar <- round(100 * attr(metaP_pca.df, "percentVar"))

#plot
metaP_ordination_plot<- ggplot(data = metaP_pca.df, aes(x = PC1, y = PC2, shape = Type))+
  geom_point(fill = "black", size = 6) +
  geom_point(aes(colour = Type), size = 4,alpha = 0.8) +
  geom_text(aes(x = PC1, y = PC2,label = SampleID), 
            nudge_y= -8,size=4)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_colour_manual(values =c("T0"="gray50",
                                "Jelly"="red",
                                "Control"="blue"))+
  coord_fixed()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(size=14),legend.position = "bottom")

#save the plot
ggsave(paste0(wd,"/R_figures/metaP_ordination.pdf"), 
       plot = metaP_ordination_plot,
       units = "cm",
       width = 15, height = 15, 
       #scale = 1,
       dpi = 300)

#test significance of clustering
df <- colData(metaP.DEseq.vsd)[,c("Sample_name","Replicate","Type")]
sample_distance <- dist(t(assay(metaP.DEseq.vsd)))
adonis_all <- adonis2(sample_distance ~ Type, df)
adonis_all

#posthoc to check which ponds are different
groups <- df[["Type"]]
mod <- betadisper(sample_distance, groups)
permutest(mod)

#dispersion is different between groups
plot(mod)
boxplot(mod)
mod.HSD <- TukeyHSD(mod)
mod.HSD
plot(mod.HSD)

###################
#Identify sig. enriched proteins between Jelly and Control bottles
###################
#run DESeq2 enrichment test
metaP.DEseq <- DESeq(metaP_obj0.ddsMat, fitType="local")

###################
#Generate heatmap of the top 100 proteins
###################
#select only the top 100 proteins and produce dataframe with metadata
select_prot <- order(rowMeans(counts(metaP.DEseq, normalized=FALSE)),
                     decreasing=TRUE)[1:100]
meta_df <- as.data.frame(colData(metaP.DEseq)[,c("Replicate","Type")])

#produce and save a heatmap
pheatmap(assay(metaP.DEseq.vsd)[select_prot,], cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=meta_df,
         main = "Top 100 proteins", filename = paste0(wd, "R_figures/metaP_heatmap.pdf"))


#export results for Jelly vs. Control
metaP.DEseq.res <- results(metaP.DEseq, name = "Type_Jelly_vs_Control")

#export results for Jelly vs. Control with shrinked LFC estimates
#(explained: http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#altshrink)
metaP.DEseq.res_LFC <- lfcShrink(metaP.DEseq, coef="Type_Jelly_vs_Control", type="apeglm")

#extract only significant proteins and merge with annotations
metaP.DEseq.res.sig <- as(metaP.DEseq.res_LFC, "data.frame") %>%  mutate(gene_caller_id = row.names(.)) %>% 
                            filter(padj < 0.1 ) %>% 
                            left_join(as.data.frame(tax_table(metaP_obj0)), by = "gene_caller_id")


###################
#summarize the significantly enr. proteins by Taxonomic orders
###################
metaP.DEseq.res.sig.TaxOrder <- metaP.DEseq.res.sig %>% 
  mutate(Type = case_when(log2FoldChange>1 ~ "Jelly",
                          log2FoldChange< -1 ~ "Control")) %>% 
  group_by(Type, Class, Order) %>% 
  filter(Type %in% c("Jelly","Control")) %>% 
  summarize(Total_prot.= length(Type))


#Alteromonadales
metaP.DEseq.res.sig.Alt<- metaP.DEseq.res.sig %>% 
  mutate(Type = case_when(log2FoldChange>1 ~ "Jelly",
                          log2FoldChange< -1 ~ "Control")
         ) %>% 
  filter(Type %in% c("Jelly","Control"),
         Order =="Alteromonadales")

#Vibrio
metaP.DEseq.res.sig.Vib<- metaP.DEseq.res.sig %>% 
  mutate(Type = case_when(log2FoldChange>1 ~ "Jelly",
                          log2FoldChange< -1 ~ "Control")
  ) %>% 
  filter(Type %in% c("Jelly","Control"),
         Order =="Vibrionales")

#Rhodobacters
metaP.DEseq.res.sig.Rhod<- metaP.DEseq.res.sig %>% 
  mutate(Type = case_when(log2FoldChange>1 ~ "Jelly",
                          log2FoldChange< -1 ~ "Control")
  ) %>% 
  filter(Type %in% c("Jelly","Control"),
         Order =="Rhodobacterales")


###################
#summarize the significantly enr. proteins by Protein families
###################
metaP.DEseq.res.sig.Pfam <- metaP.DEseq.res.sig %>% 
    mutate(Type = case_when(log2FoldChange>1 ~ "Jelly",
                          log2FoldChange< -1 ~ "Control")) %>% 
    group_by(Type, Pfam_accession, Pfam_function) %>% 
    filter(Type %in% c("Jelly","Control")) %>% 
    summarize(Total_prot.= length(Type))


###################
#Explore results by KEGG modules
###################
#combine the results with KEGG module hits
metaG_kofam_hits <- read.csv(paste(wd,"metaG_analysis/metaG_anvio/05_ANVIO/spades-Kofam_hits.txt",sep=""), sep="\t", h= T) %>% 
                      mutate(gene_caller_id= as.character(gene_caller_id))

metaP.DEseq.res_KOfam <- metaP.DEseq.res.sig %>% 
  left_join(metaG_kofam_hits, by ="gene_caller_id") %>% 
  mutate(Type = case_when(log2FoldChange>1 ~ "Jelly",
                          log2FoldChange< -1 ~ "Control")) %>% 
  filter(Type %in% c("Jelly","Control"))


#combine with the definitions of the modules and parse according to the different modules
#(that means that each protein will appear in a new row for each module!)
metaG_KEGG_modules <- read.csv(paste(wd,"metaG_analysis/metaG_anvio/05_ANVIO/spades-Kofam_modules.txt",sep=""), sep="\t", h= T)

metaP.DEseq.res_KOfam_by_module <- metaP.DEseq.res_KOfam %>%
  tidyr::separate_rows(modules_with_ko, sep = ',') %>% 
  dplyr::rename("kegg_module" = "modules_with_ko") %>% 
  left_join(metaG_KEGG_modules, by = c("contig_name","kegg_module"))

#check how many proteins were enriched per module
enr_prot_per_KEGG_module<- metaP.DEseq.res_KOfam_by_module %>% 
  group_by(Type, module_name,module_class, module_category, module_subcategory, module_definition) %>% 
  summarize(Prot.n = n())


#plot only pathway modules
enr_prot_per_KEGG_module_pathway<- enr_prot_per_KEGG_module %>% 
  filter(module_class =="Pathway modules", !is.na(module_category))

#plot
enr_prot_KEGG_modules.p <- ggplot(enr_prot_per_KEGG_module_pathway, aes(x= module_category, y= Prot.n, fill = module_subcategory))+
  geom_bar(position="stack", stat="identity")+
  facet_grid(.~Type, space= "fixed") +
  scale_fill_manual(values = tol21rainbow)+ 
  #guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Number of proteins \n")+
  geom_hline(aes(yintercept=-Inf)) + 
  geom_vline(aes(xintercept=-Inf)) +
  geom_vline(aes(xintercept=Inf))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),axis.text.x = element_text(angle = 90),
        text=element_text(size=14),legend.position = "bottom", 
        axis.title.x = element_blank())


#save the plot
ggsave(paste0(wd,"/R_figures/metaP_KEGG_modules.pdf"), 
       plot = enr_prot_KEGG_modules.p,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)

###################
#Identify which enriched proteins were associated with Bins
###################
Refined_DAS_bins <-  read.csv(paste(wd,"metaG_analysis/metaG_anvio/06_BINS/Refined_bins_collection.txt",sep=""),
                              sep="\t", h= F)

names(Refined_DAS_bins)<- c("contig","Bin")

metaP.DEseq.res.sig_Bins <- metaP.DEseq.res.sig %>% 
  mutate(Type = case_when(log2FoldChange>1 ~ "Jelly",
                          log2FoldChange< -1 ~ "Control")) %>% 
  filter(Type %in% c("Jelly","Control")) %>% 
  left_join(Refined_DAS_bins, by = c("contig")) %>% 
  filter(!is.na(Bin))

#summarize per bin
metaP.DEseq.res.sig_Bins %>% group_by(Bin, Type) %>% 
  summarize(Total = n())

#explore proteins enriched in Bin_84
Bin_84_enr_prot. <- metaP.DEseq.res.sig_Bins %>% 
  filter(Bin =="Bin_84", Type =="Jelly")

####use input from Bin_exploration.R for KEGG pathways
test <- Bins_gene_calls_KEGG_modules %>% 
  filter(gene_caller_id %in% intersect(Bins_gene_calls_KEGG_modules$gene_caller_id, Bin_84_enr_prot.$gene_caller_id))


test <- Bins_gene_calls_KEGG_modules %>%  
  filter(db_name=="Bin_84") %>% 
  mutate(colour = "yellow") %>% 
  select("ko","colour") %>% 
  unique() 


test %>% 
  mutate(colour = case_when(ko %in% Bin_84_enr_prot.$KeggGhostKoala_accession ~ "red",
                            TRUE~ colour))%>% 
  write.table("data/KEGG/Bin_84_ko.txt",
              row.names = FALSE,
              col.names = FALSE,quote = FALSE)

Bin_84_enr_prot. %>% 
  select(KeggGhostKoala_accession) %>% 
  mutate(colour = "orange",
         ko= KeggGhostKoala_accession) %>% 
  select("ko","colour") %>% 
  unique() %>% 
  bind_rows(test) %>% 
  write.table("data/KEGG/Bin_84_ko_prot.txt",
              row.names = FALSE,
              col.names = FALSE,quote = FALSE)
