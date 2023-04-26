require(pheatmap)
require(phyloseq)
require(DESeq2)
require(dplyr)
require(ggplot2)

source("scripts/extra_functions.R")

#load metaproteome phyloseq object
metaP_merged<- readRDS("data/metaproteome/metaP_runB_merged.rds")

###################
#Identify sig. enriched proteins between Jelly and Control bottles in each fraction separately
###################
DESeq_res<- data.frame()

for (i in c("Cellular","Exocellular")){
metaP_runB_sub <- subset_samples(metaP_merged, Treatment %in% c("Control","Cteno-OM") & Type == i) %>% 
                    prune_taxa(taxa_sums(.)>0,.)

#run DESeq2 enrichment test
metaP_runB.ddsMat <- phyloseq_to_deseq2(metaP_runB_sub, ~Treatment)
metaP.DEseq.vsd <- vst(metaP_runB.ddsMat,  fitType="local", blind=FALSE)

metaP.DEseq <- DESeq(metaP_runB.ddsMat, fitType="local")

###################
#Generate heatmap of the top 100 proteins
###################
#select only the top 100 proteins and produce dataframe with metadata
select_prot <- order(rowMeans(counts(metaP.DEseq, normalized=FALSE)),
                     decreasing=TRUE)[1:100]
meta_df <- as.data.frame(colData(metaP.DEseq)[,c("Sample_name","Treatment")])

#produce and save a heatmap
pheatmap(assay(metaP.DEseq.vsd)[select_prot,], cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=meta_df,
         main = "Top 100 proteins", filename = paste0("Figures/",i,"_heatmap.pdf"))


#export results for Jelly vs. Control
metaP.DEseq.res <- results(metaP.DEseq, contrast = c("Treatment","Cteno-OM","Control"))

#export results for Jelly vs. Control with shrinked LFC estimates
#(explained: http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#altshrink)
metaP.DEseq.res_LFC <- as(lfcShrink(metaP.DEseq, coef="Treatment_Cteno.OM_vs_Control", type="apeglm"), "data.frame")%>% 
                        mutate(gene_callers_id = as.numeric(row.names(.)),
                                Fraction = i) 

metaP.annotations <- as.data.frame(tax_table(metaP_runB_sub)) %>% 
                      mutate(gene_callers_id = as.numeric(gene_callers_id))

#extract only significant proteins and merge with annotations
metaP.DEseq.res.sig <- metaP.DEseq.res_LFC %>% 
                      filter(padj < 0.1 ) %>% 
                      mutate(Type = case_when(log2FoldChange>0 ~ "Cteno-OM",
                                              log2FoldChange<0 ~ "Control",
                                              TRUE ~ "Not.sig.")) %>% 
                      left_join(metaP.annotations, by = "gene_callers_id")

DESeq_res <- rbind(DESeq_res,metaP.DEseq.res.sig)
}

write.csv(DESeq_res, "data/DESEq_res.csv")

#total of sig. proteins in each treatment
DESeq_res.overview<- DESeq_res %>% 
  group_by(Fraction, Type) %>% 
  summarize(Total_prot.= length(Type))


###################
#Explore the significantly enr. proteins
###################
DESeq_res<- read.csv("data/DESEq_res.csv")

DESeq_res_Cteno_cell <- DESeq_res %>% 
                  filter(Fraction =="Cellular",
                         log2FoldChange >15) %>% 
                  select(c("InterPro_function","KeggGhostKoala_function",
                         "KeggGhostKoala_accession","KOfam_accession","KOfam_function","Family", "log2FoldChange"))


DESeq_res_Cteno_Exo <- DESeq_res %>% 
  filter(Fraction =="Exocellular",
         #log2FoldChange >15
         ) %>% 
  select(c("InterPro_function","KeggGhostKoala_function",
           "KeggGhostKoala_accession","KOfam_accession","KOfam_function","Family", "log2FoldChange"))

###################
#summarize the significantly enr. proteins by Taxonomic orders
###################
#number of proteins per family 
metaP.DEseq.res.sig.Tax_top <- DESeq_res %>% 
  group_by(Fraction, Type, Class, Order, Family, Genus) %>% 
  summarize(Total_prot.= length(Fraction))

top_fam <- as.vector(unique(metaP.DEseq.res.sig.Tax_top$Family))

#filter only families of interest and prune annotations
DESeq_res_top_fam <- DESeq_res %>% 
                      filter(Family %in% top_fam, Type %in% c("Cteno-OM","Control")) %>% 
                      group_by(Fraction, Type, Class, Order, Family, Genus, COG20_CATEGORY_function, COG20_CATEGORY_accession) %>% 
                      mutate(COG20_CATEGORY_function = gsub("!!!.*","",COG20_CATEGORY_function),
                             COG20_CATEGORY_accession = gsub("!!!.*","",COG20_CATEGORY_accession),
                             Fraction = factor(Fraction, levels = c("Cellular", "Exocellular"))) 


#aggregate by category
DESeq_res_top_fam_agg <- DESeq_res_top_fam %>% 
  group_by(Fraction, Type, Class, Order, Family, Genus,) %>% 
  summarize(log2_mean = mean(log2FoldChange), log2_se = se(log2FoldChange), n_prot= length(log2FoldChange)) %>% 
  filter(n_prot> 4, !Genus %in% grep("Unknown.*", Genus, value = T))


#plot
DESeq_res_top_fam_agg.p <- ggplot(DESeq_res_top_fam_agg, aes(x= Genus, y= log2_mean, label = n_prot))+
  geom_errorbar(aes(ymin= log2_mean-log2_se, ymax= log2_mean+log2_se), width = 0.2)+
  geom_point(aes( shape = Type), size = 6, colour = "black")+
  geom_point(aes(colour = Order, shape = Type), size = 4)+
  geom_text(nudge_x = 0.5, colour = "gray50")+
  facet_grid(.~Fraction, space= "fixed") +
  scale_colour_manual(values = tol21rainbow)+ 
  #guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Mean log2 fold change \n")+
  xlab("Genus (> 5 proteins)")+
  geom_hline(aes(yintercept=0), linetype = 2, alpha = 0.3) +
  coord_flip()+
  theme_EF+
  theme(legend.position = "bottom")

#save the plot
ggsave("./Figures/Figure_5-metaP_log2foldchange.pdf", 
       plot = DESeq_res_top_fam_agg.p,
       units = "mm",
       width = 90, height = 90, 
       scale = 3,
       dpi = 300)


#########################################################
#Explore enr. proteins by taxa in all fractions       ###
#########################################################
#Pseudoalteromonas
DESeq_res_top_Pseudoalt <- DESeq_res_top_fam %>% 
  filter(Genus == "Pseudoalteromonas")

#identify the faunction categories with most of the proteins
DESeq_res_top_Pseudoalt_total_COG20_category<- DESeq_res_top_Pseudoalt%>% 
  select(gene_callers_id, Type, COG20_CATEGORY_function, COG20_CATEGORY_accession) %>% 
  unique() %>% 
  group_by(Type, COG20_CATEGORY_function, COG20_CATEGORY_accession) %>% 
  summarize(n_prot = n())

#identify the COGs within each category
DESeq_res_top_Pseudoalt_total_COG20_function <- DESeq_res_top_Pseudoalt%>% 
  select(gene_callers_id, Fraction, Type, COG20_CATEGORY_accession, COG20_CATEGORY_function, COG20_FUNCTION_function, COG20_FUNCTION_accession,log2FoldChange) %>% 
  unique() %>% 
  group_by(Type,Fraction, COG20_CATEGORY_accession, COG20_CATEGORY_function, COG20_FUNCTION_function, COG20_FUNCTION_accession) %>% 
  summarize(log2_mean = mean(log2FoldChange), log2_se = se(log2FoldChange), n_prot= length(log2FoldChange)) 


#Alteromonas
DESeq_res_top_Alteromonas <- DESeq_res_top_fam %>% 
  filter(Genus == "Alteromonas")

DESeq_res_top_Alteromonas_total_COG20_category<- DESeq_res_top_Alteromonas%>% 
  select(gene_callers_id, Type, COG20_CATEGORY_function, COG20_CATEGORY_accession) %>% 
  unique() %>% 
  group_by(Type, COG20_CATEGORY_function, COG20_CATEGORY_accession) %>% 
  summarize(n_prot = n())

DESeq_res_top_Alteromonas_total_COG20_function <- DESeq_res_top_Alteromonas%>% 
  select(gene_callers_id, Fraction, Type, COG20_CATEGORY_accession, COG20_CATEGORY_function, COG20_FUNCTION_function, COG20_FUNCTION_accession,log2FoldChange) %>% 
  unique() %>% 
  group_by(Type, Fraction, COG20_CATEGORY_accession, COG20_CATEGORY_function, COG20_FUNCTION_function, COG20_FUNCTION_accession) %>% 
  summarize(log2_mean = mean(log2FoldChange), log2_se = se(log2FoldChange), n_prot= length(log2FoldChange)) 


#vibrio
DESeq_res_top_Vibrio <- DESeq_res_top_fam %>% 
                          filter(Genus == "Vibrio") 


DESeq_res_top_Vibrio_total_COG20_category<- DESeq_res_top_Vibrio%>% 
  select(gene_callers_id, Type,  COG20_CATEGORY_function, COG20_CATEGORY_accession) %>% 
  unique() %>% 
  group_by(Type, COG20_CATEGORY_function, COG20_CATEGORY_accession) %>% 
  summarize(n_prot = n())

DESeq_res_top_Vibrio_total_COG20_function <- DESeq_res_top_Vibrio%>% 
  select(gene_callers_id, Type,  COG20_CATEGORY_accession,COG20_FUNCTION_function, COG20_FUNCTION_accession,log2FoldChange) %>% 
  unique() %>% 
  group_by(Type, COG20_CATEGORY_accession,COG20_CATEGORY_function, COG20_FUNCTION_function, COG20_FUNCTION_accession) %>% 
  summarize(log2_mean = mean(log2FoldChange), log2_se = se(log2FoldChange), n_prot= length(log2FoldChange)) 


#Pelagibacter
DESeq_res_top_Pelagi <- DESeq_res_top_fam %>% 
  filter(Family == "Pelagibacteraceae")

DESeq_res_top_Pelagi_total_COG20_function <- DESeq_res_top_Pelagi%>% 
  select(gene_callers_id, Type, Genus, COG20_FUNCTION_function, COG20_FUNCTION_accession,log2FoldChange) %>% 
  unique() %>% 
  group_by(Type, Genus,  COG20_CATEGORY_function, COG20_FUNCTION_function, COG20_FUNCTION_accession) %>% 
  summarize(log2_mean = mean(log2FoldChange), log2_se = se(log2FoldChange), n_prot= length(log2FoldChange)) 

DESeq_res_top_Pelagi_total_COG20_category<- DESeq_res_top_Pelagi%>% 
  select(gene_callers_id, Type, Genus, COG20_CATEGORY_function, COG20_CATEGORY_accession) %>% 
  unique() %>% 
  group_by(Type, Genus,  COG20_CATEGORY_function, COG20_CATEGORY_accession) %>% 
  summarize(n_prot = n())




#Flavobacteria
DESeq_res_top_Flavo <- DESeq_res_top_fam %>% 
  filter(Family == "Flavobacteriaceae")

DESeq_res_top_Flavo_total_COG20_function <- DESeq_res_top_Flavo%>% 
  select(gene_callers_id, Type, Genus, COG20_FUNCTION_function, COG20_FUNCTION_accession,log2FoldChange) %>% 
  unique() %>% 
  group_by(Type, Genus,  COG20_FUNCTION_function, COG20_FUNCTION_accession) %>% 
  summarize(log2_mean = mean(log2FoldChange), log2_se = se(log2FoldChange), n_prot= length(log2FoldChange)) 

DESeq_res_top_Flavo_total_COG20_category<- DESeq_res_top_Flavo%>% 
  select(gene_callers_id, Type, Genus, COG20_CATEGORY_function, COG20_CATEGORY_accession) %>% 
  unique() %>% 
  group_by(Type, Genus,  COG20_CATEGORY_function, COG20_CATEGORY_accession) %>% 
  summarize(n_prot = n())



#Rhodobacters
DESeq_res_top_Rhodo <- DESeq_res_top_fam %>% 
  filter(Family == "Rhodobacteraceae")
  
DESeq_res_top_Rhodo_total_COG20_function <- DESeq_res_top_Rhodo%>% 
  select(gene_callers_id, Type, Genus, COG20_FUNCTION_function, COG20_FUNCTION_accession,log2FoldChange) %>% 
  unique() %>% 
  group_by(Type, Genus,  COG20_FUNCTION_function, COG20_FUNCTION_accession) %>% 
  summarize(log2_mean = mean(log2FoldChange), log2_se = se(log2FoldChange), n_prot= length(log2FoldChange)) 

DESeq_res_top_Rhodo_total_COG20_category<- DESeq_res_top_Rhodo%>% 
  select(gene_callers_id, Type, Genus, COG20_CATEGORY_function, COG20_CATEGORY_accession) %>% 
  unique() %>% 
  group_by(Type, Genus,  COG20_CATEGORY_function, COG20_CATEGORY_accession) %>% 
  summarize(n_prot = n())

#print session info and clean the workspace
sessionInfo()
rm(list = ls())
gc()