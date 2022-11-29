require(pheatmap)
require(phyloseq)
require(DESeq2)
require(dplyr)
require(ggplot2)

#calculate standard error
se <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}

tol21rainbow<- c("#771155", "#AA4488","#CC99BB","#114477", 
                 "#4477AA","#117744","#117777","#88CCAA", 
                 "#77CCCC","#00ffff","#44AA77","#44AAAA", 
                 "#777711","#AAAA44","#DDDD77","#774411", 
                 "#AA7744","#DDAA77","#771122","#AA4455", "#DD7788"
)

#load metaproteome phyloseq object
metaP_runB<- readRDS("data/metaproteome/metaP_ps_runB.rds")

###################
#Identify sig. enriched proteins between Jelly and Control bottles in each fraction separately
###################
DESeq_res<- data.frame()

for (i in c("MP","EH","EL")){
metaP_runB_sub <- subset_samples(metaP_runB, Treatment %in% c("Control","Jelly") & Fraction == i) %>% 
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
meta_df <- as.data.frame(colData(metaP.DEseq)[,c("Replicate","Treatment")])

#produce and save a heatmap
pheatmap(assay(metaP.DEseq.vsd)[select_prot,], cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=meta_df,
         main = "Top 100 proteins", filename = paste0("Figures/",i,"_heatmap.pdf"))


#export results for Jelly vs. Control
metaP.DEseq.res <- results(metaP.DEseq, name = "Treatment_Jelly_vs_Control")

#export results for Jelly vs. Control with shrinked LFC estimates
#(explained: http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#altshrink)
metaP.DEseq.res_LFC <- as(lfcShrink(metaP.DEseq, coef="Treatment_Jelly_vs_Control", type="apeglm"), "data.frame")%>% 
                        mutate(gene_callers_id = as.numeric(row.names(.)),
                                Fraction = i) 

metaP.annotations <- as.data.frame(tax_table(metaP_runB_sub)) %>% 
                      mutate(gene_callers_id = as.numeric(gene_callers_id))

#extract only significant proteins and merge with annotations
metaP.DEseq.res.sig <- metaP.DEseq.res_LFC %>% 
                      filter(padj < 0.1 ) %>% 
                      mutate(Type = case_when(log2FoldChange>1 ~ "Jelly",
                                              log2FoldChange< -1 ~ "Control",
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
#summarize the significantly enr. proteins by Taxonomic orders
###################
#number of proteins per family - filter out those with < 10 proteins
metaP.DEseq.res.sig.Tax_top <- DESeq_res %>% 
  group_by(Fraction, Type, Class, Order, Family) %>% 
  summarize(Total_prot.= length(Type)) %>% 
  filter(Total_prot.>10)

top_fam <- as.vector(unique(metaP.DEseq.res.sig.Tax_top$Family))

#filter only families of interest and prune annotations
DESeq_res_top_fam <- DESeq_res %>% 
                      filter(Family %in% top_fam, Type %in% c("Jelly","Control")) %>% 
                      group_by(Fraction, Type, Class, Order, Family, COG20_CATEGORY_function, COG20_CATEGORY_accession) %>% 
                      mutate(COG20_CATEGORY_function = gsub("!!!.*","",COG20_CATEGORY_function),
                             COG20_CATEGORY_accession = gsub("!!!.*","",COG20_CATEGORY_accession),
                             Fraction = factor(Fraction, levels = c("MP","EH","EL"))) 


#aggregate by category
DESeq_res_top_fam_agg <- DESeq_res_top_fam %>% 
  group_by(Fraction, Type, Class, Order, Family) %>% 
  summarize(log2_mean = mean(log2FoldChange), log2_se = se(log2FoldChange), n_prot= length(log2FoldChange)) %>% 
  filter(n_prot> 2)


#plot
ggplot(DESeq_res_top_fam_agg, aes(x= Family, y= log2_mean, label = n_prot))+
  geom_errorbar(aes(ymin= log2_mean-log2_se, ymax= log2_mean+log2_se), width = 0.2)+
  geom_point(aes(colour = Class, shape = Type), size = 3)+
  geom_text(nudge_x = 0.3, colour = "gray50")+
  facet_grid(.~Fraction, space= "fixed") +
  scale_colour_manual(values = tol21rainbow)+ 
  #guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Mean log2 fold change \n")+
  geom_hline(aes(yintercept=0), linetype = 2, alpha = 0.3) +
  geom_hline(aes(yintercept=-Inf)) + 
  geom_vline(aes(xintercept=-Inf)) +
  geom_vline(aes(xintercept=Inf))+
  coord_flip()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),axis.text.x = element_text(angle = 90),
        #text=element_text(size=14),
        #legend.position = "bottom", 
        #axis.title.x = element_blank()
        )

#save the plot
ggsave("./Figures/metaP_log2foldchange_per_taxa.pdf", 
       plot = last_plot(),
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)


#########################################################
#Explore KEGG pathways with enr. proteins             ###
#########################################################
DESeq_res_sub <- DESeq_res %>%  filter(Type =="Jelly", KOfam_accession!= "Unk") %>% 
  mutate(KOfam_accession = gsub("\\|.*","",KOfam_accession)) %>% select(KOfam_accession,log2FoldChange)%>% 
  tibble::deframe()

#pathways of interest
pathways<- c("00010", #Glycolysis / Gluconeogenesis
             "00020", #TCA cycle
             "01200", #carbon metabolism
             "00190", #Oxidative phosphorylation
             "00071", #Fatty acid degradation
             "00250", #Alanine, aspartate and glutamate metabolism
             "00260", #Glycine, serine and threonine metabolism
             "00270", #Cysteine and methionine metabolism
             "00280", #Valine, leucine and isoleucine degradation
             "00310", #Lysine degradation
             "00330", #Arginine and proline metabolism
             "00340", #Histidine metabolism
             "00350", #Tyrosine metabolism
             "00360", #Phenylalanine metabolism
             "00380", #Tryptophan metabolism
             "00400" #Phenylalanine, tyrosine and tryptophan biosynthesis
)

#plot the pathways of interest for all enriched proteins
pathview(gene.data = DESeq_res_sub, 
         pathway.id =pathways,
         species = "ko", 
         keys.align = "y", 
         kegg.native = T,
         #map.null=FALSE,
         #both.dirs = TRUE,
         #discrete	=list(gene=TRUE),
         res = 300, cex = 0.25,
         out.suffix = "all_proteins",
         kegg.dir="Figures/KEGG/")
