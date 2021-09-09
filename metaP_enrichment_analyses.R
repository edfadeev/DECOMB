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

#load metaproteome phyloseq object
metaP_obj0<- readRDS("data/metaP_ps_raw.rds")
exoP_obj0<- readRDS("data/exoP_ps_raw.rds")

###################
#Identify enriched proteins in MetaP
###################
#run DESeq2
metaP_obj0.ddsMat <- phyloseq_to_deseq2(metaP_obj0, ~Type)
metaP.DEseq <- DESeq(metaP_obj0.ddsMat, fitType="local")

#conduct variance transformation
metaP.DEseq.vsd <- vst(metaP_obj0.ddsMat,  fitType="local", blind=FALSE)

#select only the top 50 proteins and produce dataframe with metadata
select_prot <- order(rowMeans(counts(metaP.DEseq, normalized=TRUE)),
                decreasing=TRUE)[1:50]
meta_df <- as.data.frame(colData(metaP.DEseq)[,c("Replicate","Type")])

#produce and save a heatmap
pheatmap(assay(metaP.DEseq.vsd)[select_prot,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=meta_df,
         main = "Top 50 proteins", filename = paste0(wd, "R_figures/metaP_heatmap.pdf"))






metaP_obj0_no_T0<- subset_samples(metaP_obj0, Type != "T0")
#remove proteins that were not observed 
metaP_obj0_no_T0<- prune_taxa(taxa_sums(metaP_obj0_no_T0)>0,metaP_obj0_no_T0)









sampleDists <- dist(t(assay(vsd)))


library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$Type
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)



plotPCA(vsd, intgroup=c("Type"),  returnData=TRUE)




#extract only significant proteins
metaP.DEseq.res.sig <- as(metaP.DEseq.res, "data.frame") %>%  filter(padj < 0.1)

metaP.DEseq.res.sig <- cbind(as(metaP.DEseq.res.sig, "data.frame"),
                             as(tax_table(metaP_obj0)[rownames(metaP.DEseq.res.sig), ], "matrix"))



metaP_obj0_sub <- prune_taxa(rownames(log2FoldChange_proteins), metaP_obj.vst)

#plot heatmap
pheatmap(otu_table(metaP_obj0_sub), labels_row = tax_table(metaP_obj0_sub)[,"Class"])




metaP.DEseq.res <- results(metaP.DEseq, name = "Type_Jelly_vs_Control")

###################
#Import KEGG modules estimates for each gene
###################
metaG_kofam_hits <- read.csv(paste(wd,"metaG_analysis/metaG_anvio/05_ANVIO/spades_kofam_hits.txt",sep=""), sep="\t", h= T)

metaG_KEGG_modules <- read.csv(paste(wd,"metaG_analysis/metaG_anvio/05_ANVIO/spades_modules.txt",sep=""), sep="\t", h= T)

metaP_sig_prot_by_module <- merge(metaP.DEseq.res.sig, metaG_kofam_hits, by ="gene_caller_id") %>%
  tidyr::separate_rows(modules_with_ko, sep = ',') %>% 
  dplyr::rename("kegg_module" = "modules_with_ko") %>% 
  merge(metaG_KEGG_modules, by = "kegg_module")


enriched_prot.p <- ggplot(metaP_sig_prot_by_module, aes(x=module_name, y=log2FoldChange, 
                                     colour= module_category))+
  geom_point(size = 2)+
  scale_color_manual(values = tol21rainbow)+
  coord_flip()+
  theme_bw()+
  theme(text=element_text(size=14),legend.position = "right")



#save the plot
ggsave(paste0(wd,"/R_figures/enriched_prot.pdf"), 
       plot = enriched_prot.p,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)


###################
#Import metabolism estimates for the bins
###################
modules_info <- read.csv(paste(wd,"metaG_analysis/metaG_anvio/07_METABOLISM/modules_info.txt",sep=""), sep="\t", h= T)

Bins_kofam_hits <- read.csv(paste(wd,"metaG_analysis/metaG_anvio/07_METABOLISM/Bins_kofam_hits.txt",sep=""), sep="\t", h= T)

Bins_modules <- read.csv(paste(wd,"metaG_analysis/metaG_anvio/07_METABOLISM/Bins_modules.txt",sep=""), sep="\t", h= T)


#merge the identified modules and the genes calls
Bins_gene_calls_KEGG_modules<-Bins_modules %>% 
  tidyr::separate_rows(gene_caller_ids_in_module, sep = ',') %>% 
  dplyr::rename("gene_caller_id" = "gene_caller_ids_in_module") %>% 
  mutate(gene_caller_id=as.integer(gene_caller_id),
         genome_name = factor(genome_name, levels = c("Bin_84_1","Bin_76_1","Bin_5_2",
                                                      "Bin_2_1","Bin_5_3","Bin_38_1","Bin_2_2",
                                                      "Bin_179_1","Bin_115_2","Bin_115_1","Bin_102_1",
                                                      "Bin_12_1","Bin_150_1_1"))) %>% 
  left_join(Bins_kofam_hits_parsed[,c("genome_name","gene_caller_id","contig","kegg_module","ko","ko_definition")],
            by = c("genome_name","kegg_module","gene_caller_id"))


#merge the enriched proteins with the bins
metaP.DEseq.res.sig <- metaP.DEseq.res.sig  %>%
  mutate(gene_caller_id = as.integer(gene_caller_id)) %>% 
  left_join(., Bins_gene_calls_KEGG_modules, by = "gene_caller_id") %>% 
  select(-c("start","stop", "direction", "call_type", "source","unique_id","version","contig.y"))


write.table(metaP.DEseq.res.sig, "metaP/metaP_DEseq_res_sig.txt")

#see how many enriched proteins were associated with bins
metaP.DEseq.res.sig %>% group_by(genome_name) %>% 
  dplyr::summarize(num_of_prot= n()) %>% 
  dplyr::arrange(desc(num_of_prot))


#subset each bins proteins and summarize according to KO modules 
#Bin 84_1 - Pseudoalteromonas phenolica - length 3.92Mbp (C93/R0)
metaP.DEseq.res.sig_Bin_P_phenolica<-metaP.DEseq.res.sig %>% filter(genome_name == "Bin_84_1") 


Bins_gene_calls_KEGG_modules_Bin_84_1<- Bins_gene_calls_KEGG_modules %>% 
  mutate(gene_caller_id = as.character(gene_caller_id)) %>% 
  left_join(metaP.DEseq.res.sig, by = "gene_caller_id") %>% 
  filter(genome_name == "Bin_84_1") 

#export for KEGG website
Bins_gene_calls_KEGG_modules_Bin_84_1_for_KEGG <- Bins_gene_calls_KEGG_modules_Bin_84_1 %>%  
  mutate(colour = case_when(log2FoldChange>0 ~"red", 
                            log2FoldChange<0 ~"blue",
                            is.na(log2FoldChange)~"orange")) %>% 
  select("ko","colour")
  
write.table(Bins_gene_calls_KEGG_modules_Bin_84_1_for_KEGG, "Bin_84_KEGG.txt",
            row.names = FALSE,
            col.names = FALSE,quote = FALSE)







###################
#Identify enriched proteins in MetaP
###################
metaP_obj0_no_T0<- subset_samples(exoP_obj0, Type != "T0")
#remove proteins that were not observed 
metaP_obj0_no_T0<- prune_taxa(taxa_sums(metaP_obj0_no_T0)>0,metaP_obj0_no_T0)

metaP_obj0.ddsMat <- phyloseq_to_deseq2(metaP_obj0_no_T0, ~Type)
metaP_obj0.ddsMat = estimateSizeFactors(metaP_obj0.ddsMat)
metaP_obj0.ddsMat <- estimateDispersions(metaP_obj0.ddsMat)
metaP.DEseq <- DESeq(metaP_obj0.ddsMat, fitType="local")
metaP.DEseq.res <- results(metaP.DEseq)

#extract only significant proteins
metaP.DEseq.res.sig <- as(metaP.DEseq.res, "data.frame") %>%  filter(padj < 0.1)

metaP.DEseq.res.sig <- cbind(as(metaP.DEseq.res.sig, "data.frame"),
                             as(tax_table(exoP_obj0)[rownames(metaP.DEseq.res.sig), ], "matrix"))




metaP.DEseq.res.sig_Bin_84_1<-metaP.DEseq.res.sig %>% filter(genome_name == "Bin_84_1")

Bin_84_1_logFC <- metaP.DEseq.res.sig_Bin_84_1$log2FoldChange
names(Bin_84_1_logFC)<- metaP.DEseq.res.sig_Bin_84_1$ko


pv.out <- pathview(gene.data = Bin_84_1_logFC, 
                   pathway.id ="01230",
                   species = "hsa",
                   kegg.native = T,
                   keys.align = "y", 
                   gene.idtype="KEGG",
                   low = "blue", mid = "gray", high = "#ffa500", bin = 30,
                   out.suffix = "Bin_84_1")




#explore metabolic capacities
test <- Bins_modules %>% filter(module_completeness > 0.5) %>% group_by(genome_name, module_subcategory) %>% 
  summarize(total= n())


ggplot(test, aes(x = genome_name, y = total))+
  geom_col()+
  facet_wrap(module_subcategory~., scales = "free_y")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(size=14),legend.position = "bottom", 
        axis.title.x = element_blank(), axis.text.x = element_text(angle=90))





#subset each bins proteins and summarize according to KO modules 
#Bin 84_1 - Pseudoalteromonas phenolica - length 3.92Mbp (C93/R0)
metaP.DEseq.res.sig_Bin_P_phenolica<-metaP.DEseq.res.sig %>% filter(genome_name == "Bin_84_1") 

#check how many proteins in each KO module
metaP.DEseq.res.sig_Bin_P_phenolica %>% 
  separate_rows(modules_with_ko, sep = ',') %>% 
  group_by(modules_with_ko) %>% summarize(num_of_prot= n()) %>% 
  arrange(desc(num_of_prot))

#check for potentially duplicated rows as a result of different potential KO annotations (will be resolved at a later stage)
metaP.DEseq.res.sig_Bin_P_phenolica[duplicated(metaP.DEseq.res.sig_Bin_P_phenolica$start),]




metaP.DEseq.res.sig_Bin_Vibrio<-metaP.DEseq.res.sig %>% filter(genome_name == "Bin_38_1") 

metaP.DEseq.res.sig_Bin_Bermanella<-metaP.DEseq.res.sig %>% filter(genome_name == "Bin_115_1") 

metaP.DEseq.res.sig_Bin_Glaciecola<-metaP.DEseq.res.sig %>% filter(genome_name == "Bin_115_2") 

metaP.DEseq.res.sig_Bin_Alteromonas<-metaP.DEseq.res.sig %>% filter(genome_name %in% c("Bin_2_1","Bin_2_2"))

#not assigned to any bin
metaP.DEseq.res.sig_no_bin<-metaP.DEseq.res.sig %>% filter(is.na(genome_name) == TRUE) 


#split posible annotations into separated rows
test <- prot_frac_agg_metaP.DEseq.res.sig_no_bin 
