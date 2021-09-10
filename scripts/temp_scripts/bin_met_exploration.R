



###################
#Combine the enriched proteins with the metabolism estimates for the bins
###################
modules_info <- read.csv(paste(wd,"metaG_analysis/metaG_anvio/07_METABOLISM/modules_info.txt",sep=""), sep="\t", h= T)

Bins_kofam_hits <- read.csv(paste(wd,"metaG_analysis/metaG_anvio/07_METABOLISM/Bins_kofam_hits.txt",sep=""), sep="\t", h= T)

Bins_modules <- read.csv(paste(wd,"metaG_analysis/metaG_anvio/07_METABOLISM/Bins_modules.txt",sep=""), sep="\t", h= T)




Bins.DEseq.enr.prot <- Bins_kofam_hits %>% 
                        mutate(gene_caller_id = as.character(gene_caller_id)) %>% 
                        left_join(metaP.DEseq.res.sig, by = "gene_caller_id")


#parse by modules and add grouping according to enr. protein/genomic data
Bins.DEseq.enr.prot_by_module <- Bins.DEseq.enr.prot %>% 
  tidyr::separate_rows(modules_with_ko, sep = ',') %>% 
  mutate(Type = case_when(log2FoldChange>0 ~ "Jelly",
                          log2FoldChange<0 ~ "Control",
                          is.na(log2FoldChange) ~ "gene"))







#########################################################################
#work in progress
#########################################################################



############################################################################
#plot the pathways of interest for each bin
############################################################################
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
             "00400", #Phenylalanine, tyrosine and tryptophan biosynthesis
             "00430",# Taurine and hypotaurine metabolism
)



#plot in a loop all the pathways of interest for each bin
for (bin in unique(Bins.DEseq.enr.prot_by_module$genome_name)) {
  Bins_gene_calls_KEGG_modules_complete<- Bins.DEseq.enr.prot_by_module %>%  
    filter(genome_name == bin & module_completeness >0.7) 
  
  pathview(gene.data = Bins_gene_calls_KEGG_modules_complete$ko, 
           pathway.id =pathways,
           species = "ko", 
           keys.align = "y", 
           kegg.native = T, both.dirs = TRUE,
           discrete	=list(gene=TRUE),
           res = 300, cex = 0.25,
           out.suffix = bin)
}






























Bins.DEseq.enr.prot_by_module_bin<-Bins.DEseq.enr.prot_by_module %>% 
                      filter(genome_name == "Bin_115_2") %>% 
                        mutate(color = case_when(Type =="Jelly" ~ "red",
                                                 Type =="Control" ~ "blue",
                                                 Type =="gene" ~ "yellow")) %>% 
                        select(ko, color) %>% 
                        unique()


pv.out <- pathview(gene.data = Bins.DEseq.enr.prot_by_module_bin, 
                   pathway.id ="01230",
                   species = "ko",
                   kegg.native = T,
                   keys.align = "y", 
                   gene.idtype="KEGG",
                   low = "blue", mid = "gray", high = "#ffa500", bin = 30,
                   out.suffix = "Bin_84_1")





#merge the identified modules and the genes calls
Bins_gene_calls_KEGG_modules<-Bins_modules %>% 
  tidyr::separate_rows(gene_caller_ids_in_module, sep = ',') %>% 
  dplyr::rename("gene_caller_id" = "gene_caller_ids_in_module") %>% 
  filter(module_completeness >0.7) %>% 
  mutate(genome_name = factor(genome_name, levels = c("Bin_84_1","Bin_76_1","Bin_5_2",
                                                      "Bin_2_1","Bin_5_3","Bin_38_1","Bin_2_2",
                                                      "Bin_179_1","Bin_115_2","Bin_115_1","Bin_102_1",
                                                      "Bin_12_1","Bin_150_1_1"))) 




Bins_gene_calls_KEGG_modules_enr_prot <- Bins_gene_calls_KEGG_modules %>% 
                        left_join(metaP.DEseq.res_KOfam, by = "gene_caller_id")
  
  
