

Refined_DAS_bins <-  read.csv(paste(wd,"metaG_analysis/metaG_anvio/06_BINS/Refined_bins_collection.txt",sep=""),
                              sep="\t", h= F)

names(Refined_DAS_bins)<- c("contig","Bin")

metaP.DEseq.res.sig_Bins <- metaP.DEseq.res.sig %>% 
  left_join(Refined_DAS_bins, by = c("contig")) %>% 
  filter(!is.na(Bin))



metaP.DEseq.res <- as(metaP.DEseq.res, "data.frame") %>%  mutate(gene_caller_id = row.names(.)) %>% 
  left_join(as.data.frame(tax_table(metaP_obj0)), by = "gene_caller_id")


test <- Bins_gene_calls_KEGG_modules %>% 
  mutate(gene_caller_id = as.character(gene_caller_id)) %>% 
  left_join(metaP.DEseq.res.sig[,c("gene_caller_id","contig","log2FoldChange",
                                   "COG20_FUNCTION_function","Pfam_function","InterPro_function")],
            by = c("gene_caller_id","contig")) %>% 
  filter(!is.na(log2FoldChange))
