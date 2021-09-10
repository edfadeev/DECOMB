library(pathview)
#simulate data first, gene hits are 1's, others are NA's
gnames <- sim.mol.data(nmol = 10000, discrete=T)
gdata=matrix(NA, ncol=4, nrow=10000)
rownames(gdata)=gnames
idx1=sample(1:10000, 4000)
idx2=sample(1:10000, 4000)
idx3=sample(1:10000, 4000)
idx4=sample(1:10000, 4000)
gdata[idx1,1]=1
gdata[idx2,2]=1
gdata[idx3,3]=1
gdata[idx4,4]=1


test <- Bins_kofam_hits%>% 
  filter(genome_name == "Bin_115_2")

#parse by modules and add grouping according to enr. protein/genomic data
Bins.DEseq.enr.prot_by_module <- Bins.DEseq.enr.prot %>% 
  tidyr::separate_rows(modules_with_ko, sep = ',') %>% 
  mutate(Control = ifelse(log2FoldChange<0, 1,NA),
         Jelly = ifelse(log2FoldChange>0, 1,NA),
         gene = 1)


Bins.DEseq.enr.prot_by_module_sub<- Bins.DEseq.enr.prot_by_module %>% 
                                      filter(modules_with_ko == "M00001",
                                             genome_name == "Bin_84_1") %>% 
                                      select(-c("unique_id","contig.x","gene_caller_id")) %>% 
                                        unique()



#donwload, parse and map your data
Bins.DEseq.enr.prot_sub<-Bins.DEseq.enr.prot %>% 
  filter(genome_name == "Bin_115_2") %>% 
  select(ko, log2FoldChange) %>% 
  mutate(Control = ifelse(log2FoldChange<0, 1,NA),
         Jelly = ifelse(log2FoldChange>0, 1,NA),
         gene = 1)
                             


download.kegg(pathway.id = "00010", species = "ko", kegg.dir = "./data/KEGG/",
              file.type=c("xml", "png"))

xml.file="./data/KEGG/ko00010.xml"
node.data=node.info(xml.file)
names(node.data)

plot.data.gene=node.map(mol.data=Bins.DEseq.enr.prot_sub, node.data, node.types="gene")

head(plot.data.gene)

#generate and customize node colors
cols.ts.gene=node.color(plot.data.gene, limit=2, bins=4)
head(cols.ts.gene)
cols=rainbow(4)
bg.col="#FFFFFF"
for(i in 1:4){
  cols.ts.gene[cols.ts.gene[,i]!=bg.col,i]=cols[i]
}

#KEGG view
pv.pars= keggview.native(plot.data.gene=plot.data.gene,
                         cols.ts.gene=cols.ts.gene, node.data=node.data, pathway.name="hsa04110",
                         same.layer=T, plot.col.key=F, out.suffix = "custom.4cols")
