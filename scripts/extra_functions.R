require(DESeq2)
library(phyloseq); packageVersion("phyloseq")

#define function for geometric mean
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}


#counts table geometric mean transformation
phyloseq_gm_mean_trans <- function(physeq){
  
  Dor_ps.dds <- phyloseq_to_deseq2(physeq, ~1)
  geoMeans = apply(counts(Dor_ps.dds), 1, gm_mean)
  Dor_ps.dds = estimateSizeFactors(Dor_ps.dds, geoMeans = geoMeans)
  Dor_ps.dds <- estimateDispersions(Dor_ps.dds)
  otu.vst <- getVarianceStabilizedData(Dor_ps.dds)
  
  Dor_ps.prev.vst<-physeq
  otu_table(Dor_ps.prev.vst)<- otu_table(otu.vst, taxa_are_rows = TRUE)
  
  Dor_ps.prev.vst
}

#scale parameters
scale_par <- function(x) scale(x, center = FALSE, scale = TRUE)[,1]

#define presence-absence function #https://github.com/hms-dbmi/UpSetR/issues/85
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

#calculate standard error
se <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}

#conduct NSAF transformation
#https://github.com/moldach/proteomics-spectralCount-normalization/blob/master/nsaf.R
#https://rdrr.io/github/DanielSprockett/reltools/man/add_nsaf.html
add_nsaf=function(ps, prot_length){
  if(ps@otu_table@taxa_are_rows == TRUE){
    mat <- (otu_table(ps))
  }else{
    mat <- t((otu_table(ps)))
  }
  prot_len <- as.numeric(unlist(as.numeric(tax_table(ps)[,prot_length]))) # Unlist your protein lengths before you sweep
  mat_prop <- sweep(mat,1,prot_len,"/") # Divide spectral counts (SpC) for a protein by its length (L)
  mat_sum <- as.data.frame(colSums(mat_prop)) # Get the column sums for each cell-line/treatment
  mat_sum <- mat_sum[,1]
  mat_nsaf <- sweep(mat_prop,2,mat_sum,"/") # Normalize by dividing by the sum of all SpC/L for all proteins identified 
  otu_table(ps) <- otu_table(mat_nsaf, taxa_are_rows = TRUE)
  return(ps)
}


#ggplot related functions
#theme
theme_EF <- theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        axis.text = element_text(size = 18),
        strip.text = element_text(size = 22, face ="bold"),
        axis.title = element_text(size = 20))

#colours
tol21rainbow<- c("#771155", "#AA4488","#CC99BB","#114477", 
                 "#4477AA","#117744","#117777","#88CCAA", 
                 "#77CCCC","#00ffff","#44AA77","#44AAAA", 
                 "#777711","#AAAA44","#DDDD77","#774411", 
                 "#AA7744","#DDAA77","#771122","#AA4455", "#DD7788"
)