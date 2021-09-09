##########################################
# Import and preprocess cell counts tables
##########################################
library(dplyr);packageVersion("dplyr")
library(tidyr);packageVersion("tidyr")

#microscope calculation factor
calc.factor <- 99515.5458411807

## Import raw counts
raw.counts <- read.csv("./FOV_all_groups_SH.csv",
                       sep = ",", dec = ".", header = TRUE)

#split the sample name column into different columns 



#calculate counts per FOV
raw.counts.SH <- raw.counts %>% 
  mutate(DAPI_conc= (DAPI_Nr_Set*calc.factor)/Volume,
         FISH_conc = (DAPI_Nr_SubSet*calc.factor)/Volume)


#define standard error function
se <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}

#Calculate cell abundances per sample
counts_all <-  raw.counts.SH %>% 
  group_by(Sample, Domain) %>% 
  summarise(DAPI.conc.mn = mean(DAPI_conc),
            DAPI.conc.md =  median(DAPI_conc), 
            DAPI.conc.sd = sd(DAPI_conc),
            FISH.conc.mn = mean(FISH_conc),
            FISH.conc.md =  median(FISH_conc), 
            FISH.conc.sd = sd(FISH_conc),
            n = length(DAPI_conc))%>%
  ungroup %>%
  as.data.frame()
