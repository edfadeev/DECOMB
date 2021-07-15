### microscopy calculations Jen ###
### percentages only ####

library(dplyr);packageVersion("dplyr")
library(tidyr);packageVersion("tidyr")
setwd("~/Jennifer/University/MSc/Master Thesis - MIDAS/Data/Microscopy")

#------- import raw counts -------#
library(readxl)
raw.counts <- read_excel(file.choose())
View(raw.counts)


###------rename column 1 ------###        
install.packages("tidyverse")
install.packages("tidyselect")
library("tidyverse")

raw.counts$sample_name <- str_replace(raw.counts$sample_name,"(C1|C2|C3|J1|J2|J3)(T4|T7)(ALT|EUB|GV|PSU)(HPG|RSG)",
    "\\1_\\2_\\3_\\4")
raw.counts$sample_name <- str_replace(raw.counts$sample_name,"(T0)(ALT|EUB|GV|PSU)(HPG|RSG)",
                                 "\\C0_\\1_\\2_\\3")
   
                              
###------ split colum 1 in serveral columns ------###   
raw.counts <- separate(raw.counts, col = 1, 
                   into = c("bottle", "timepoint", "probe", "trtment"), sep = "_", remove = FALSE,)

###------ microscope calculation factor ------###
filtration.surface.area <- 176714437.5 #um
counting.surface.area <- 14560 # lxw = 140*104 um
calc.factor <- 12136.9806 #filtration surface area/counting surface area um
volume <- 1 #mL

###------ calculate counts per FOV per mL ------###
raw.counts.conc <- raw.counts %>% 
  mutate(DAPI.conc = (DAPI*calc.factor)/volume,
        all.positive.conc = (DAPI_FITC_FISH*calc.factor)/volume,
         DAPI.FISH.conc = (DAPI_FISH*calc.factor)/volume,
         DAPI.FITC.conc = (DAPI_FITC*calc.factor)/volume,)

###------ define standard error function ------### NOT SURE WHAT I NEED IT FOR
se <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}

###------ Calculate cell abundances per sample (averages of FOV's) ------###
# averages of FOV's = conc (per mL) of cells grouped by sample bottle and probe
average_FOV <- raw.counts.conc %>% group_by(sample_nr, sample_name, bottle, timepoint, probe, trtment) %>% 
  summarise(DAPI.mn = mean(DAPI.conc),
            DAPI.md = median(DAPI.conc), # dropped median later on because if I wanna make Whisker Plots he does it anyway automatically
            DAPI.sd = sd(DAPI.conc),
            all.mn = mean(all.positive.conc),
            all.md = median(all.positive.conc),
            all.sd = sd(all.positive.conc),
            FITC.mn = mean(DAPI.FITC.conc),
            FITC.md = median(DAPI.FITC.conc),
            FITC.sd = sd(DAPI.FITC.conc),
            FISH.mn = mean(DAPI.FISH.conc),
            FISH.md = median(DAPI.FISH.conc),
            FISH.sd = sd(DAPI.FISH.conc)) %>%
  ungroup %>%
  as.data.frame()


# percentages of probe/trtm positive of all DAPI (already from averages) | ((100/DAPI.mn)*all.mn) = (all.mn/DAPI.mn)*100

av.percentages <- average_FOV %>% group_by(sample_nr, sample_name) %>%
  mutate(perc.all.pos = (100/DAPI.mn)*all.mn,
         perc.probe.pos = (100/DAPI.mn)*FISH.mn,
         perc.trt.pos = (100/DAPI.mn)*FITC.mn) %>%
  ungroup %>%
  as.data.frame()


###------ rename for combining bottles ------###

combined.bottles <- av.percentages

combined.bottles$sample_name <- str_replace(combined.bottles$sample_name, 
                              "(C1|C2|C3)(_)(T4|T7)(_)(ALT|EUB|GV|PSU)(_)(HPG|RSG)", "\\control\\2\\3\\4\\5\\6\\7")

combined.bottles$sample_name <- str_replace(combined.bottles$sample_name, 
              "(J1|J2|J3)(_)(T4|T7)(_)(ALT|EUB|GV|PSU)(_)(HPG|RSG)", "\\cteno\\2\\3\\4\\5\\6\\7")
  

###------ averages of bottles C1-C3, J1-J3, T0 per each probe ------###

install.packages("zoo")
library("zoo")

av.combined.bottles <- combined.bottles %>% group_by(sample_name, timepoint, probe, trtment) %>% 
  summarise(DAPI.mn = mean(DAPI.mn),
            DAPI.sd = sd(DAPI.sd),
            all.mn = mean(all.mn),
            all.sd = sd(all.sd),
            FITC.mn = mean(FITC.mn),
            FITC.sd = sd(FITC.sd),
            FISH.mn = mean(FISH.mn),
            FISH.sd = sd(FISH.sd),
            perc.all.mn = mean(perc.all.pos),
            perc.all.sd = sd(perc.all.pos),
            perc.probe.mn = mean(perc.probe.pos),
            perc.probe.sd = sd(perc.probe.pos),
            perc.trt.mn = mean(perc.trt.pos),
            perc.trt.sd = sd(perc.trt.pos)) %>%
  ungroup %>%
  as.data.frame()
# --> NA because some samples do not have sd since there is only one sample sometimes


# na.locf(av.combined.bottles, na.rm = FALSE, fromLast = TRUE) 

