library(dplyr)
library(purrr)
library(readr)

source("scripts/helperFunctions.R")

AllDataAllTraitsCors <- read_csv("data/AllDataAllTraitsCors.csv")
uniqueStringDB <- read_csv("data/uniqueStringDB.csv")
pathwaySubset<-subset(AllDataAllTraitsCors, AllDataAllTraitsCors$Symbol %in% uniqueStringDB$x)

pheno<- colnames(AllDataAllTraitsCors)[3:ncol(AllDataAllTraitsCors)]

zpTable<- pmap(list(pathwaySubset[,3:ncol(AllDataAllTraitsCors)], AllDataAllTraitsCors[,3:ncol(AllDataAllTraitsCors)], pheno), function(p,a,ph) {
    getZ_stats(pathCorrelat = p, allCorrelat = a,phenotype = ph, 
                  abs = T)
})

#I know this is ugly coding. I need to think of a cleaner way to clean this data
temp<- as.data.frame(zpTable)
temp<- t(temp)
colnames(temp)<- c("Phenotype", "Zscore", "P-value")
write.csv(temp, "output/summaryStats_pathwayeQTLenrich.csv",row.names = F)
