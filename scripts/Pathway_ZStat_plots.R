library(dplyr)
library(purrr)
library(readr)

source("scripts/helperFunctions.R")

AllDataAllTraitsCors <- read_csv("data/AllDataAllTraitsCors.csv")
uniqueStringDB <- read_csv("data/uniqueStringDB.csv")
pathwaySubset<-subset(AllDataAllTraitsCors, AllDataAllTraitsCors$Symbol %in% uniqueStringDB$x)


titles<- list()
for(n in 3:ncol(AllDataAllTraitsCors)){
  titles[n-2]<- paste("Pathway Enrichment Z Statistic of ", colnames(AllDataAllTraitsCors)[n], sep = ": ")
}

pdf(file = "plots/ZScore_surveyplot_phenotypes.pdf")
pmap(list(pathwaySubset[,3:ncol(AllDataAllTraitsCors)], AllDataAllTraitsCors[,3:ncol(AllDataAllTraitsCors)], titles), function(p, a,t) {
  plotZStatistics(pathCorrelat = p, allCorrelat = a, 
                  abs = T, histcols = "lightgrey",linecol = "midnightblue",lwd = 3,
                  histTitle = t)
})
dev.off()
