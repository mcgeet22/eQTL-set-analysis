
for(phenot in 3:ncol(AllDataAllTraitsCors))
  fit<-plotZStatistics(phenoCorrelat = pathwaySubset[,phenot], allCorrelat = AllDataAllTraitsCors[,phenot], 
                       abs = T, histcols = "lightgrey",linecol = "darkgreen",lwd = 3,
                       histTitle = "Pathway Enrichment Z Statistic of Lung Phenotype")
text(x = 2.6, y = .55, labels = paste("pvalue", pnorm(q = 1.945061, mean = lungFit$estimate[1], sd = lungFit$estimate[2],lower.tail = F), sep = " : "), col = "darkgreen")




lungZ<-calculateZ(pathCorr = pathCor$GS.Lung, allCorr = AllDataAllTraitsCors$GS.Lung, abs = F)
leftVentZ<- calculateZ(pathCorr = pathCor$`GS.Right Ventricle`, allCorr = AllDataAllTraitsCors$`GS.Right Ventricle`, abs = T)
rightVZ<- calculateZ(pathCorr = pathCor$`GS.Right Ventricle`, allCorr = AllDataAllTraitsCors$`GS.Right Ventricle`, abs = T)

simulatedZ<- rep(NA, 10000)
set.seed(199)

for(n in 1:length(simulatedZ)){
  #randomly subseting the correlation matrix to get a random set of genes
  randomSubset<- sample_n(AllDataAllTraitsCors, size = nrow(pathCor),replace = F)
  simulatedZ[n]<- calculateZ(randomSubset$`GS.Right Ventricle`, allCorr = AllDataAllTraitsCors$`GS.Right Ventricle`, abs = F)
  #if(simulatedZ[n]>4){
  # print(randomSubset)
  #}
}

hist(simulatedZ, freq = F,xlim = c(-5,5))
fitdistrplus::fitdist(data = simulatedZ,"norm")
abline(v = rightVZ)
