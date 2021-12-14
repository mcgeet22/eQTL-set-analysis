

#FUNCTION TO CALCULATE THE Z STATISTIC
# import a single column of each, should be the same length
calculateZ <- function(pathCorrelat, allCorrelat, abs = F){
  # if you dont want to use the absolute value - 
  # by using the absolute value you are including the negative correlation and the positive correlations
  # so if there was one gene with a high positive correlation, and one with a low negative correlation, they won't cancel eachother out
  if(abs == F){
    subAvg<- mean(pathCorrelat[1])
    allAvg<- mean(allCorrelat)
  }else{
    subAvg<- mean(abs(pathCorrelat))
    allAvg<- mean(abs(allCorrelat))
  }
  #Xi - mu
  numerator<- subAvg- allAvg
  #sigma - sqrt(n)
  denom <- sd(allCorrelat) / sqrt(length(pathCorrelat))
  return(numerator / denom)  
}

#FUNCTION TO PLOT THE SIMULATED Z AND Z STATISTIC

plotZStatistics<- function(pathCorrelat, allCorrelat, abs = F, histcols = "steelblue2", linecol = "black", lwd = 3, histTitle = "Histogram of Z Statistics"){
    #calculate the Z for the pathway of interest
  pathwayZstat<- calculateZ(pathCorrelat = pathCorrelat, allCorrelat = allCorrelat, abs = abs)
  
  ### calculating the simulated Z values to get the distribution ###
  #generate the list
  simulatedZ<- rep(NA, 10000)
  
  #set seed
  set.seed(199)
  
  #iterate through n randomized samples
  for(n in 1:length(simulatedZ)){
    
    #randomly subseting the data to get a random set of genes
    randomSubset<- sample(x = allCorrelat, size = length(pathCorrelat), replace = F)

    # calculating Z for this randomized subset
    simulatedZ[n]<- calculateZ(pathCorrelat = randomSubset, allCorrelat = allCorrelat, abs = abs)
  }
  
  #plot the simulated data
  hist(simulatedZ, freq = F, main = histTitle, col = histcols)
  
  #visualize the pathway Zstat
  abline(v = pathwayZstat, col = linecol, lwd = lwd)
  
  #print out the Z statistic of our pathway genes
  fit<- fitdistrplus::fitdist(data = simulatedZ,"norm")
  pval<-pnorm(q = pathwayZstat, mean = fit$estimate[1], 
              sd = fit$estimate[2],lower.tail = F)
  text(x = pathwayZstat + .6, y = .5,
       labels = paste("pvalue", round(pval,digits = 4), sep = " : "), 
       col = linecol)

  
  #return the distribution so that we can calculate the p
  return(fit)
}

getZ_stats<- function(pathCorrelat, allCorrelat, phenotype, abs = F){
  #calculate the Z for the pathway of interest
  pathwayZstat<- calculateZ(pathCorrelat = pathCorrelat, allCorrelat = allCorrelat, abs = abs)
  
  ### calculating the simulated Z values to get the distribution ###
  #generate the list
  simulatedZ<- rep(NA, 10000)
  
  #set seed
  set.seed(199)
  
  #iterate through n randomized samples
  for(n in 1:length(simulatedZ)){
    
    #randomly subseting the data to get a random set of genes
    randomSubset<- sample(x = allCorrelat, size = length(pathCorrelat), replace = F)
    
    # calculating Z for this randomized subset
    simulatedZ[n]<- calculateZ(pathCorrelat = randomSubset, allCorrelat = allCorrelat, abs = abs)
  }

  fit<- fitdistrplus::fitdist(data = simulatedZ,"norm")
  
  pval<-pnorm(q = pathwayZstat, mean = fit$estimate[1], 
              sd = fit$estimate[2],lower.tail = F)

  statsList<- c(phenotype, pathwayZstat, pval)
  
  #return the distribution so that we can calculate the p
  return(statsList)
}

