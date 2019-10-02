# Authors:        Tim Pollington/Emily S Nightingale
# Institutions:   London Schoool of Hygiene and Tropical Medicine, London, UK
# Date Published: XX July 2019
################################################################################
# This script contains premilinary analysis to explore the extent of spatial and 
# temporal correlation in the block-level time series. This is used to inform 
# the primary model selection.
################################################################################
################################################################################


nblock<-ncol(cases1)
nmth<-nrow(cases1)

#Investigate AR lags
critical.value <- 1.96/sqrt(nmth)
cols<-2+2*nmth

# For each neighbouring block-pair, append their time series in a single row.
# This list contains a row for every pair of first-order neighbours.
adj.list <- matrix(0,nrow = 1,ncol = cols) #134=2(i,j)+72(mo)+72(mo)
for (i in 1:nblock) {
  for (j in 1:nblock) {
    if (i < j) {
      if (neighbourhood(stsobj1)[i,j] == 1) {
        adj.list <- rbind(adj.list,c(i,j,rep.int(0,2*nmth)))
      }
    }
  }
}
adj.list <- adj.list[-1,] #remove redundant 1st row
for (i in 1:nrow(adj.list)) {
  adj.list[i,3:74]  <- t(stsobj1@observed[,adj.list[i,1]])
  adj.list[i,75:146] <- t(stsobj1@observed[,adj.list[i,2]])
}

# For each neighbouring pair, calculate the pacf of their bivariate time series
# and assess the significance of correlation at up to an 18 month lag. 
betw.pacf.result <- matrix(data = 0,nrow = 1, ncol = 19,dimnames = list(1,paste0('lag',0:18)))
for (i in 1:nrow(adj.list)) {
  tmp <- pacf(ts.union(ts(adj.list[i,3:(nmth+2)]),ts(adj.list[i,(nmth+3):cols])), plot = F, lag.max = 18)$acf[,,1][,2]
  print(i)
  lag <- 0
  while (lag != -1 & lag < 19) {
    if(!is.na(abs(tmp[lag+1])) & abs(tmp[lag+1]) > critical.value){
      betw.pacf.result[lag + 1] <- betw.pacf.result[lag + 1] + 1 # keep a tally of each block-block border connection
      lag <- lag + 1
    }
    else{
      lag <- -1
    }
  }
}
betw.pacf.result.percent <- betw.pacf.result/nrow(adj.list)*100

within.pacf.result <- matrix(data = 0,nrow = 1, ncol = 18,dimnames = list(1,paste0('lag',1:18)))
n_0<-0
for (i in 1:nblock) {
  tmp <- pacf(stsobj@observed[,i], plot = F, max.lag = 18)$acf
  if (sum(is.nan(tmp))>0){n_0<-n_0+1}
  lag <- 1
  while (lag != -1) {
    if(!is.nan(abs(tmp[lag])) & abs(tmp[lag]) > critical.value){
      within.pacf.result[lag] <- within.pacf.result[lag] + 1 # keep a tally of each district that passes for a continuous run of signif ACFs until one fails
      lag <- lag + 1
    }
    else{
      lag <- -1
    }
  }
}
within.pacf.result.percent <- within.pacf.result/nblock*100

round(betw.pacf.result.percent, 1)
# lag0 lag1 lag2 lag3 lag4 lag5 lag6 lag7 lag8 lag9 lag10 lag11 lag12 lag13 lag14 lag15 lag16 lag17 lag18
# 22.8  9.1  4.4  2.9    2  1.5  1.2  1.1  0.9  0.7   0.7   0.6   0.5   0.5   0.5   0.4   0.3   0.2     0
round(within.pacf.result.percent, 1)
# lag1 lag2 lag3 lag4 lag5 lag6 lag7 lag8 lag9 lag10 lag11 lag12 lag13 lag14 lag15 lag16 lag17 lag18
# 31.3  6.4  1.6  0.5    0    0    0    0    0     0     0     0     0     0     0     0     0     0

par(mfrow=c(2,1))
plot(c(0:18),betw.pacf.result.percent,type="l", xlab="Temporal Lag",ylab="Percentage",main="% with significant correlation between blocks")
plot(c(1:18),within.pacf.result.percent, type="l", xlab="Temporal Lag",ylab="Percentage",main="% with significant correlation within block")

