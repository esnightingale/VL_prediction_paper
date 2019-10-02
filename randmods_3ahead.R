library(surveillance)
library(hhh4addon)
library(dplyr)
library(reshape2)


# Calculate 3 ahead predictions
rand_3ahd<-lapply(randmods.list,FUN=stepaheadN,48,3)
save(rand_3ahd, file="randmods30_3ahead.RData")