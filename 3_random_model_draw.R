library(surveillance)
library(hhh4addon)
library(rlist)
library(MASS)
library(ggplot2)

setwd("C:/Users/phpuenig/Dropbox/VL/Monthly prediction/surveillance forecasting/Paper/Update Feb 2020/exploratory")

set.seed(101)
SCORES <- c("logs", "rps", "dss", "ses")

# Run all models on a subset in order to later compare with up to 12th order lagged model.
subset<-13:72

################################################################################
# Generate random models

randmods<-random.model(stsobj, 30)
save(randmods,file="randmods.Rdata")

################################################################################
# Assessment

#load("randmods.RData")
randmods.list<-randmods[[1]]
OSAfirst<-randmods[[2]]
OSAroll<-randmods[[3]]

# Update models to just the training period
randmods.train<-lapply(randmods.list,update,subset=13:48)
save(randmods.train,file="randmods_train.Rdata")

scores.rand.first <- lapply(OSAfirst, scores, which = SCORES, individual = T)
scores.rand.roll <- lapply(OSAroll, scores, which = SCORES, individual = T)

save(scores.rand.first,file="randmods_scorefirst.Rdata")
save(scores.rand.roll,file="randmods_scoreroll.Rdata")

# Number of parameters
nparm<-sapply(randmods.list,FUN=function(x){k<-length(x$coefficients)
return(k)})

# AIC for test and total period
AIC_all<-sapply(randmods.list,FUN=AIC)
AIC_train<-sapply(randmods.train,FUN=AIC)

# Proper scoring rules calculated from rolling and fixed fit 
scores_rolling<-array(dim=c(length(scores.rand.roll),dim(scores.rand.roll[[1]])))
for (i in 1:length(scores.rand.roll)){scores_rolling[i,,,]<-scores.rand.roll[[i]]}
scores_first<-array(dim=c(length(scores.rand.first),dim(scores.rand.first[[1]])))
for (i in 1:length(scores.rand.first)){scores_first[i,,,]<-scores.rand.first[[i]]}

scores_rolling_avg<-apply(scores_rolling,MARGIN=c(1,4),FUN="mean")
scores_first_avg<-apply(scores_first,MARGIN=c(1,4),FUN="mean")

# Calibration (based on RPS)
calibr.first <- lapply(OSAfirst, calibrationTest, which = "rps", individual = T)
calibr.rolling <- lapply(OSAroll, calibrationTest, which = "rps", individual = T)
calib_all<-setNames(data.frame(array(dim=c(length(calibr.first),2))),c("calib.stat.first","calib.p.first")) #,"calib.stat.roll","calib.p.roll"))
for (i in 1:length(calibr.first)){calib_all[i,]<-c(calibr.first[[i]]$statistic,calibr.first[[i]]$p.value)} #,calibr.rolling[[i]]$statistic,calibr.rolling[[i]]$p.value)}

# Combine into one results table
nmod<-nrow(scores_first_avg)
parmlim<-c(min(nparm),max(nparm))

results_random<-data.frame(Model=c(1:nmod),k=nparm,AIC_total=AIC_all,AIC_train=AIC_train)
results_random<-cbind(results_random,scores_first_avg,scores_rolling_avg,calib_all) 
names(results_random)[5:8]<-paste0(SCORES,".first")
names(results_random)[9:12]<-paste0(SCORES,".roll")
results_random$Miscal.first[results_random$calib.p.first>=0.1]<-"No evidence"
results_random$Miscal.first[results_random$calib.p.first<0.1]<-"Borderline"
results_random$Miscal.first[results_random$calib.p.first<0.05]<-"Strong"

save(results_random,file="results_random.Rdata")


