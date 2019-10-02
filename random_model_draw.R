library(surveillance)
library(hhh4addon)
library(rlist)
library(MASS)
library(ggplot2)

setwd("C:/Users/phpuenig/Dropbox/VL/Monthly prediction/surveillance forecasting/model selection - exploratory/random model draw")

#set.seed(101)
SCORES <- c("logs", "rps", "dss", "ses")

# Run all models on a subset in order to later compare with up to 12th order lagged model.
subset<-13:72

################################################################################
# Generate random models

randmods<-random.model(stsobj, 30)
save(randmods,file="randmods30_wcov.Rdata")

################################################################################
# Assessment

load("randmods30_wcov.RData")
randmods.list<-randmods[[1]]
OSAfirst<-randmods[[2]]
OSAroll<-randmods[[3]]

# Update models to just the training period
randmods.train<-lapply(randmods.list,update,subset=13:48)
save(randmods.train,file="randmods30_train.Rdata")

# Calculate n-step-ahead predictions (takes a long time - not used in paper)
# rand_3ahd<-lapply(randmods,FUN=stepaheadN,48,3)
# save(rand_3ahd, file="randmods30_3ahead.RData")
# rand_4ahd<-lapply(randmods,FUN=stepaheadN,48,4)
# save(rand_4ahd, file="randmods30_4ahead.RData")

scores.rand.first <- lapply(OSAfirst, scores, which = SCORES, individual = T)
scores.rand.roll <- lapply(OSAroll, scores, which = SCORES, individual = T)
# scores.rand.3ahd <- lapply(rand_3ahd, scores, which = SCORES, individual = T)
# scores.rand.4ahd <- lapply(rand_4ahd, scores, which = SCORES, individual = T)

save(scores.rand.first,file="randmods30_scorefirst.Rdata")
save(scores.rand.roll,file="randmods30_scoreroll.Rdata")

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
# scores_3ahd<-array(dim=c(length(scores.rand.3ahd),dim(scores.rand.3ahd[[1]])))
# for (i in 1:length(scores.rand.3ahd)){scores_first[i,,,]<-scores.rand.3ahd[[i]]}

scores_rolling_avg<-apply(scores_rolling,MARGIN=c(1,4),FUN="mean")
scores_first_avg<-apply(scores_first,MARGIN=c(1,4),FUN="mean")
# scores_3ahd_avg<-apply(scores_3ahd,MARGIN=c(1,4),FUN="mean")


# Calibration (based on RPS)
calibr.first <- lapply(OSAfirst, calibrationTest, which = "rps", individual = T)
calibr.rolling <- lapply(OSAroll, calibrationTest, which = "rps", individual = T)
# calibr.3ahd <- lapply(rand_3ahd, calibrationTest, which = "rps", individual = T)
# calib_all<-setNames(data.frame(array(dim=c(length(calibr.first),6))),c("calib.stat.first","calib.p.first","calib.stat.rolling","calib.p.rolling","calib.stat.3ahd","calib.p.3ahd"))
# for (i in 1:length(calibr.first)){calib_all[i,]<-c(calibr.first[[i]]$statistic,calibr.first[[i]]$p.value,calibr.rolling[[i]]$statistic,calibr.rolling[[i]]$p.value,calibr.3ahd[[i]]$statistic,calibr.3ahd[[i]]$p.value)}
calib_all<-setNames(data.frame(array(dim=c(length(calibr.first),4))),c("calib.stat.first","calib.p.first","calib.stat.roll","calib.p.roll"))
for (i in 1:length(calibr.first)){calib_all[i,]<-c(calibr.first[[i]]$statistic,calibr.first[[i]]$p.value,calibr.rolling[[i]]$statistic,calibr.rolling[[i]]$p.value)}

# Utility score
quants <- mapply(predquants,model=randmods.list,pred.output=OSAfirst,probs=list(0.1,0.25,0.45,0.55,0.75,0.9))
save(quants, file="quants_randmods.Rdata")
U <- lapply(quants, utility, cases[49:72,])
save(U, file="U_randmods.RData")
Uvec <- vector(length = 30)
for (i in 1:nmod){Uvec[i]<-U[[i]][[1]][2]} # extract score for 25-75% interval

# Combine into one results table
nmod<-nrow(scores_first_avg)
parmlim<-c(min(nparm),max(nparm))

result_table<-data.frame(Model=c(1:nmod),k=nparm,AIC_total=AIC_all,AIC_train=AIC_train)
result_table<-cbind(result_table,scores_first_avg,scores_rolling_avg,calib_all,Uvec)
names(result_table)[5:8]<-paste0(SCORES,".first")
names(result_table)[9:12]<-paste0(SCORES,".roll")
result_table$Miscal.first[result_table$calib.p.first>=0.1]<-"No evidence"
result_table$Miscal.roll[result_table$calib.p.roll>=0.1]<-"No evidence"
result_table$Miscal.first[result_table$calib.p.first<0.1]<-"Borderline"
result_table$Miscal.roll[result_table$calib.p.roll<0.1]<-"Borderline"
result_table$Miscal.first[result_table$calib.p.first<0.05]<-"Strong"
result_table$Miscal.roll[result_table$calib.p.roll<0.05]<-"Strong"

save(result_table,file="resulttable_randmods.Rdata")


################################################################################
# summary(randmods.list[[20]])
# for (i in 1:30){print(head(randmods.list[[i]]$control$ne$lag))}
# 
# ctl <- randmods.list[[1]]$control
# ctl$max_lag=1
# ctl$ar$lag=1
# ctl$ne$lag=7
# mod<-hhh4(stsobj,ctl)
# dim(ctl$end$offset)
# plot(mod,type="neweights")
