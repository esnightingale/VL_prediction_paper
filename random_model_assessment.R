setwd("C:/Users/phpuenig/Dropbox/VL/Monthly Prediction/surveillance forecasting/model selection 0119/Results")
library(hhh4addon)
library(rlist)
library(MASS)
library(ggplot2)

source("C:/Users/phpuenig/Dropbox/VL/Monthly Prediction/analysis data/data_run.R")

#set.seed(101)
SCORES <- c("logs", "rps", "dss", "ses")

# Run all models on a subset in order to later compare with up to 4th order lagged model.
subset<-13:72

# 
# END.form.list<-c("1","1+t","1+t+logpopdens","addSeason2formula(~1,S=1, period=sts.object@freq)","addSeason2formula(~1+t,S=1, period=sts.object@freq)","addSeason2formula(~1+t+logpopdens,S=1, period=sts.object@freq)")
# END.offset.list<-c("population(sts.object)","logpopdens")
# AR.form.list<-c("1","1+t","1+t+logpopdens","addSeason2formula(~1,S=1, period=sts.object@freq)","addSeason2formula(~1+t,S=1, period=sts.object@freq)","addSeason2formula(~1+t+logpopdens,S=1, period=sts.object@freq)")
# AR.lag.list<-c(1:12)
# NE.form.list<-c("1","1+t","-1+t+log(population(sts.object))","addSeason2formula(~1,S=1, period=sts.object@freq)","addSeason2formula(~1+t,S=1, period=sts.object@freq)","addSeason2formula(~-1+t+log(population(sts.object)),S=1, period=sts.object@freq)")
# NE.lag.list<-c(1:7)  #using powerlaw syntax with maxlag=1 same as neighbourhood(sts.object)==1?
# fam.list <- c("Poisson","NegBin1","as.factor(df_wide$State)","as.factor(df_wide$District)")


# END.form.list<-c("~1","~1+t","addSeason2formula(~1,S=1, period=sts.object@freq)","addSeason2formula(~1+t,S=1, period=sts.object@freq)")
# END.offset.list<-c("population(sts.object)","logpopdens")
# AR.form.list<-c("1","1+t","addSeason2formula(~1,S=1, period=sts.object@freq)","addSeason2formula(~1+t,S=1, period=sts.object@freq)")
# AR.lag.list<-c(1:12)
# NE.form.list<-c("1","1+t","addSeason2formula(~1,S=1, period=sts.object@freq)","addSeason2formula(~1+t,S=1, period=sts.object@freq)")
# NE.lag.list<-c(1:7)  #using powerlaw syntax with maxlag=1 same as neighbourhood(sts.object)==1?
# fam.list <- c("Poisson","NegBin1",as.factor(df_wide$State),as.factor(df_wide$District))
# 

  


rand.mod.list<-random.model(stsobj, 15)
save(rand.mod.list,file="randmodels_13to72.Rdata")
rand.mod.list2<-random.model(stsobj, 10)
save(rand.mod.list2,file="randmodels2_13to72.Rdata")

#Add covariate options
rand.mod.list3<-random.model(stsobj, 5)
save(rand.mod.list3,file="randmodels3_13to72.Rdata")
# Output of first run of random model generator just included the models themselves, in a list. For the second run (10 models), calculated the OSA predictions alongside to 
# ensure necessary predictions converged for all models drawn. Therefore the second output contains a list of 3: the model list, OSA first and OSA rolling. 

rand.mod.train<-lapply(rand.mod.list,update,subset=13:48)
save(rand.mod.train,file="randmodels_trainperiod.Rdata")

randmods2<-rand.mod.list2[[1]]
rand.mod.train2<-lapply(randmods2,update,subset=13:48)
save(rand.mod.train2,file="randmodels2_trainperiod.Rdata")

osa.rand.first<-lapply(rand.mod.list,oneStepAhead_hhh4lag,tp=c(48,71),type="first",which.start="current",keep.estimates=T)
osa.rand.roll<-lapply(rand.mod.list,oneStepAhead_hhh4lag,tp=c(48,71),type="rolling",which.start="current",keep.estimates=T)
save(osa.rand.first, file="randmodels_osa_first.Rdata")
save(osa.rand.roll, file="randmodels_osa_roll.Rdata")

osa.rand.first2<-rand.mod.list2[[2]]
osa.rand.roll2<-rand.mod.list2[[3]]
osa.rand.first3<-rand.mod.list3[[2]]
osa.rand.roll3<-rand.mod.list3[[3]]


scores.rand.first <- lapply(osa.rand.first, scores, which = SCORES, individual = T)
scores.rand.roll <- lapply(osa.rand.roll, scores, which = SCORES, individual = T)
save(scores.rand.first, file="randmodels_scores_first.RData")
save(scores.rand.roll,file="randmodels_scores_rolling.RData")

scores.rand.first2 <- lapply(osa.rand.first2, scores, which = SCORES, individual = T)
scores.rand.roll2 <- lapply(osa.rand.roll2, scores, which = SCORES, individual = T)
save(scores.rand.first2, file="randmodels2_scores_first.RData")
save(scores.rand.roll2,file="randmodels2_scores_rolling.RData")

#assess.rand<-modelassess(rand.mod.list,c(48,71))
#save(assess.rand,file="./Results/randmodels_assess.Rdata")

# Number of parameters
nparm<-sapply(rand.mod.list,FUN=function(x){k<-length(x$coefficients)
return(k)})
nparm<-sapply(randmods2,FUN=function(x){k<-length(x$coefficients)
return(k)})

# AIC for test and total period
AIC_all<-sapply(randmods2,FUN=AIC)
AIC_train<-sapply(randmods2,FUN=AIC)

# Proper scoring rules calculated from rolling and fixed fit
scores_rolling<-array(dim=c(length(scores.rand.roll2),dim(scores.rand.roll2[[1]])))
for (i in 1:length(scores.rand.roll2)){scores_rolling[i,,,]<-scores.rand.roll2[[i]]}
scores_first<-array(dim=c(length(scores.rand.first2),dim(scores.rand.first2[[1]])))
for (i in 1:length(scores.rand.first2)){scores_first[i,,,]<-scores.rand.first2[[i]]}

scores_rolling_avg<-apply(scores_rolling,MARGIN=c(1,4),FUN="mean")
scores_first_avg<-apply(scores_first,MARGIN=c(1,4),FUN="mean")

# Calibration (based on RPS)
calibr.first <- lapply(osa.rand.first2, calibrationTest, which = "rps", individual = T)
calibr.rolling <- lapply(osa.rand.roll2, calibrationTest, which = "rps", individual = T)
# calib_all<-setNames(data.frame(array(dim=c(length(calibr.rolling),2))),c("calib.stat.rolling","calib.p.rolling"))
# for (i in 1:length(calibr.rolling)){calib_all[i,]<-c(calibr.rolling[[i]]$statistic,calibr.rolling[[i]]$p.value)}
calib_all<-setNames(data.frame(array(dim=c(length(calibr.first),4))),c("calib.stat.first","calib.p.first","calib.stat.rolling","calib.p.rolling"))
for (i in 1:length(calibr.first)){calib_all[i,]<-c(calibr.first[[i]]$statistic,calibr.first[[i]]$p.value,calibr.rolling[[i]]$statistic,calibr.rolling[[i]]$p.value)}


nmod<-nrow(scores_first_avg)
parmlim<-c(min(nparm),max(nparm))

result_table2<-data.frame(Model=c(1:nmod),k=nparm,AIC_total=AIC_all,AIC_train=AIC_train)
result_table2<-cbind(result_table2,scores_first_avg,scores_rolling_avg,calib_all) #,ncms_all)
names(result_table2)[5:8]<-paste0(SCORES,".first")
names(result_table2)[9:12]<-paste0(SCORES,".rolling")
result_table2$Miscal.first<-(result_table2$calib.p.first<0.1)
result_table2$Miscal.rolling<-(result_table2$calib.p.rolling<0.1)
#save(result_table,file="result_table_randmods.Rdata")
save(result_table2,file="result_table_randmods2.Rdata")

load("result_table_randmods.Rdata")
result_table<-result_table[,-(15:16)]
result_table$Miscal.rolling<-(result_table$calib.p.rolling<0.1)

result_table25<-rbind(result_table,result_table2[,-c(13,14,17)])
save(result_table25,file="result_table_randmods25.Rdata")

for (i in 1:15){print(head(rand.mod.list[[i]]$control$family))}

