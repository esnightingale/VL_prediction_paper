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

  
random.model<-function(sts.object, n, tp=c(48,71)){
  
  #First version had no covariate options. Second added pop and pop density:
  #cov<-sample(list(population(sts.object),logpopdens))
  
  seas1<-addSeason2formula(~1,S=1, period=sts.object@freq)
  seas1t<-addSeason2formula(~1+t,S=1, period=sts.object@freq)
  seas1cov<-addSeason2formula(~1+cov,S=1, period=sts.object@freq)
  seas1tcov<-addSeason2formula(~1+t+cov,S=1, period=sts.object@freq)
  #END.form.list<-list(~1,~1+t, seas1,seas1t)
  END.form.list<-list(~1,~1+t,~1+t+cov,~1+cov,seas1,seas1t,seas1cov,seas1tcov)
  END.offset.list<-list(NULL,population(sts.object))
  #AR.form.list<-list(~1,~1+t,seas1,seas1t)
  AR.form.list<-list(~1,~1+t,~1+t+cov,~1+cov,seas1,seas1t,seas1cov,seas1tcov)
  AR.lag.list<-c(1:12)
  #NE.form.list<-list(~1,~1+t,seas1,seas1t)
  NE.form.list<-list(~1,~1+t,~1+t+cov,~1+cov,seas1,seas1t,seas1cov,seas1tcov)
  NE.lag.list<-c(1:7)  #using powerlaw syntax with maxlag=1 same as neighbourhood(sts.object)==1?
  fam.list <- list("NegBin1",as.factor(df_wide$State))# ,"Poisson",as.factor(df_wide$District))
  
  n.opt<-sapply(list(END.form.list,END.offset.list,AR.form.list,AR.lag.list,NE.form.list,NE.lag.list,fam.list),length)
  
  models <- list()
  osa.first <- list()
  osa.roll <- list()
  
  while (length(models) < n){
  
  rand.opt <- sapply(n.opt,FUN=sample,size=1)
  
  END.form <- END.form.list[[rand.opt[1]]]
  END.offset <- END.offset.list[[rand.opt[2]]]
  AR.form <- AR.form.list[[rand.opt[3]]]
  AR.lag <- AR.lag.list[[rand.opt[4]]]
  NE.form <- NE.form.list[[rand.opt[5]]]
  NE.lag <- NE.lag.list[[rand.opt[6]]]
  fam <- fam.list[[rand.opt[7]]]
 

  # print(paste0("END = ",END.form,", offset = ",END.offset,", AR = ",AR.form,", lag = ",AR.lag,", NE = ",NE.form,", sp.lag = ",NE.lag,", family = ",fam))
  # print(fam)
  

  # print(paste0("list(end = list(f = ~",END.form,", offset=population(sts.object)), 
  #                 ar = list(f = ~",AR.form,"),
  #                 ne = list(f = ~",NE.form,", weights=W_powerlaw(maxlag = ",NE.lag,")),
  #                 subset=subset,
  #                 max_lag=",AR.lag,",
  #                 family = ",fam,")"))
  
  control <- list(end = list(f = END.form, offset=END.offset), 
                  ar = list(f = AR.form),
                  ne = list(f = NE.form, weights=W_powerlaw(maxlag = NE.lag)),
                  subset=subset,
                  max_lag=AR.lag,
                  family = fam)

  # assign(paste("model",i, sep = ""),profile_par_lag(sts.object,control=c1))
  model <- tryCatch(profile_par_lag(sts.object, control=control), error = function(e) e, warning = function(w) w) #, finally = function(){ models <- list.append(models, model) })
  osa1<-tryCatch(oneStepAhead_hhh4lag(model,tp=tp,type="first",which.start="current",keep.estimates=T), error = function(e) e, warning = function(w) w)
  osa2<-tryCatch(oneStepAhead_hhh4lag(model,tp=tp,type="rolling",which.start="current",keep.estimates=T), error = function(e) e, warning = function(w) w)
  # if(inherits(model, "error")) next 
  # if(inherits(model, "warning")) next   
  if(!(inherits(model, "error")|inherits(model, "warning")|inherits(osa1, "error")|inherits(osa1, "warning")|inherits(osa2, "error")|inherits(osa2, "warning"))){
    print(summary(model))
    models <- list.append(models, model) 
    osa.first<-list.append(osa.first,osa1)
    osa.roll<-list.append(osa.roll,osa2)
  }
  
  }
  return(list(models,osa.first,osa.roll))
}

# rand.mod.list<-random.model(stsobj, 15)
# save(rand.mod.list,file="randmodels_13to72.Rdata")
# rand.mod.list2<-random.model(stsobj, 10)
# save(rand.mod.list2,file="randmodels2_13to72.Rdata")
# osa.rand.first<-lapply(rand.mod.list,oneStepAhead_hhh4lag,tp=c(48,71),type="first",which.start="current",keep.estimates=T)
# osa.rand.roll<-lapply(rand.mod.list,oneStepAhead_hhh4lag,tp=c(48,71),type="rolling",which.start="current",keep.estimates=T)
# save(osa.rand.first, file="randmodels_osa_first.Rdata")
# save(osa.rand.roll, file="randmodels_osa_roll.Rdata")

# Output of first run of random model generator just included the models themselves, in a list. For the second run (10 models), calculated the OSA predictions alongside to 
# ensure necessary predictions converged for all models drawn. Therefore the second output contains a list of 3: the model list, OSA first and OSA rolling. 

#Five additional draws from version with covariate options
# rand.mod.list3<-random.model(stsobj, 5)
# save(rand.mod.list3,file="randmodels3_13to72.Rdata")



# Load 20 models + 5 with covariate options and join elements
load("randmods20.RData")
load("randmodels3_13to72.Rdata")
randmods25<-c(randmods20[[1]],rand.mod.list3[[1]])
OSAfirst25<-c(randmods20[[2]],rand.mod.list3[[2]])
OSAroll25<-c(randmods20[[3]],rand.mod.list3[[3]])


# Update models to just the training period
randmods25.train<-lapply(randmods25,update,subset=13:48)
save(randmods25.train,file="randmods25_train.Rdata")


# Calculate 3 ahead predictions
# rand25_3ahd<-lapply(randmods25,FUN=stepaheadN,48,3)
# save(rand25_3ahd, file="randmods_3ahead.RData")

scores.rand.first <- lapply(OSAfirst25, scores, which = SCORES, individual = T)
scores.rand.roll <- lapply(OSAroll25, scores, which = SCORES, individual = T)
scores.rand.3ahd <- lapply(rand25_3ahd, scores, which = SCORES, individual = T)
save(scores.rand.first,file="randmodels_scorefirst.Rdata")
save(scores.rand.roll,file="randmodels_scoreroll.Rdata")

#assess.rand<-modelassess(randmods25,c(48,71))
#save(assess.rand,file="./Results/randmodels_assess.Rdata")


# Number of parameters
nparm<-sapply(randmods25,FUN=function(x){k<-length(x$coefficients)
return(k)})

# AIC for test and total period
AIC_all<-sapply(randmods25,FUN=AIC)
AIC_train<-sapply(randmods25.train,FUN=AIC)


# Proper scoring rules calculated from rolling, fixed fit and 3 ahead predictions
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
calibr.first <- lapply(OSAfirst25, calibrationTest, which = "rps", individual = T)
calibr.rolling <- lapply(OSAroll25, calibrationTest, which = "rps", individual = T)
# calibr.3ahd <- lapply(rand25_3ahd, calibrationTest, which = "rps", individual = T)
# calib_all<-setNames(data.frame(array(dim=c(length(calibr.first),6))),c("calib.stat.first","calib.p.first","calib.stat.rolling","calib.p.rolling","calib.stat.3ahd","calib.p.3ahd"))
# for (i in 1:length(calibr.first)){calib_all[i,]<-c(calibr.first[[i]]$statistic,calibr.first[[i]]$p.value,calibr.rolling[[i]]$statistic,calibr.rolling[[i]]$p.value,calibr.3ahd[[i]]$statistic,calibr.3ahd[[i]]$p.value)}
calib_all<-setNames(data.frame(array(dim=c(length(calibr.first),4))),c("calib.stat.first","calib.p.first","calib.stat.rolling","calib.p.rolling"))
for (i in 1:length(calibr.first)){calib_all[i,]<-c(calibr.first[[i]]$statistic,calibr.first[[i]]$p.value,calibr.rolling[[i]]$statistic,calibr.rolling[[i]]$p.value)}


# Combine into one results table
nmod<-nrow(scores_first_avg)
parmlim<-c(min(nparm),max(nparm))

result_table25<-data.frame(Model=c(1:nmod),k=nparm,AIC_total=AIC_all,AIC_train=AIC_train)
result_table25<-cbind(result_table,scores_first_avg,scores_rolling_avg,calib_all) #,ncms_all)
names(result_table25)[5:8]<-paste0(SCORES,".first")
names(result_table25)[9:12]<-paste0(SCORES,".roll")
# names(result_table25)[13:14]<-c("calib.stat.roll","calib.p.roll")
# result_table25<-result_table25[,-15]
# result_table25$Model<-1:25
result_table25$Miscal.first[result_table25$calib.p.first>=0.1]<-"No evidence"
result_table25$Miscal.roll[result_table25$calib.p.roll>=0.1]<-"No evidence"
result_table25$Miscal.first[result_table25$calib.p.first<0.1]<-"Borderline"
result_table25$Miscal.roll[result_table25$calib.p.roll<0.1]<-"Borderline"
result_table25$Miscal.first[result_table25$calib.p.first<0.05]<-"Strong"
result_table25$Miscal.roll[result_table25$calib.p.roll<0.05]<-"Strong"

save(result_table25,file="result_table_randmods25.Rdata")


load("result_table_randmods.Rdata")
result_table<-result_table[,-(15:16)]
result_table$Miscal.rolling<-(result_table$calib.p.rolling<0.1)

result_table25<-rbind(result_table,result_table2[,-c(13,14,17)])
save(result_table25,file="result_table_randmods25.Rdata")

for (i in 1:15){print(head(rand.mod.list[[i]]$control$family))}

