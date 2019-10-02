library(rlist)

# Run all models on a subset in order to later compare with up to 12th order lagged model.
subset<-13:obsmth

# Basic model: AR 1, END pop offset, NE 1, NegBin1
AIC_basic<-vector(length=12)
BIC_basic<-vector(length=12)
scores_basic<-matrix(nrow=12,ncol=4)
#calibp_lag<-matrix(nrow=12,ncol=4)
fits<-list()
preds<-list()
for (i in 1:12){
  ctl <- list(end = list(f = ~1, offset=population(stsobj)), 
                        ar = list(f = ~1, lag=1), 
                        ne = list(f = ~1, weights = neighbourhood(stsobj) == 1),
                        max_lag=i,
                        subset=subset,
                        family = "NegBin1")
  mod<-profile_par_lag(stsobj,control=ctl)
  osa<-oneStepAhead_hhh4lag(mod, tp=c(48,65), type = "rolling", which.start = "current")
  # pit(osa)
  scores_basic[i,]<-colMeans(scores(osa,which=SCORES))
  #calib<-calibrationTest(osa,which="rps", individual=T)
  #calibp_lag[i,]<-calib[["p.value"]]
  AIC_basic[i]<-AIC(mod)
  BIC_basic[i]<-BIC(mod)
  
  fits<-list.append(fits,mod)
  preds<-list.append(preds,osa)
}
save(fits, file="fits_basic_lag1_12.RData")
save(preds, file="preds_basic_lag1_12.RData")

load(file="fits_basic_lag1_12.RData")
load(file="preds_basic_lag1_12.RData")

calibp_lag<-vector(length=12)
for (i in 1:12){
  AIC_basic[i]<-AIC(fits[[i]])
  BIC_basic[i]<-BIC(fits[[i]])
  calib<-calibrationTest(preds[[i]],which="rps", individual=T)
  calibp_lag[i]<-calib[["p.value"]]
}
par(mfrow=c(1,1))
plot(c(1:12),BIC_basic,type="l", col="blue", ylab=" ", main="AIC (red) and BIC (blue) for different temporal lags", xlab="max lag")
lines(c(1:12),AIC_basic,type="l",col="red")
plot(c(1:12),calibp_lag, type="l", col="orange", main="Calibration test for different temporal lags", xlab="max lag", ylab="P-value")

par(mfrow=c(2,2))
for (i in 1:4){
  plot(c(1:12),scores_lag[,i], type="l", main=SCORES[i], xlab="max lag", ylab="")
}

# Model coeffients versus lags
coefs.basic<-matrix(nrow=12,ncol=4)
parmnames<-names(coefficients(fits[[1]]))
for (i in 1:12){
  coefs.basic[i,]<-coefficients(fits[[i]])
}
par(mfrow=c(2,2))
for (i in 1:4){
  plot(c(1:12),coefs.basic[,i], type="l", main=parmnames[i], xlab="max lag", ylab="")
}

par(mfrow=c(1,1))
for (i in 1:12){
  parlags[i]<-distr_lag(fits[[i]])$par_lag
}
library(boot)
plot(c(2:12),parlags[2:12],type="l", xlab="max lag", ylab="par lag", main="Par lag = logit(alpha), unnormalised weight of first lag")

hhhlag.results<-list(AIC_lag, BIC_lag, RPS_lag, calibp_lag)
save(lag.results, file="AR lag test - basic model.RData")


par(mfrow=c(2,2))
fit.avgerr<-vector(length=12)
osa.avgerr<-vector(length = 12)
for (i in 1:12){
  error.fit<-abs(fitted(fits[[i]])-cases[13:66,])
  error.osa<-abs(preds[[i]]$pred-preds[[i]]$observed)
  plot(c(13:66),rowSums(error.fit), type="l", col="blue", main=paste0("Lag ",i))
  lines(c(13:66),c(rep(NA,36),rowSums(error.osa)), col="red")
  abline(h=0,lty="dashed",col="grey")
  fit.avgerr[i]<-mean(error.fit)
  osa.avgerr[i]<-mean(error.osa)
}
par(mfrow=c(1,2))
plot(fit.avgerr,type="l", xlab="max lag", ylab="Absolute error", main="Fitted")
plot(osa.avgerr,type="l", xlab="max lag", ylab="", main="OSA")

par(mfrow=c(2,2))
for (i in c(1,4,8,12)){
  #plot(fits[[i]],total=T)
  error.fit<-abs(fitted(fits[[i]])-cases[13:66,])
  error.osa<-abs(preds[[i]]$pred-preds[[i]]$observed)
  plot(c(13:66),rowSums(error.fit), type="l", col="blue", main=paste0("Lag ",i))
  lines(c(13:66),c(rep(NA,36),rowSums(error.osa)), col="red")
  abline(h=0,lty="dashed",col="grey")
  fit.avgerr[i]<-mean(error.fit[])
  osa.avgerr[i]<-mean(error.osa[])
}
mean(rowSums(error.fit))
mean(rowSums(error.osa))



# Base model without neighbour component
AIC_noNE<-vector(length=12)
BIC_noNE<-vector(length=12)
scores_noNE<-matrix(nrow=12,ncol=4)
fits.noNE<-list()
preds.noNE<-list()
for (i in 1:12){
  ctl <- list(end = list(f = ~1, offset=population(stsobj)), 
              ar = list(f = ~1, lag=1), 
              max_lag=i,
              subset=subset,
              family = "NegBin1")
  mod<-profile_par_lag(stsobj,control=ctl)
  osa<-oneStepAhead_hhh4lag(mod, tp=c(48,65), type = "rolling", which.start = "current")
  # pit(osa)
  scores_noNE[i,]<-colMeans(scores(osa,which=SCORES))
  #calib<-calibrationTest(osa,which="rps", individual=T)
  #calibp_lag[i,]<-calib[["p.value"]]
  # AIC_lag[i]<-AIC(mod)
  # BIC_lag[i]<-BIC(mod)
  
  fits.noNE<-list.append(fits.noNE,mod)
  preds.noNE<-list.append(preds.noNE,osa)
}
save(fits.noNE, file="fits_basic_lag1_12.RData")
save(preds.noNE, file="preds_basic_lag1_12.RData")

for(i in c(1,4,8,12)){
  plot(fits.noNE[[i]],total=T)
  }
for (i in c(1,4,8,12)){
  plot(fits.noNE[[i]],total=T)
  error.fit<-fitted(fits.noNE[[i]])-cases[13:66,]
  error.osa<-preds.noNE[[i]]$pred-preds.noNE[[i]]$observed
  plot(c(13:66),rowSums(error.fit), type="l", col="blue", main=paste0("Lag ",i))
  lines(c(13:66),c(rep(NA,36),rowSums(error.osa)), col="red")
  abline(h=0,lty="dashed",col="grey")
  # fit.avgerr[i]<-mean(error.fit)
  # osa.avgerr[i]<-mean(error.osa)
}

# Scores versus lags
par(mfrow=c(2,2))
for (i in 1:4){
  plot(c(1:12),scores_noNE[,i], type="l", main=SCORES[i], xlab="max lag", ylab="")
}

# Model coeffients versus lags
coefs<-matrix(nrow=12,ncol=3)
parmnames<-names(coefficients(fits.noNE[[1]]))
for (i in 1:12){
  coefs[i,]<-coefficients(fits.noNE[[i]])
}
par(mfrow=c(1,3))
for (i in 1:3){
  plot(c(1:12),coefs[,i], type="l", main=parmnames[i], xlab="max lag", ylab="")
}

par(mfrow=c(1,1))
parlags.noNE<-vector(length=12)
for (i in 1:12){
  parlags.noNE[i]<-distr_lag(fits.noNE[[i]])$par_lag
}
plot(c(2:12),parlags.noNE[2:12],type="l", xlab="max lag", ylab="par lag", main="Par lag = logit(alpha), unnormalised weight of first lag")



# Model with seasonality
AIC_seas<-vector(length=12)
BIC_seas<-vector(length=12)
scores_seas<-matrix(nrow=12,ncol=4)
#calibp_lag<-matrix(nrow=12,ncol=4)
fits.seas<-list()
preds.seas<-list()
for (i in 1:12){
  ctl <- list(end = list(f = ~1, offset=population(stsobj)), 
              ar = list(f = addSeason2formula(~1+t,S=1,period=stsobj@freq)), 
              ne = list(f = addSeason2formula(~1+t,S=1,period=stsobj@freq), weights = neighbourhood(stsobj) == 1),
              max_lag=i,
              subset=subset,
              family = "NegBin1")
  mod<-profile_par_lag(stsobj,control=ctl)
  osa<-oneStepAhead_hhh4lag(mod, tp=c(48,65), type = "rolling", which.start = "current")
  # pit(osa)
  scores_seas[i,]<-colMeans(scores(osa,which=SCORES))
  #calib<-calibrationTest(osa,which="rps", individual=T)
  #calibp_lag[i,]<-calib[["p.value"]]
  AIC_seas[i]<-AIC(mod)
  BIC_seas[i]<-BIC(mod)
  
  fits.seas<-list.append(fits.seas,mod)
  preds.seas<-list.append(preds.seas,osa)
}
save(fits.seas, file="fits_seas_lag1_12.RData")
save(preds.seas, file="preds_seas_lag1_12.RData")
#load(file="fits_seas_lag1_12.RData")
#load(file="preds_seas_lag1_12.RData")

par(mfrow=c(2,2))
for (i in c(1,4,8,12)){
  #plot(fits[[i]],total=T)
  error.fit<-abs(fitted(fits.seas[[i]])-cases[13:66,])
  error.osa<-abs(preds.seas[[i]]$pred-preds.seas[[i]]$observed)
  plot(c(13:66),rowSums(error.fit), type="l", col="blue", main=paste0("Lag ",i))
  lines(c(49:66),rowSums(error.osa), col="red")
  abline(h=0,lty="dashed",col="grey")
  fit.avgerr[i]<-mean(error.fit)
  osa.avgerr[i]<-mean(error.osa)
}

par(mfrow=c(1,1))
plot(c(1:12),BIC_seas,type="l", col="blue", ylab=" ", main="AIC (red) and BIC (blue) - basic with AR/NE seasonality", xlab="max lag")
lines(c(1:12),AIC_seas,type="l",col="red")
#plot(c(1:12),calibp_lag, type="l", col="orange", main="Calibration test for different temporal lags", xlab="max lag", ylab="P-value")

par(mfrow=c(2,2))
for (i in 1:4){
  plot(c(1:12),scores_seas[,i], type="l", main=SCORES[i], xlab="max lag", ylab="")
}

# Model coeffients versus lags
coefs.seas<-matrix(nrow=12,ncol=10)
parmnames<-names(coefficients(fits.seas[[1]]))
for (i in 1:12){
  coefs.seas[i,]<-coefficients(fits.seas[[i]])
}
par(mfrow=c(2,5))
for (i in 1:10){
  plot(c(1:12),coefs.seas[,i], type="l", main=parmnames[i], xlab="max lag", ylab="")
}
par(mfrow=c(2,2))
for (i in 1:4){
  plot(c(1:12),scores_seas[,i], type="l", main=SCORES[i], xlab="max lag", ylab="")
}

par(mfrow=c(1,1))
parlags.seas<-vector(length=12)
for (i in 1:12){
  parlags.seas[i]<-distr_lag(fits.seas[[i]])$par_lag
}
plot(c(2:12),parlags.seas[2:12],type="l", xlab="max lag", ylab="par lag",main="Par lag = logit(alpha), unnormalised weight of first lag")

lag.results<-list(AIC_lag, BIC_lag, RPS_lag, calibp_lag)
save(lag.results, file="AR lag test - basic model.RData")

parlags.noNE<-vector(length=12)
for (i in 1:12){
  parlags.noNE[i]<-distr_lag(fits.noNE[[i]])$par_lag
}
plot(c(2:12),parlags.noNE[2:12],type="l", xlab="max lag", ylab="par lag", main="Par lag = logit(alpha), unnormalised weight of first lag")




models.rps <- lapply(preds, scores, which = "rps", individual = T)
plot(c(1:12),sapply(models.rps, mean))

par(mfrow=c(1,2))
plot(c(1:12),lapply(fits,BIC),type="l", col="blue", ylab=" ", main="AIC (red) and BIC (blue) for different temporal lags", xlab="max lag")
lines(c(1:12),lapply(fits,AIC),type="l",col="red")
plot(c(1:12),RPS_lag, type="l", col="forestgreen", main="RPS for different temporal lags", xlab="max lag", ylab="Ranked Probability Score (RPS)")
plot(c(1:12),calibp_lag, type="l", col="orange", main="Calibration test for different temporal lags", xlab="max lag", ylab="P-value")


fitsp<-list()
predsp<-list()
for (i in 2:12){ #error if 1 included
  ctlp <- list(end = list(f = ~1, offset=population(stsobj)), 
               ar = list(f = ~1, lag=1), 
               ne = list(f = ~1, weights = neighbourhood(stsobj) == 1),
               max_lag=i,
               funct_lag=poisson_lag,
               subset=subset,
               family = "NegBin1")
  ctlp$funct_lag=poisson_lag
  modp<-profile_par_lag(stsobj,control=ctlp)
  osap<-oneStepAhead_hhh4lag(modp, tp=c(48,65), type = "rolling", which.start = "current")
  fitsp<-list.append(fitsp,modp)
  predsp<-list.append(predsp,osap)
}
for (i in c(1,3,7,10)){
  error.fit<-abs(fitted(fitsp[[i]])-cases[13:66,])
  error.osa<-abs(predsp[[i]]$pred-predsp[[i]]$observed)
  plot(c(13:66),rowSums(error.fit), type="l", col="blue", main=paste0("Lag ",i+1))
  lines(c(13:66),c(rep(NA,36),rowSums(error.osa)), col="red")
  abline(h=0,lty="dashed",col="grey")
  # fit.avgerr[i]<-mean(error.fit)
  # osa.avgerr[i]<-mean(error.osa)
}


ctl <- list(end = list(f = ~1+t, offset=population(stsobj)), 
            ar = list(f = addSeason2formula(~1+t,S=2,period=stsobj@freq)), 
            ne = list(f=~1+t, weights = neighbourhood(stsobj) == 1),
            family = "NegBin1")
mod.seas2<-hhh4(stsobj,control=ctl)
summary(mod.seas2) # 70537.71 

ctl <- list(end = list(f = addSeason2formula(~1+t,S=2,period=stsobj@freq), offset=population(stsobj)), 
            ar = list(f = ~1+t), 
            ne = list(f=~1+t, weights = neighbourhood(stsobj) == 1),
            family = "NegBin1")
mod.seas2.1<-hhh4(stsobj,control=ctl)
summary(mod.seas2.1) # 70557.64

ctl <- list(end = list(f = ~1+t, offset=population(stsobj)), 
            ar = list(f = ~1+t), 
            ne = list(f=addSeason2formula(~1+t,S=2,period=stsobj@freq), weights = neighbourhood(stsobj) == 1),
            family = "NegBin1")
mod.seas2.2<-hhh4(stsobj,control=ctl)
summary(mod.seas2.2)
# AIC:              70430.49 
# BIC:              70522.57 

ctl <- list(end = list(f = addSeason2formula(~1+t,S=2,period=stsobj@freq), offset=population(stsobj)), 
            ar = list(f = ~1+t), 
            ne = list(f=addSeason2formula(~1+t,S=2,period=stsobj@freq), weights = neighbourhood(stsobj) == 1),
            family = "NegBin1")
mod.seas2.3<-hhh4(stsobj,control=ctl)
summary(mod.seas2.3) 
# AIC:              70415.57 
# BIC:              70541.13

ctl <- list(end = list(f = ~1+t, offset=population(stsobj)), 
            ar = list(f = addSeason2formula(~1+t,S=2,period=stsobj@freq)), 
            ne = list(f=addSeason2formula(~1+t,S=2,period=stsobj@freq), weights = neighbourhood(stsobj) == 1),
            family = "NegBin1")
mod.seas2.4<-hhh4(stsobj,control=ctl)
summary(mod.seas2.4) 
# AIC:              70393.7 
# BIC:              70519.26 **

# AR/NE seasonality is still best with S=2

ctl <- list(end = list(f = ~1+t, offset=population(stsobj)), 
            ar = list(f = addSeason2formula(~1+t,S=1,period=stsobj@freq)), 
            ne = list(f=addSeason2formula(~1+t,S=2,period=stsobj@freq), weights = neighbourhood(stsobj) == 1),
            family = "NegBin1")
mod.seas2.5<-hhh4(stsobj,control=ctl)
summary(mod.seas2.5) 
# AIC:              70401.84 
# BIC:              70510.66 **

ctl <- list(end = list(f = ~1+t, offset=population(stsobj)), 
            ar = list(f = addSeason2formula(~1+t,S=1,period=stsobj@freq)), 
            ne = list(f=addSeason2formula(~1+t,S=1,period=stsobj@freq), weights = neighbourhood(stsobj) == 1),
            family = "NegBin1")
mod.seas2.6<-hhh4(stsobj,control=ctl)
summary(mod.seas2.6) 
# AIC:              70531.26 
# BIC:              70623.34

ctl <- list(end = list(f = ~1+t, offset=population(stsobj)), 
            ar = list(f = addSeason2formula(~1+t,S=2,period=stsobj@freq)), 
            ne = list(f=addSeason2formula(~1+t,S=1,period=stsobj@freq), weights = neighbourhood(stsobj) == 1),
            family = "NegBin1")
mod.seas2.7<-hhh4(stsobj,control=ctl)
summary(mod.seas2.7) 
# AIC:              70475.46 
# BIC:              70584.29 

fits.seas2<-list()
preds.seas2<-list()
for (i in 1:12){
  ctl <- list(end = list(f = ~1, offset=population(stsobj)), 
              ar = list(f = addSeason2formula(~1+t,S=1,period=stsobj@freq)), 
              ne = list(f = addSeason2formula(~1+t,S=2,period=stsobj@freq), weights = neighbourhood(stsobj) == 1),
              max_lag=i,
              subset=subset,
              family = "NegBin1")
  mod<-profile_par_lag(stsobj,control=ctl)
  osa<-oneStepAhead_hhh4lag(mod, tp=c(48,65), type = "rolling", which.start = "current")
  # pit(osa)
  # scores_seas[i,]<-colMeans(scores(osa,which=SCORES))
  # #calib<-calibrationTest(osa,which="rps", individual=T)
  # #calibp_lag[i,]<-calib[["p.value"]]
  # AIC_seas[i]<-AIC(mod)
  # BIC_seas[i]<-BIC(mod)
  
  fits.seas2<-list.append(fits.seas2,mod)
  preds.seas2<-list.append(preds.seas2,osa)
}

for (i in 1:12){plot(fits.seas2[[i]],total=T)}
for (i in c(1,4,8,12)){
  error.fit<-abs(fitted(fits.seas2[[i]])-cases[13:66,])
  error.osa<-abs(preds.seas2[[i]]$pred-preds.seas2[[i]]$observed)
  plot(c(13:66),rowSums(error.fit), type="l", col="blue", main=paste0("Lag ",i))
  lines(c(13:66),c(rep(NA,36),rowSums(error.osa)), col="red")
  abline(h=0,lty="dashed",col="grey")
  # fit.avgerr[i]<-mean(error.fit)
  # osa.avgerr[i]<-mean(error.osa)
}

AIC_seas2<-vector(length=12)
BIC_seas2<-vector(length=12)
scores.seas2<-matrix(nrow=12,ncol=4)
for (i in 1:12){BIC_seas2[i]<-BIC(fits.seas2[[i]])
                AIC_seas2[i]<-AIC(fits.seas2[[i]])
                scores.seas2[i,]<-colMeans(scores(preds.seas2[[i]], which=SCORES))}
par(mfrow=c(1,1))
plot(c(1:12),BIC_seas,type="l", col="blue", ylab=" ", main="AIC (red) and BIC (blue) - basic with AR/NE seasonality", xlab="max lag")
lines(c(1:12),AIC_seas,type="l",col="red")
par(mfrow=c(2,2))
for (i in 1:4){
  plot(c(1:12),scores.seas2[,i], type="l", main=SCORES[i], xlab="max lag", ylab="")
}


### Final model ###


ctl.final <- list(end = list(f = ~1, offset=population(stsobj)), 
            ar = list(f = addSeason2formula(~1,S=2,period=stsobj@freq)), 
            ne = list(f = addSeason2formula(~1,S=2,period=stsobj@freq), weights = W_powerlaw(maxlag=2)),
            max_lag=4,
            subset=5:66,
            family = "NegBin1")
m.final<-profile_par_lag(stsobj,control=ctl.final)
 
# compare with basic model on same subset
ctl.basic <- list(end = list(f = ~1, offset=population(stsobj)), 
                  ar = list(f = ~1, lag=1), 
                  ne = list(f = ~1, weights = neighbourhood(stsobj)==1),
                  subset=5:66,
                  family = "NegBin1")
m.basic<-hhh4(stsobj,control=ctl.basic)


# Loop through all possible combinations
END.options<-list(a = ~1, 
                  b = ~1 + t, 
                  c = ~addSeason2formula(~1 + t, S=1, period=stsobj@freq), 
                  d = ~addSeason2formula(~1, S=1, period=stsobj@freq))
AR.options<-list(a = ~1, 
                 b = ~1 + t, 
                 c = ~addSeason2formula(~1 + t, S=1, period=stsobj@freq), 
                 d = ~addSeason2formula(~1, S=1, period=stsobj@freq))
NE.options<-list(a = ~1, 
                 b = ~1 + t, 
                 c = ~addSeason2formula(~1 + t, S=1, period=stsobj@freq), 
                 d = ~addSeason2formula(~1, S=1, period=stsobj@freq),
                 e = ~-1+log(pop), 
                 f = ~-1+log(pop) + t, 
                 g = ~addSeason2formula(~-1 + log(pop) + t, S=1, period=stsobj@freq), 
                 h = ~addSeason2formula(~-1 + log(pop), S=1, period=stsobj@freq))
NE.weight<-c(1,2)
#disp.options<-list(D1 = "NegBin1", DS = as.factor(df_wide$State))
disp.options<-list(D1 = "NegBin1", DS = "State")
disp.options<-list(D1 = "NegBin1", DS = "State")

subset=c(13:66)
osa.tp<-c(49:65)
model.index=0
output=list() 
for (lag in 1:3){
  for (lag.sp in 1:7){
    for (END.ctl in END.options){
      for (AR.ctl in AR.options){
        for (NE.ctl in NE.options){
          for (disp in disp.options){
            for (wt in NE.weight){
              
              model.index<-model.index+1
              print(model.index)
              # print(END.ctl)
              # print(AR.ctl)
              # print(NE.ctl)
              # print(disp)
              # print(wt)
              
              #model.id<-paste0("m-",)
              
              # if (lag == 1){
              #   if (wt == 1){
              #       control<-list(end = list(END.ctl, offset=population(stsobj))), 
              #                     ar = list(AR.ctl, lag=1), 
              #                     ne = list(NE.ctl, weights = neighbourhood(stsobj) == 1),
              #                     data = list(pop = popfrac),
              #                     subset=subset,
              #                     family = disp)
              #   print(control)
              #   fit<-hhh4(stsobj, control)
              #   osa<-oneStepAhead(fit, tp=osa.tp, which.start = "current", type="rolling", keep.estimates = T)
              #   mod.out<-list(fit=fit,osa=osa)}
              #   if (wt == 2){
              #   for (lag.sp in 2:7){
              #       control<-list(end = list(END.ctl, offset=population(stsobj))), 
              #                     ar = list(AR.ctl, lag=1), 
              #                     ne = list(NE.ctl, weights = W_powerlaw(max_lag=lag.sp)),
              #                     data = list(pop = popfrac),
              #                     subset=subset,
              #                     family = disp)
              #       print(control)
              #       fit<-hhh4(stsobj, control)
              #       osa<-oneStepAhead(fit, tp=osa.tp, which.start = "current", type="rolling", keep.estimates = T)
              #       mod.out<-list(fit=fit,osa=osa)
              #   }}}
              # 
              # 
              # if (lag > 1){
              #   if (wt == 1){
              #   control<-list(end = list(END.ctl, offset=population(stsobj))), 
              #                 ar = list(AR.ctl), 
              #                 ne = list(NE.ctl, weights = neighbourhood(stsobj) == 1),
              #                 data = list(pop = popfrac),
              #                 max_lag=lag,
              #                 subset=subset,
              #                 family = disp)
              #   print(control)
              #   fit<-hhh4_lag(stsobj, control)
              #   osa<-oneStepAhead_hhh4lag(fit, tp=osa.tp, which.start = "current", type="rolling", keep.estimates = T)
              #   mod.out<-list(fit=fit,osa=osa)}
              #   if (wt == 2){
              #   for (lag.sp in 2:7){
              #       control<-list(end = list(END.ctl, offset=population(stsobj))), 
              #                     ar = list(AR.ctl, lag=1), 
              #                     ne = list(NE.ctl, weights = W_powerlaw(max_lag=lag.sp)),
              #                     data = list(pop = popfrac),
              #                     subset=subset,
              #                     family = disp)
              #       print(control)
              #       fit<-hhh4(stsobj, control)
              #       osa<-oneStepAhead(fit, tp=osa.tp, which.start = "current", type="rolling", keep.estimates = T)
              #       mod.out<-list(fit=fit,osa=osa)
              #   }}}
              # 
              # output<-list.append(output,mod.out)
              # 
              }
            }
          }
        }
      }
    }
  }
}


# testctl <- list(end = list(f = ~1, offset=population(stsobj)), 
#                       ar = list(f = ~1, lag=1), 
#                       ne = list(f = ~-1+log(pop), weights = neighbourhood(stsobj) == 1),
#                       data=list(pop=popfrac),
#                       subset=subset,
#                       family = "NegBin1")
# fittest<-hhh4(stsobj,testctl)
# summary(fittest)
