# Author:         Emily S Nightingale
# Institutions:   London Schoool of Hygiene and Tropical Medicine, London, UK
# Date Published: October 2019
################################################################################
# Consolidate fit and prediction metrics for all fitted models for comparison. 
# Calculate and evaluate three- and four-step-ahead predictions for the final 
# model.
################################################################################
################################################################################

probs <- list(0.1,0.25,0.45,0.55,0.75,0.9)

# set.seed(101)
# SCORES <- c("logs", "rps", "dss", "ses")
# 
# # Run all models on a subset in order to later compare with up to 4th order lagged model.
# subset<-5:72
# nblock <- 502

# Location of output from selection process
setwd("C:/Users/phpuenig/Dropbox/VL/Monthly Prediction/surveillance forecasting/Paper/results")
load("all_models.Rdata")
load("all_models_trainper.Rdata")
# load("osa_first.Rdata")
# load("osa_rolling.Rdata")
# load("scores_first.RData")
# load("scores_rolling.RData")
# load("coverage.RData")

# Location to save results of model evaluation
setwd("C:/Users/phpuenig/Dropbox/VL/Monthly Prediction/surveillance forecasting/Paper/results/evaluation/")

### ------------------ MODEL EVALUATION ------------------- ###

nmod <- length(model.list.all)

# Highlight selected model at each stage
selected <- c(1,6,23,33,42,52)
stage <- c(rep(1,12), rep(2,13), rep(3,10), rep(4,12), rep(5,11))

# length(model.list1) #12
# length(model.list2) #14
# length(model.list3) #11
# length(model.list4) #13 -> model 8/13 is final model
# length(model.list5) #12

# Model 42 is the final selected model
f <- 42

## Number of parameters
nparm <- sapply(model.list.all, FUN=function(x){k <- length(x$coefficients)
return(k)})
parmlim <- c(min(nparm), max(nparm))

## AIC for test and total period
AIC_all <- sapply(model.list.all, FUN=AIC)
AIC_train <- sapply(model.list.train, FUN=AIC)

# ---------------------------- #
#   Evaluate OSA Predictions   #
# ---------------------------- #

## Make OSA predictions based on rolling fit updates and on only the fit from the test period. 
osa.final.first <- lapply(model.list.all,oneStepAhead_hhh4lag,tp=test,type="first",which.start="current",keep.estimates=T)
osa.final.rolling <- lapply(model.list.all,oneStepAhead_hhh4lag,tp=test,type="rolling",which.start="current",keep.estimates=T)
save(osa.final.first, file="osa_first.Rdata")
save(osa.final.rolling, file="osa_rolling.Rdata")

## Proper scoring rules
models.scores.first <- lapply(osa.final.first, scores, which = SCORES, individual = T)
models.scores.rolling <- lapply(osa.final.rolling, scores, which = SCORES, individual = T)
save(models.scores.first, file="scores_first.RData")
save(models.scores.rolling, file="scores_rolling.RData")

# Rearrange into a 4D array of each model's four scores for each block-month
scores_first <- array(dim=c(length(models.scores.first), dim(models.scores.first[[1]])))
for (i in 1:length(models.scores.first)){scores_first[i,,,] <- models.scores.first[[i]]}
scores_rolling <- array(dim=c(length(models.scores.rolling), dim(models.scores.rolling[[1]])))
for (i in 1:length(models.scores.rolling)){scores_rolling[i,,,] <- models.scores.rolling[[i]]}

# Average each score over all blocks and months
scores_first_avg <- apply(scores_first, MARGIN=c(1,4), FUN="mean")
scores_rolling_avg <- apply(scores_rolling, MARGIN=c(1,4), FUN="mean")

## Calibration (based on RPS)
calibr.first <- lapply(osa.final.first, calibrationTest, which = "rps", individual = T)
calib_all <- setNames(data.frame(array(dim=c(length(calibr.first), 2))), c("calib.stat","calib.p"))
for (i in 1:length(calibr.first)){calib_all[i,] <- c(calibr.first[[i]]$statistic, calibr.first[[i]]$p.value)}

## Predicted quantiles
quants <- lapply(osa.final.first,predquants,probs=probs)
save(quants, file="quants.Rdata")

## Empirical coverage
coverage <- lapply(quants, covprob, cases[49:72,])
save(coverage, file="coverage.RData")

C1090 <- matrix(nrow = nmod, ncol=2)
C2575 <- matrix(nrow = nmod, ncol=2)
C4555 <- matrix(nrow = nmod, ncol=2)
for(i in 1:nmod){
  C1090[i,1] <- coverage[[i]][[1]][1]
  C2575[i,1] <- coverage[[i]][[1]][2]
  C4555[i,1] <- coverage[[i]][[1]][3]
  C1090[i,2] <- coverage[[i]][[2]][1]
  C2575[i,2] <- coverage[[i]][[2]][2]
  C4555[i,2] <- coverage[[i]][[2]][3]}
C_all <- data.frame(C1090=C1090[,1], C1090_qwd=C1090[,2], C2575=C2575[,1],
                    C2575_qwd=C2575[,2], C4555=C4555[,1], C4555_qwd=C4555[,2])

# MSE
models.mse.first <- sapply(osa.final.first, FUN = function(x){mean((x$pred - x$observed)^2)})
models.mse.rolling <- sapply(osa.final.rolling, FUN = function(x){mean((x$pred - x$observed)^2)})

summary(models.mse.first)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.498   1.603   1.655   1.851   1.883   3.875 
summary(models.mse.rolling)  
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.493   1.604   1.646   1.847   1.870   3.743 

models.mse.first[42]

1-(models.mse.first[42]/models.mse.first[1]) # 0.5363314
1-(models.mse.rolling[42]/models.mse.rolling[1]) # 0.5411541

1-(models.mse.first[23]/models.mse.first[1]) # 0.4876999

# ----------------------------- #
#    Construct result tables    #
# ----------------------------- #

result_table <- data.frame(Model=c(1:nmod), k=nparm, AIC_total=AIC_all,
                           AIC_train=AIC_train)
result_table <- cbind(result_table, scores_first_avg, #scores_rolling_avg,
                      calib_all, C_all)
names(result_table)[5:8] <- SCORES #paste0(SCORES, ".first")
#names(result_table)[9:12] <- paste0(SCORES, ".roll")
result_table$Miscal[result_table$calib.p >= 0.1] <- "No evidence"
result_table$Miscal[result_table$calib.p < 0.1] <- "Borderline"
result_table$Miscal[result_table$calib.p < 0.05] <- "Strong"
result_table$selected <- 0
result_table$selected[result_table$Model%in%selected] <- 1
result_table$stage <- stage
save(result_table, file="result_table.Rdata")
write.csv(result_table, file="result_table.csv", row.names = F)


# ------ Selected model table ----- #

result_selected <- result_table[selected,]
save(result_selected, file="result_table_selected.Rdata")
write.csv(result_selected, file="result_table_selected.csv")


# --------------------------------------------------- #
#   Evaluate 3/4SA Predictions for final model only   #
# --------------------------------------------------- #

# Predicted distributions for each block and time point in test period
model.final <- model.list.all[[f]]
pred3.final <- stepaheadN(model.final, start=48, type="first", n=3)
save(pred3.final, file="pred3_final.RData")
pred4.final <- stepaheadN(model.final, start=48, type="first", n=4)
save(pred4.final, file="pred4_final.RData")

# RPS of predictions
for(i in 1:length(pred3.final)){pred3.final[[i]][[2]]}
size <- 1/matrix(rep(pred3.final[[2]],nblock), nrow=length(pred3.final[[2]]))
scores_final_3ahead <- scores(x=cases[51:72,], mu=pred3.final[[1]], size=size, 
                              which = "rps", individual = T)[,,1]
for(i in 1:length(pred4.final)){pred4.final[[i]][[2]]}
size <- 1/matrix(rep(pred4.final[[2]],nblock), nrow=length(pred4.final[[2]]))
scores_final_4ahead <- scores(x=cases[52:72,], mu=pred4.final[[1]], size=size, 
                              which = "rps", individual = T)[,,1]

# Average RPS for 1/3/4SA over all blocks and over months 52:72 
# for comparison  
scores_mf_first <- colMeans(scores_first[f,-c(1:3),,2])
scores_mf_rolling <- colMeans(scores_rolling[f,-c(1:3),,2])
scores_3ahead_avg <- colMeans(scores_final_3ahead[-1,])
scores_4ahead_avg <- colMeans(scores_final_4ahead)
save(scores_mf_first, file="scores_mf_first.RData")
save(scores_mf_rolling, file="scores_mf_rolling.RData")
save(scores_3ahead_avg, file="scores_3ahead_final.RData")
save(scores_4ahead_avg, file="scores_4ahead_final.RData")

# hist(scores_3ahead_avg, breaks=30)
# hist(scores_4ahead_avg, breaks=30)

# Predicted quantiles
quants3_final <- predquants(pred3.final, probs)
quants4_final <- predquants(pred4.final, probs)
save(quants3_final, file="quants3_final.RData")
save(quants4_final, file="quants4_final.RData")

# Empirical coverage 
C.3ahd <- covprob(quants3_final, cases[51:72,])
C.4ahd <- covprob(quants4_final, cases[52:72,])
save(C.3ahd, file="C_3ahd.RData")
save(C.4ahd, file="C_4ahd.RData")


# Compare RPS of predictions 1/3/4 steps ahead (OSA calculated from rolling and "first"/training fit)
mean(models.scores.first[[42]][,,2]) # 0.4201441
mean(models.scores.rolling[[42]][,,2]) # 0.4197563
mean(scores_3ahead_avg) # 0.4413778
mean(scores_4ahead_avg) # 0.438052

permut.test2(models.scores.first[[42]][-c(1,2),,2],scores_final_3ahead)
# $`diffObs`
# [1] -0.02444261
# 
# $pVal.permut
# [1] 9.999e-05
# 
# $pVal.t
# [1] 1.722715e-19
permut.test2(models.scores.first[[42]][-c(1:3),,2],scores_final_4ahead)
# $`diffObs`
# [1] -0.0275056
# 
# $pVal.permut
# [1] 9.999e-05
# 
# $pVal.t
# [1] 8.024757e-17
permut.test2(scores_final_3ahead[-1,],scores_final_4ahead)
# $`diffObs`
# [1] -0.004445747
# 
# $pVal.permut
# [1] 0.02389761
# 
# $pVal.t
# [1] 0.02258403


