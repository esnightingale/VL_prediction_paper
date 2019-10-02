# Author:         Emily S Nightingale
# Institutions:   London Schoool of Hygiene and Tropical Medicine, London, UK
# Date Published: XX July 2019
################################################################################
# Consolidate fit and prediction metrics for all fitted models for comparison
################################################################################
################################################################################

setwd("C:/Users/phpuenig/Dropbox/VL/Monthly Prediction/surveillance forecasting/final model selection/Results")

set.seed(101)
SCORES <- c("logs", "rps", "dss", "ses")

# Run all models on a subset in order to later compare with up to 4th order lagged model.
subset<-5:72
nblock <- 502

### ------------------ MODEL EVALUATION ------------------- ###

load("models_osafirst.Rdata")
load("models_trainper_osafirst.Rdata")
load("osa_first_osafirst.Rdata")
load("osa_rolling_osafirst.Rdata")
load("scores_first_osafirst.RData")
load("scores_rolling_osafirst.RData")
load("U_osafirst.RData")

## Consolidate model metrics ##

nmod <- length(model.list.all)

# Number of parameters
nparm <- sapply(model.list.all, FUN=function(x){k <- length(x$coefficients)
                                                return(k)})

# AIC for test and total period
AIC_all <- sapply(model.list.all, FUN=AIC)
AIC_train <- sapply(model.list.train, FUN=AIC)

# Proper scoring rules calculated from rolling and fixed fit
scores_rolling <- array(dim=c(length(models.scores.rolling), dim(models.scores.rolling[[1]])))
for (i in 1:length(models.scores.rolling)){scores_rolling[i,,,] <- models.scores.rolling[[i]]}
scores_first <- array(dim=c(length(models.scores.first), dim(models.scores.first[[1]])))
for (i in 1:length(models.scores.first)){scores_first[i,,,] <- models.scores.first[[i]]}

scores_rolling_avg <- apply(scores_rolling, MARGIN=c(1,4), FUN="mean")
scores_first_avg <- apply(scores_first, MARGIN=c(1,4), FUN="mean")

save(scores_rolling, file="scores_rolling.RData")
save(scores_first, file="scores_first.RData")

# Calibration (based on RPS)
calibr.first <- lapply(osa.final.first, calibrationTest, which = "rps", individual = T)
calibr.rolling <- lapply(osa.final.rolling, calibrationTest, which = "rps", individual = T)
calib_all <- setNames(data.frame(array(dim=c(length(calibr.first), 4))), c("calib.stat.first","calib.p.first","calib.stat.roll","calib.p.roll"))
for (i in 1:length(calibr.first)){calib_all[i,] <- c(calibr.first[[i]]$statistic, calibr.first[[i]]$p.value, calibr.rolling[[i]]$statistic, calibr.rolling[[i]]$p.value)}

U1090 <- matrix(nrow = nmod, ncol=2)
U2575 <- matrix(nrow = nmod, ncol=2)
U4555 <- matrix(nrow = nmod, ncol=2)
for(i in 1:nmod){
  U1090[i,1] <- U.osafirst[[i]][[1]][1]
  U2575[i,1] <- U.osafirst[[i]][[1]][2]
  U4555[i,1] <- U.osafirst[[i]][[1]][3]
  U1090[i,2] <- U.osafirst[[i]][[2]][1]
  U2575[i,2] <- U.osafirst[[i]][[2]][2]
  U4555[i,2] <- U.osafirst[[i]][[2]][3]}


U_all <- data.frame(U1090=U1090[,1], U1090_qwd=U1090[,2], U2575=U2575[,1],
                    U2575_qwd=U2575[,2], U4555=U4555[,1], U4555_qwd=U4555[,2])

# Highlight selected model at each stage
selected <- c(1,6,23,33,42,52)
stage <- c(rep(1,12), rep(2,13), rep(3,10), rep(4,12), rep(5,11))
length(stage)

# length(model.list1) #12
# length(model.list2) #14
# length(model.list3) #11
# length(model.list4) #13 -> model 8/13 is final model
# length(model.list5) #12

nmod <- nrow(scores_first_avg)
parmlim <- c(min(nparm), max(nparm))

result_table <- data.frame(Model=c(1:nmod), k=nparm, AIC_total=AIC_all,
                           AIC_train=AIC_train)
result_table <- cbind(result_table, scores_first_avg, scores_rolling_avg,
                      calib_all,U_all)
names(result_table)[5:8] <- paste0(SCORES, ".first")
names(result_table)[9:12] <- paste0(SCORES, ".roll")
result_table$Miscal.first[result_table$calib.p.first >= 0.1] <- "No evidence"
result_table$Miscal.roll[result_table$calib.p.roll >= 0.1] <- "No evidence"
result_table$Miscal.first[result_table$calib.p.first < 0.1] <- "Borderline"
result_table$Miscal.roll[result_table$calib.p.roll < 0.1] <- "Borderline"
result_table$Miscal.first[result_table$calib.p.first < 0.05] <- "Strong"
result_table$Miscal.roll[result_table$calib.p.roll < 0.05] <- "Strong"
result_table$selected <- 0
result_table$selected[result_table$Model%in%selected] <- 1
result_table$stage <- stage
save(result_table, file="result_table_osafirst.Rdata")


# Correlation between AIC and RPS 
# (Spearman more robust if not bivariate normal): 

library(energy)
mvnorm.etest(results[,c(4,6)],R=10000)
# Energy test of multivariate normality: estimated parameters
# 
# data:  x, sample size 58, dimension 2, replicates 10000
# E-statistic = 6.547, p-value < 2.2e-16

cor(results$AIC_train, results$rps.first, method ="spearman")
#[1] 0.9454109


# ------ Selected model table ----- #

result_selected <- result_table[selected,]
save(result_selected, file="result_selected_osafirst.Rdata")
write.csv(result_selected, file="result_selected_osafirst.csv")

par(mfrow=c(4,4))
lapply(osa.final.first, FUN=pit, J=50)
lapply(osa.final.rolling, FUN=pit, J=50)

# Model 42 is the final selected model
f <- 42

# Calculate and evaluate 3-/4-ahead predictions for final model
model.final <- model.list.all[[f]]
pred3.final <- stepaheadN(model.final, start=48, type="first", n=3)
save(pred3.final, file="pred3_final_osafirst.RData")
pred4.final <- stepaheadN(model.final, start=48, type="first", n=4)
save(pred4.final, file="pred4_final_osafirst.RData")

load("pred3_final_osafirst.RData")
load("pred4_final_osafirst.RData")
load("scores_rolling.RData")
load("scores_first.RData")

probs <- list(0.1,0.25,0.45,0.55,0.75,0.9)
quants3_final <- predquants(model.final, pred3.final, probs)
quants4_final <- predquants(model.final, pred4.final, probs)
save(quants3_final, file="quants3_final.RData")
save(quants4_final, file="quants4_final.RData")


for(i in 1:length(pred3.final)){pred3.final[[i]][[2]]}
size <- 1/matrix(rep(pred3.final[[2]],nblock), nrow=length(pred3.final[[2]]))
scores_final_3ahead <- scores(x=cases[51:72,], mu=pred3.final[[1]], size=size, 
                              which = "rps", individual = T)[,,1]
for(i in 1:length(pred4.final)){pred4.final[[i]][[2]]}
size <- 1/matrix(rep(pred4.final[[2]],nblock), nrow=length(pred4.final[[2]]))
scores_final_4ahead <- scores(x=cases[52:72,], mu=pred4.final[[1]], size=size, 
                              which = "rps", individual = T)[,,1]


scores_mf_first <- colMeans(scores_first[f,-c(1:3),,2])
scores_mf_rolling <- colMeans(scores_rolling[f,-c(1:3),,2])
scores_3ahead_avg <- colMeans(scores_final_3ahead)
scores_4ahead_avg <- colMeans(scores_final_4ahead)
save(scores_mf_first, file="scores_mf_first.RData")
save(scores_mf_rolling, file="scores_mf_rolling.RData")
save(scores_3ahead_avg, file="scores_3ahead_final.RData")
save(scores_4ahead_avg, file="scores_4ahead_final.RData")

hist(scores_3ahead_avg, breaks=30)
hist(scores_4ahead_avg, breaks=30)

## Utility for 3/4 ahead ##

U.3ahd <- utility(quants3_final, cases[51:72,])
U.4ahd <- utility(quants4_final, cases[52:72,])
save(U.3ahd, file="U_3ahd.RData")
save(U.4ahd, file="U_4ahd.RData")

