# Author:         Emily S Nightingale
# Institutions:   London Schoool of Hygiene and Tropical Medicine, London, UK
# Date Published: October 2019
################################################################################
# Starting from a basic model with only an endemic linear trend (equal across 
# all blocks), test possible additions to the model and select that which yields 
# the greatest reduction in RPS of one-step-ahead predictions. These predictions
# are made by estimating parameters from training period and then projecting 
# forward based on that fit. 
################################################################################
################################################################################

# Location for saving output
setwd("C:/Users/phpuenig/Dropbox/VL/Monthly Prediction/surveillance forecasting/Paper")

set.seed(101)
probs <- list(0.1,0.25,0.45,0.55,0.75,0.9)
SCORES <- c("logs", "rps", "dss", "ses")

# Run all models on a subset in order to later compare with up to 4th order lagged model.
# (NB. within comp.arlag function all are refit to 13:72 subset for comparison
# up to 12th. 
subset <- 5:72

# Define test period for evaluating predictive power
test <- c(48,71)

#------------------------------------------------------------------------------#

#----------------#
# M1: Basic model
#----------------#
basic.control <- list(end=list(f = ~1+t, offset=population(stsobj)), 
                      subset=subset,
                      family="NegBin1")
m1 <- hhh4(stsobj, control=basic.control)
# m1_nb<-update(m1,ne=list(f=~1, weights=neighbourhood(stsobj)==1))
# m1_nblag<-comp.nblag(m1_nb)

# Double check suitability of NB versus Poisson
mP <- update(m1, family="Poisson")

# # Compare temporal lags up to four months
# m1_ar <- update(m1, ar=list(~1))
# m1_arlag <- comp.arlag(m1,tp=test, type="first")
# save(m1_arlag, file = "./results/m1_arlag.Rdata")
# plot.ar(m1_arlag, file = "./lagcomp_final/m1_arcomp")
# 
# # Compare spatial lags
# m1_nb <- update(m2, ne=list(f = ~1, weights=neighbourhood(stsobj)==1))
# m1_nblag <- comp.nblag(m1_nb, test, type="first")
# save(m1_nblag, file = "./results/m1_nblag.Rdata")
# plot.nb(m1_nblag, file = "./lagcomp_final/m1_nbcomp")


#------------------------------------------------#
# M2: Test possible components with basic model
#------------------------------------------------#

m1a <- update(m1, end=list(f=addSeason2formula(~1+t, S=1, period=stsobj@freq), offset=population(stsobj)))
m1a2 <- update(m1, end=list(f=~1, offset=population(stsobj)))
m1a3 <- update(m1, end=list(f=~1+t+logpopdens, offset=1))
m1b <- update(m1, ar=list(f=~1, lag=1)) 
m1b2 <- update(m1, ar=list(f=addSeason2formula(~1, S=1, period=stsobj@freq), lag=1))
m1b3 <- update(m1, ar=list(f=addSeason2formula(~1+t, S=1, period=stsobj@freq), lag=1))
m1c <- update(m1, ne=list(f=~1, weights=W_powerlaw(maxlag=2)))
m1c2 <- update(m1, ne=list(f=~-1+logpopdens, weights=W_powerlaw(maxlag=2)))
m1c3 <- update(m1, ne=list(f=addSeason2formula(~1, S=1, period=stsobj@freq), weights=W_powerlaw(maxlag=2)))
m1c4 <- update(m1, ne=list(f=addSeason2formula(~1+t, S=1, period=stsobj@freq), weights=W_powerlaw(maxlag=2))) 
m1d <- update(m1, family=as.factor(df_wide$State))

# RE model
#m1e <- update(m1, end=list(f=~-1+ri(type="iid",corr="all"), offset=population(stsobj)))
# [13,] 1.029526 0.5425977 0.4651303 2.461902
#                         z 
#9.284141e-31 -1.153028e+01

model.list1 <- list(m1,m1a,m1a2,m1a3,m1b,m1b2,m1b3,m1c,m1c2,m1c3,m1c4,m1d)
save(model.list1, file="./results/modellist1.Rdata")

select1 <- modelassess(model.list1, test,type="first")
# logs       rps      dss      ses
# [1,] 1.208445 0.6569231 2.644396 3.263767
# [2,] 1.203510 0.6540003 2.542672 3.249254
# [3,] 1.226083 0.6978473 2.373225 3.873827
# [4,] 1.209090 0.6617311 2.403268 3.250645
# [5,] 1.039168 0.4950514 1.034439 1.883609
# [6,] 1.037969 0.4933758 1.031731 1.864375 ***
# [7,] 1.038272 0.4960356 1.000509 1.883507
# [8,] 1.023068 0.5155403 1.344840 2.307995
# [9,] 1.023272 0.5156969 1.350205 2.309949
# [10,] 1.022030 0.5146495 1.346082 2.293508
# [11,] 1.021963 0.5156422 1.273479 2.301376
# [12,] 1.205094 0.6591315 2.526938 3.266119
# z 
# 4.310509e-75 1.833550e+01 
# z 
# 3.297320e-57 1.594087e+01 
# z 
# 1.386028e-57 -1.599494e+01 
# z 
# 3.543340e-59 1.622169e+01 
# z 
# 0.1364016 1.4893256 
# z 
# 0.1444309 1.4594880 
# z 
# 0.0004326095 3.5193475732 
# z 
# 0.002699159 -3.000071852 
# z 
# 0.00151371 -3.17204229 
# z 
# 0.0009712016 -3.2987390026 
# z 
# 0.2009203 -1.2789340 
# z 
# 6.382924e-83 1.929108e+01 

# scores1 <- score_tidy(select1)

# Plot PIT histograms for all potential models
png(filename = "PIT_select1.png", height=700, width=1000)
par(mfrow=c(3,5))
for (m in model.list1){
  osa <- oneStepAhead_hhh4lag(m, tp=test, type="first", which.start = "current")
  pit(osa,J=50)
}
dev.off()

# Define best predicting model as M2
m2 <- m1b2
summary(m2, idx2Exp=T)
# Coefficients:
# Estimate   Std. Error
# exp(ar.1)                   6.685e-01  8.699e-03 
# exp(ar.sin(2 * pi * t/12))  1.161e+00  2.002e-02 
# exp(ar.cos(2 * pi * t/12))  1.028e+00  1.831e-02 
# exp(end.1)                  2.610e+02  7.461e+00 
# exp(end.t)                  9.908e-01  6.704e-04 
# overdisp                    7.408e-01  1.726e-02 
# 
# Log-likelihood:   -41009.4 
# AIC:              82030.8 
# BIC:              82081.43 
# 
# Number of units:        502 
# Number of time points:  68 

# Permutation test of improvement in scores (focus on RPS)
permut.test(m1,m2,test, type="first")
# logs          rps           dss          ses         
# diffObs     0.1704754     0.1635473     1.612665     1.399392    
# pVal.permut 9.999e-05     9.999e-05     9.999e-05    9.999e-05   
# pVal.t      1.749376e-121 1.524884e-109 6.818378e-16 7.519585e-38


#------#
#  M3
#------#

# # Test temporal and spatial lags in updated model
# m2_arlag <- comp.arlag(m2, tp=test, type="first")
# save(m2_arlag, file="./results/m2_arlag.Rdata")
# plot.ar(m2_arlag, file="./lagcomp_final/m2_arcomp")
# 
# m2_nb <- update(m2, ne=list(f=~1, weights=neighbourhood(stsobj)==1))
# m2_nblag <- comp.nblag(m2_nb,test,type="first")
# save(m2_nblag, file="./results/m2_nblag.Rdata")
# plot.nb(m2_nblag, file="./lagcomp_final/m2_nbcomp")


m2a <- update(m2, end=list(f=addSeason2formula(~1, S=1, period=stsobj@freq)))
m2a2 <- update(m2, end=list(f = ~1, offset=population(stsobj)))
m2a3 <- update(m2, end=list(f = ~1+logpopdens, offset=1))
#m2a4 <- update(m2, end=list(f=~-1+ri(type="iid",corr="all"), offset=population(stsobj)))
#[5,] 0.9828792 0.4863920 0.2442499 1.908734
#                        z 
#3.813094e-37 -1.273429e+01

m2b <- update(m2, ar=list(f=addSeason2formula(~1+t, S=1, period=stsobj@freq)))
m2b2 <- update(m2, ar=list(f=addSeason2formula(~1+t, S=2, period=stsobj@freq)))
#m2b3 <- update(m2, ar=list(f=addSeason2formula(~-1+ri(type="iid",corr="all"), S=1, period=stsobj@freq)))
#[8,] 1.0357375 0.4862850 1.0091507 1.969483
#                      z 
#0.002846629 2.983832166 

m2c <- update(m2, ne=list(f=~1, weights=W_powerlaw(maxlag = 2)))
m2c2 <- update(m2, ne=list(f=addSeason2formula(~1, S=1, period=stsobj@freq), weights=W_powerlaw(maxlag = 2)))
m2c3 <- update(m2, ne=list(f=addSeason2formula(~1+t, S=1, period=stsobj@freq), weights=W_powerlaw(maxlag = 2)))
m2c4 <- update(m2, ne=list(f=~-1+logpopdens, weights=W_powerlaw(maxlag = 2)))

m2d <- update(m2, family=as.factor(df_wide$State))

c2e <- list(end = list(f = ~1+t, offset=population(stsobj)),
            ar=list(f=addSeason2formula(~1, S=1, period=stsobj@freq)),
            max_lag = 2,
            subset = subset,
            family = "NegBin1")
m2e <- profile_par_lag(stsobj,control=c2e)
c2f <- list(end = list(f = ~1+t, offset=population(stsobj)),
            ar=list(f=addSeason2formula(~1, S=1, period=stsobj@freq)),
            max_lag = 3,
            subset = subset,
            family = "NegBin1")
m2f <- profile_par_lag(stsobj,control=c2f)
c2g <- list(end = list(f = ~1+t, offset=population(stsobj)),
            ar=list(f=addSeason2formula(~1, S=1, period=stsobj@freq)),
            max_lag = 4,
            subset = subset,
            family = "NegBin1")
m2g <- profile_par_lag(stsobj,control=c2g)

model.list2 <- list(m2,m2a,m2a2,m2a3,m2b,m2b2,m2c,m2c2,m2c3,m2c4,m2d,m2e,m2f,m2g)
save(model.list2, file="./results/modellist2.Rdata")

select2 <- modelassess(model.list2, test, type="first")
# logs       rps       dss      ses
# [1,] 1.0379691 0.4933758 1.0317307 1.864375
# [2,] 1.0470322 0.5022538 0.9886852 1.911420
# [3,] 1.0463911 0.5018110 0.9809106 1.907204
# [4,] 1.0472418 0.4991506 0.9816050 1.881122
# [5,] 1.0382724 0.4960356 1.0005095 1.883507
# [6,] 1.0384075 0.4960851 1.0014541 1.883513
# [7,] 0.9636050 0.4576658 0.4452471 1.734608
# [8,] 0.9632727 0.4574278 0.4463472 1.733270
# [9,] 0.9632998 0.4575646 0.4471027 1.734213
# [10,] 0.9641162 0.4581613 0.4461834 1.739059
# [11,] 1.0379414 0.4932968 1.0329081 1.864929
# [12,] 0.9831597 0.4547311 0.7147635 1.672341 ***
# [13,] 0.9552730 0.4389159 0.5607740 1.615055
# [14,] 0.9388063 0.4277159 0.5459888 1.541338
# z 
# 0.1444309 1.4594880 
# z 
# 3.791389e-19 -8.942885e+00 
# z 
# 8.182675e-18 -8.597000e+00 
# z 
# 4.305667e-24 -1.012441e+01 
# z 
# 0.0004326095 3.5193475732 
# z 
# 0.000279131 3.633935219 
# z 
# 0.2026591 -1.2740123 
# z 
# 0.287748 -1.063075 
# z 
# 0.1134269 -1.5829760 
# z 
# 0.1379056 -1.4836355 
# z 
# 0.1536718 1.4266814 
# z 
# 0.2186326 1.2301723 
# z 
# 0.00594747 2.75066353 
# z 
# 3.416377e-07 5.098897e+00 

# scores2 <- score_tidy(select2)

png(filename = "PIT_select2.png", height=700, width=1000)
par(mfrow=c(3,5))
for (m in model.list2){
  osa <- oneStepAhead_hhh4lag(m, tp=test, type="first", which.start = "current")
  pit(osa,J=50)
}
dev.off()

m3 <- m2e
summary(m3, idx2Exp=T)
# Estimate   Std. Error
# exp(ar.1)                   7.967e-01  7.958e-03 
# exp(ar.sin(2 * pi * t/12))  1.179e+00  1.552e-02 
# exp(ar.cos(2 * pi * t/12))  1.080e+00  1.458e-02 
# exp(end.1)                  1.408e+02  5.423e+00 
# exp(end.t)                  9.932e-01  8.826e-04 
# overdisp                    4.709e-01  1.298e-02 
# 
# Distributed lags used (max_lag = 2). Weights: 0.54; 0.46
# Use distr_lag() to check the applied lag distribution and parameters.
# 
# Log-likelihood:   -38731.62 
# AIC:              77477.24 
# BIC:              77536.31 
# 
# Number of units:        502 
# Number of time points:  68 

permut.test(m2,m3,tp=test, type="first")
# logs         rps          dss          ses         
# diffObs     0.05480942   0.03864467   0.3169673    0.1920342   
# pVal.permut 9.999e-05    9.999e-05    9.999e-05    9.999e-05   
# pVal.t      9.751091e-44 1.108821e-54 5.637425e-10 2.874575e-12


#------#
#  M4
#------#

# ctl.m3nb <- list(end = list(f = ~1+t, offset=population(stsobj)),
#                  ar = list(f = addSeason2formula(~1, S=1, period=stsobj@freq)),
#                  ne=list(f = ~1, weights=neighbourhood(stsobj)==1),
#                  max_lag = 2,
#                  subset = subset,
#                  family = "NegBin1")
# m3_nb <- profile_par_lag(stsobj,control=ctl.m3nb)
# m3_nblag <- comp.nblag(m3_nb, tp=test, type="first")
# save(m3_nblag, file="./results/m3_nblag.Rdata")
# plot.nb(m3_nblag, file="./lagcomp_final/m3_nbcomp")
# 
# m3_arlag <- comp.arlag(m2, tp=test, type="first")
# save(m3_arlag, file="./results/m3_arlag.Rdata")
# plot.ar(m3_arlag, file="./lagcomp_final/m3_arcomp")

c3a <- m3$control
c3a$end <- list(f = addSeason2formula(~1, S=1, period=stsobj@freq))
c3a2 <- m3$control
c3a2$end <- list(f = ~1, offset=population(stsobj))
c3a3 <- m3$control
c3a3$end <- list(f = ~1+logpopdens, offset=1)
c3a4 <- m3$control
c3a4$end <- list(f = ~-1 + ri(type="iid",corr="all"), offset=population(stsobj))
#[1,] 0.9641210 0.4657732 0.1876271 1.764208
# z 
# 1.022487e-36 -1.265707e+01 

c3b <- m3$control
c3b$ar <- list(f=addSeason2formula(~1+t, S=1, period=stsobj@freq))
c3b2 <- m3$control
c3b2$ar <- list(f=addSeason2formula(~1+t, S=2, period=stsobj@freq))
c3b3 <- m3$control
c3b3$max_lag <- 3
c3b4 <- m3$control
c3b4$ar <- list(f=addSeason2formula(~-1+ri(type="iid",corr="all")+t, S=1, period=stsobj@freq))
#[2,] 0.9897456 0.4549736 0.6855733 1.653120
# z 
# 0.0002544396 3.6577494163 

c3c <- m3$control
c3c$ne <- list(f=~1, weights=neighbourhood(stsobj)==1)
c3c2 <- m3$control
c3c2$ne <- list(f=addSeason2formula(~1, S=1, period=stsobj@freq), weights=neighbourhood(stsobj)==1)
c3c3 <- m3$control
c3c3$ne <- list(f=addSeason2formula(~1+t, S=1, period=stsobj@freq), weights=W_powerlaw(maxlag=3))
c3c4 <- m3$control
c3c4$ne <- list(f=~-1+logpopdens, weights=W_powerlaw(maxlag=3))

c3d <- m3$control
c3d$family <- as.factor(df_wide$State)

m3a <- profile_par_lag(stsobj,c3a)
m3a2 <- profile_par_lag(stsobj,c3a2)
m3a3 <- profile_par_lag(stsobj,c3a3)
m3a4 <- profile_par_lag(stsobj,c3a4)

m3b <- profile_par_lag(stsobj,c3b)
m3b2 <- profile_par_lag(stsobj,c3b2)
m3b3 <- profile_par_lag(stsobj,c3b3)
m3b4 <- profile_par_lag(stsobj,c3b4)

m3c <- profile_par_lag(stsobj,c3c)
m3c2 <- profile_par_lag(stsobj,c3c2)
m3c3 <- profile_par_lag(stsobj,c3c3)
m3c4 <- profile_par_lag(stsobj,c3c4)
m3d <- profile_par_lag(stsobj,c3d)
# m3d2<-profile_par_lag(stsobj,c3d2)

model.list3<-list(m3,m3a,m3a2,m3a3,m3b,m3b2,m3b3,m3c,m3c2,m3c3,m3d)
save(model.list3, file="./results/modellist3.Rdata")

select3<-modelassess(model.list3, test, type="first")
# logs       rps       dss      ses
# [1,] 0.9831597 0.4547311 0.7147635 1.672341
# [2,] 0.9877129 0.4568302 0.6824240 1.677588
# [3,] 0.9857774 0.4575140 0.6738022 1.683615
# [4,] 0.9860230 0.4562113 0.6724634 1.677049
# [5,] 0.9827643 0.4550524 0.6935611 1.672581
# [6,] 0.9830838 0.4548923 0.6965776 1.666718
# [7,] 0.9552730 0.4389159 0.5607740 1.615055
# [8,] 0.9400932 0.4367946 0.3495299 1.605271 
# [9,] 0.9422698 0.4366803 0.4042345 1.602907 ***
# [10,] 0.9401506 0.4371778 0.3677159 1.607890
# [11,] 0.9833730 0.4546964 0.7176195 1.672962
# z 
# 0.2186326 1.2301723 
# z 
# 8.371617e-08 -5.358930e+00 
# z 
# 1.504362e-05 -4.328049e+00 
# z 
# 1.884316e-08 -5.622300e+00 
# z 
# 0.00229939 3.04856255 
# z 
# 0.001178032 3.244146238 
# z 
# 0.00594747 2.75066353 
# z 
# 0.3476048 0.9392454 
# z 
# 0.1269068 1.5264138 
# z 
# 0.6307009 -0.4807405 
# z 
# 0.2227663 1.2192054 

# scores3 <- score_tidy(select3)

png(filename="PIT_select3.png", height=700, width=1000)
par(mfrow=c(3,5))
for (m in model.list3){
  osa <- oneStepAhead_hhh4lag(m, tp=test, type="first", which.start = "current")
  pit(osa,J=50)
}
dev.off()

m4 <- m3c2
summary(m4,idx2Exp=T)
# Estimate   Std. Error
# exp(ar.1)                    0.584349   0.008506 
# exp(ar.sin(2 * pi * t/12))   1.128417   0.022617 
# exp(ar.cos(2 * pi * t/12))   1.027358   0.021930 
# exp(ne.1)                    0.052826   0.001184 
# exp(ne.sin(2 * pi * t/12))   1.164880   0.032851 
# exp(ne.cos(2 * pi * t/12))   1.287080   0.038856 
# exp(end.1)                  25.030967   2.993279 
# exp(end.t)                   1.001409   0.002502 
# overdisp                     0.383121   0.011245 
# 
# Distributed lags used (max_lag = 2). Weights: 0.55; 0.45
# Use distr_lag() to check the applied lag distribution and parameters.
# 
# Log-likelihood:   -37154.76 
# AIC:              74329.51 
# BIC:              74413.89 
# 
# Number of units:        502 
# Number of time points:  68 

permut.test(m3, m4, tp=test, type="first")
# logs         rps         dss          ses        
# diffObs     0.04088989   0.01805087  0.3105289    0.06943311 
# pVal.permut 9.999e-05    9.999e-05   9.999e-05    9.999e-05  
# pVal.t      4.485867e-41 3.80822e-33 3.351044e-09 2.10258e-05



#------#
#  M5
#------#

# m4_nblag <- comp.nblag(m4, tp=test, type="first")
# save(m4_nblag, file="./results/m4_nblag.Rdata")
# plot.nb(m4_nblag, file="./lagcomp_final/m4_nbcomp")
# 
# m4_arlag <- comp.arlag(m4, tp=test, type="first")
# save(m4_arlag, file="./results/m4_arlag.Rdata")
# plot.ar(m4_arlag, file="./lagcomp_final/m4_arcomp")

# RPS and AIC worsen for spatial lags > 1 therefore adjust all
# models in this stage to include only direct neighbours. 
# Also test removal of endemic trend.

c4a <- m4$control
c4a$end <- list(f=addSeason2formula(~1, S=1, period=stsobj@freq))
c4a2 <- m4$control
c4a2$end <- list(f=~1, offset=population(stsobj))
c4a3 <- m4$control
c4a3$end <- list(f = ~1+logpopdens, offset=1)

# c4a2<-m4$control
# c4a2$end<-list(f = ~1+t, offset=log(pop_t[1:obsmth,]))

# c4a4<-m4$control
# c4a4$end<-list(f = ~1+t, offset=E)

c4b <- m4$control
c4b$ar <- list(f=~1+t)
c4b2 <- m4$control
c4b2$ar <- list(f=addSeason2formula(~1, S=2, period=stsobj@freq))
c4b3 <- m4$control
c4b3$ar <- list(f=addSeason2formula(~1 + t, S=2, period=stsobj@freq))
c4b3$ne <- list(f=addSeason2formula(~1 + t, S=2, period=stsobj@freq), weights=neighbourhood(stsobj)==1)
c4b4 <- m4$control
c4b4$max_lag <- 4
c4b4$end <- list(f=~1, offset=population(stsobj))

c4c <- m4$control
c4c$ne <- list(f=~1+t, weights=neighbourhood(stsobj)==1)
c4c2 <- m4$control
c4c2$ne <- list(f=addSeason2formula(~-1+logpopdens, S=1, period=stsobj@freq), weights=neighbourhood(stsobj)==1)
c4c3 <- m4$control
c4c3$ne <- list(f=~1, weights=W_powerlaw(maxlag=3))

c4d <- m4$control
c4d$family <- as.factor(df_wide$State)
c4d2 <- m4$control
c4d2$family <- as.factor(df_wide$State)
c4d2$ne <- list(f=addSeason2formula(~-1+logpopdens+t, S=1, period=stsobj@freq), weights=neighbourhood(stsobj)==1)

m4a <- profile_par_lag(stsobj,c4a)
m4a2 <- profile_par_lag(stsobj,c4a2)
m4a3 <- profile_par_lag(stsobj,c4a3)

m4b <- profile_par_lag(stsobj,c4b)
m4b2 <- profile_par_lag(stsobj,c4b2)
m4b3 <- profile_par_lag(stsobj,c4b3)
m4b4 <- profile_par_lag(stsobj,c4b4)

m4c <- profile_par_lag(stsobj,c4c)
m4c2 <- profile_par_lag(stsobj,c4c2)
m4c3 <- profile_par_lag(stsobj,c4c3)

m4d <- profile_par_lag(stsobj,c4d)
m4d2 <- profile_par_lag(stsobj,c4d2)

model.list4 <- list(m4,m4a,m4a2,m4a3,m4b,m4b2,m4b3,m4b4,m4c,m4c2,m4c3,m4d,m4d2)
save(model.list4, file="./results/modellist4.Rdata")

select4 <- modelassess(model.list4, test, type="first")
# logs       rps       dss      ses
# [1,] 0.9422698 0.4366803 0.4042345 1.602907
# [2,] 0.9438882 0.4371305 0.4228717 1.604054
# [3,] 0.9422393 0.4366907 0.4019233 1.602918
# [4,] 0.9434967 0.4370103 0.4196181 1.603861
# [5,] 0.9425955 0.4389949 0.3962532 1.639704
# [6,] 0.9423163 0.4367971 0.4017847 1.599133
# [7,] 0.9433994 0.4407658 0.3898617 1.643676 
# [8,] 0.9147223 0.4201441 0.2345101 1.513178 ***
# [9,] 0.9425624 0.4372183 0.3981955 1.603862
# [10,] 0.9438728 0.4380014 0.4024035 1.611063
# [11,] 0.9340098 0.4339387 0.3305142 1.596528
# [12,] 0.9424021 0.4366801 0.4053367 1.602882
# [13,] 0.9440710 0.4381503 0.4045313 1.611868
# z 
# 0.1269068 1.5264138 
# z 
# 0.07819457 1.76125903 
# z 
# 0.1814269 1.3363744 
# z 
# 0.141713 1.469442 
# z 
# 0.001138448 3.253868891 
# z 
# 0.1199699 1.5549001 
# z 
# 0.5874113 0.5425911 
# z 
# 0.3313647 0.9713687 
# z 
# 0.5829715 0.5490494 
# z 
# 0.2091886 1.2557998 
# z 
# 0.4159225 0.8135156 
# z 
# 0.1225734 1.5440616 
# z 
# 0.3312366 0.9716260 

# scores4 <- score_tidy(select4)

png(filename="PIT_select4.png", height=700, width=1000)
(mfrow=c(4,5))
for (m in model.list4){
  osa <- oneStepAhead_hhh4lag(m,tp=test,type="first",which.start = "current")
  pit(osa,J=50)
}
dev.off()

m5 <- m4b4
summary(m5, idx2Exp=T)
# Coefficients:
#   Estimate   Std. Error
# exp(ar.1)                    0.701931   0.008715 
# exp(ar.sin(2 * pi * t/12))   1.192184   0.021526 
# exp(ar.cos(2 * pi * t/12))   1.026418   0.018223 
# exp(ne.1)                    0.034791   0.001157 
# exp(ne.sin(2 * pi * t/12))   1.490941   0.066871 
# exp(ne.cos(2 * pi * t/12))   1.376839   0.058187 
# exp(end.1)                  16.499375   1.281805 
# overdisp                     0.304396   0.009870 
# 
# Distributed lags used (max_lag = 4). Weights: 0.32; 0.27; 0.23; 0.19
# Use distr_lag() to check the applied lag distribution and parameters.
# 
# Log-likelihood:   -36145.06 
# AIC:              72308.11 
# BIC:              72384.05 
# 
# Number of units:        502 
# Number of time points:  68  

permut.test(m4,m5,tp=test, type="first")
# logs         rps          dss         ses         
# diffObs     0.02754751   0.01653612   0.1697244   0.0897299   
# pVal.permut 9.999e-05    9.999e-05    0.00039996  9.999e-05   
# pVal.t      5.836594e-26 1.195758e-22 0.002376891 4.107036e-06


#------#
#  M6
#------#

# m5_nblag <- comp.nblag(m5, tp=test, type="first")
# save(m5_nblag, file="./results/m5_nblag.Rdata")
# plot.nb(m5_nblag, file="./lagcomp_final/m5_nbcomp")
# 
# m5_arlag <- comp.arlag(m5, tp=test, type="first")
# save(m5_arlag, file="./results/m5_arlag.Rdata")
# plot.ar(m5_arlag, file="./lagcomp_final/m5_arcomp")

c5a <- m5$control
c5a$end <- list(f=addSeason2formula(~1, S=1, period=stsobj@freq))
c5a2 <- m5$control
c5a2$ne <- list(f=addSeason2formula(~1+t, S=1,period=stsobj@freq), weights=neighbourhood(stsobj)==1)
c5a2$ar <- list(f=addSeason2formula(~1+t, S=1,period=stsobj@freq))
c5a3 <- m5$control
c5a3$end <- list(f = ~1+logpopdens, offset=1)

c5b <- m5$control
c5b$ar <- list(f=~1+t)
c5b2<- m5$control
c5b2$ar <- list(f=addSeason2formula(~1, S=2, period=stsobj@freq))
c5b3 <- m5$control
c5b3$ar <- list(f=addSeason2formula(~1 + t, S=2, period=stsobj@freq))
c5b3$ne <- list(f=addSeason2formula(~1 + t, S=2, period=stsobj@freq), weights=neighbourhood(stsobj)==1)

c5c <- m5$control
c5c$ne <- list(f=addSeason2formula(~1+t, S=1,period=stsobj@freq), weights=neighbourhood(stsobj)==1)
c5c2 <- m5$control
c5c2$ne <- list(f=addSeason2formula(~-1+logpopdens, S=1, period=stsobj@freq), weights=neighbourhood(stsobj)==1)
c5c3 <- m5$control
c5c3$ne <- list(f=~1, weights=neighbourhood(stsobj)==1)

c5d <- m5$control
c5d$family <- as.factor(df_wide$State)
c5d2 <- m5$control
c5d2$family <- as.factor(df_wide$State)
c5d2$ne <- list(f=addSeason2formula(~-1+logpopdens+t, S=1, period=stsobj@freq), weights=neighbourhood(stsobj)==1)

m5a <- profile_par_lag(stsobj,c5a)
m5a2 <- profile_par_lag(stsobj,c5a2)
m5a3 <- profile_par_lag(stsobj,c5a3)

m5b <- profile_par_lag(stsobj,c5b)
m5b2 <- profile_par_lag(stsobj,c5b2)
m5b3 <- profile_par_lag(stsobj,c5b3)

m5c <- profile_par_lag(stsobj,c5c)
m5c2 <- profile_par_lag(stsobj,c5c2)
m5c3 <- profile_par_lag(stsobj,c5c3)

m5d <- profile_par_lag(stsobj,c5d)
m5d2 <- profile_par_lag(stsobj,c5d2)


model.list5 <- list(m5,m5a,m5a2,m5a3,m5b,m5b2,m5b3,m5c,m5c2,m5c3,m5d,m5d2)
select5 <- modelassess(model.list5, test, type="first")
# logs       rps       dss      ses
# [1,] 0.9147223 0.4201441 0.2345101 1.513178
# [2,] 0.9155772 0.4203650 0.2414208 1.513562
# [3,] 0.9168291 0.4241301 0.2148901 1.558883
# [4,] 0.9152828 0.4202791 0.2436774 1.513418
# [5,] 0.9178918 0.4242995 0.2589604 1.562598
# [6,] 0.9150189 0.4194217 0.2352305 1.497541 ***
# [7,] 0.9163155 0.4231419 0.1973734 1.543941
# [8,] 0.9147330 0.4203081 0.2260580 1.512697
# [9,] 0.9159693 0.4208667 0.2341504 1.516199
# [10,] 0.9151533 0.4200215 0.2290569 1.512473 **
# [11,] 0.9147684 0.4201513 0.2348105 1.513194
# [12,] 0.9160168 0.4209004 0.2330223 1.516180
# z 
# 0.3313647 0.9713687 
# z 
# 0.2832815 1.0729770 
# z 
# 0.5833546 -0.5484913 
# z 
# 0.4218879 0.8031503 
# z 
# 1.582157e-06 4.800569e+00 
# z 
# 0.183114 1.331228 ***
# z 
# 0.8006071 0.2525615 
# z 
# 0.6088091 -0.5117741 
# z 
# 0.4076939 0.8279587 
# z 
# 0.2343799 1.1891518 **
# z 
# 0.3273144 0.9795376 
# z 
# 0.5475247 0.6014735 

# scores5 <- score_tidy(select5)

png(filename="PIT_select5.png", height=700, width=1000)
par(mfrow=c(3,4))
for (m in model.list5){
  osa <- oneStepAhead_hhh4lag(m,tp=test,type="first",which.start = "current")
  pit(osa,J=50)
}
dev.off()
save(model.list5, file="./results/modellist5.Rdata")

m6 <- m5b2
summary(m6, idx2Exp=T)
# Estimate   Std. Error
# exp(ar.1)                    0.700675   0.008748 
# exp(ar.sin(2 * pi * t/12))   1.343539   0.018082 
# exp(ar.cos(2 * pi * t/12))   1.131274   0.015126 
# exp(ne.1)                    0.035242   0.001098 
# exp(end.1)                  16.777846   1.292085 
# overdisp                     0.307730   0.009932 
# 
# Distributed lags used (max_lag = 4). Weights: 0.32; 0.27; 0.22; 0.19
# Use distr_lag() to check the applied lag distribution and parameters.
# 
# Log-likelihood:   -36199.19 
# AIC:              72412.39 
# BIC:              72471.45 
# 
# Number of units:        502 
# Number of time points:  68 

permut.test(m5,m6,tp=test, type="first")
# logs          rps          dss         ses         
# diffObs     -0.0004310566 0.0001226873 0.005453177 0.0007048315
# pVal.permut 0.6191381     0.8277172    0.8617138   0.9254075   
# pVal.t      0.6249124     0.8262726    0.8054116   0.9217112  




# Consolidate results -----------------------------------------------------


selected.models <- list(m1,m2,m3,m4,m5,m6)
save(selected.models, file="./results/selected_models.RData")

length(model.list1) #12
length(model.list2) #14
length(model.list3) #11
length(model.list4) #13 -> model 8/13 is final model
length(model.list5) #12 

model.list.all <- c(model.list1,model.list2[-1], model.list3[-1], model.list4[-1],model.list5[-1])
save(model.list.all, file="./results/all_models.Rdata")

# Model 42 is the final selected model.                 


# Refit all models to just the training period (first four years of data)
model.list.train <- lapply(model.list.all, update, subset=5:48)
save(model.list.train, file="./results/all_models_trainper.Rdata")

