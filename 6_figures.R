# Author:         Emily S Nightingale
# Institutions:   London Schoool of Hygiene and Tropical Medicine, London, UK
# Date Published: XX July 2019
################################################################################
# Produce figures in manuscript
################################################################################
################################################################################

library(viridis)

source("C:/Users/phpuenig/Documents/VL/VL_prediction_paper/1_data_run.R")
source("C:/Users/phpuenig/Documents/VL/VL_prediction_paper/plot_hhh4_amended.R")

set.seed(101)
SCORES <- c("logs", "rps", "dss", "ses")

path <- "C:/Users/phpuenig/Dropbox/VL/Monthly Prediction/surveillance forecasting/Paper/results/"
load(paste0(path,"all_models.Rdata"))
load(paste0(path,"./evaluation/osa_first.Rdata"))
load(paste0(path,"./evaluation/osa_rolling.Rdata"))
load(paste0(path,"./evaluation/scores_first.RData"))
load(paste0(path,"./evaluation/scores_rolling.RData"))
load(paste0(path,"./evaluation/quants.Rdata"))
load(paste0(path,"./evaluation/coverage.RData"))
load(paste0(path,"./evaluation/pred3_final.RData"))
load(paste0(path,"./evaluation/pred4_final.RData"))
load(paste0(path,"./evaluation/quants3_final.RData"))
load(paste0(path,"./evaluation/quants4_final.RData"))
load(paste0(path,"./evaluation/scores_mf_first.RData"))
load(paste0(path,"./evaluation/scores_mf_rolling.RData"))
load(paste0(path,"./evaluation/scores_3ahead_final.RData"))
load(paste0(path,"./evaluation/scores_4ahead_final.RData"))

load("C:/Users/phpuenig/Dropbox/VL/Monthly prediction/surveillance forecasting/exploratory/random model draw/result_table_randmods30.RData")
results_random <- result_table
load("C:/Users/phpuenig/Dropbox/VL/Monthly Prediction/surveillance forecasting/Paper/results/evaluation/result_table.Rdata")

setwd("C:/Users/phpuenig/Dropbox/VL/Monthly prediction/surveillance forecasting/Paper/figures")

#-----------------------DESCRIPTIVE-----------------------# 

timeall <- format(seq(as.Date("2013-1-1"), by="month", length.out = 73), "%Y-%m-%d")
yr <- seq(1,73,12)

wide2 <- dcast(input, OBJECTID+State+District+Block~interval_start, value.var = "count")
rownames(wide2) <- wide2$OBJECTID
pops2 <- dcast(input, OBJECTID+State+District+Block~interval_start, value.var = "pop")
rownames(pops2) <- pops2$OBJECTID
cases2 <- t(wide2[,-1:-4])
pops2 <- t(pops2[,-1:-4])
inc18 <- data.frame(id=wide2$OBJECTID, value=colSums(cases2[61:72,])*1e4/colMeans(pops2[61:72,]))

# Block-level incidence per 10,000 in 2018
#png(filename = "./png/Fig1.png", height = 400, width = 600)
tiff(filename = "./tif/Fig1.tif", height = 5, width = 7, units = "in", res = 300)
mapplot(shapefile=VL, data=inc18, legend_title="Cases per 10,000", value = value) +
  ggtitle("Block-level incidence per 10,000 in 2018", )
dev.off()

# Average monthly incidence per 10,000 in 2018
inc18mth <- data.frame(id=wide2$OBJECTID, value=colMeans(cases2[61:72,]*1e4/pops2[61:72,]))
#png(filename = "./png/Fig1a.png", height = 400, width = 600)
tiff(filename = "./tif/Fig1a.tif", height = 5, width = 7, units = "in", res = 300)
mapplot(shapefile=VL, data=inc18mth, legend_title="Cases per 10,000", value = value) #+ 
  #ggtitle("Block-level, average monthly incidence per 10,000 over 2018")
dev.off()

par(mfrow=c(1,1))
#png(filename="./png/Fig2.png", height=300, width=400)
tiff(filename="Fig2.tif",height=4,width=5, units = "in", res = 300)
barplot(rowSums(cases), xaxt="n", xlab="", ylab="No. reported cases", space=0)
axis(1, at=yr-1, labels=timeall[yr], las=2 , cex.axis=0.8)
dev.off()

blocktotal <- colSums(cases)
#png(filename="./png/Fig2b.png", height=450, width=600)
tiff(filename = "./tif/Fig2b.tif", height = 4.5, width = 6, units = "in", res = 300)
hist(blocktotal, breaks=30, xlab="Total cases Jan 2013-Dec 2018", 
     main="Block-level distribution of cases over whole time period (N = 502)")
dev.off()


#----------------- RANDOM MODEL ASSESSMENT ---------------------#

fitvfirst <- ggplot(results_random,aes(AIC_train,rps.first, col=k)) +
  theme(legend.position = c(0.6,0.2)) +
  labs(colour="No. parameters") +
  geom_point() +
 # geom_label_repel(aes(label = Model)) +
  xlab("AIC") + ylab("RPS") + 
  ggtitle("Training period fit") +
  scale_color_viridis_c(option = "C") 
firstvroll <- ggplot(results_random,aes(rps.first,rps.roll, col=k)) +
  theme(legend.position = "none") +
  geom_point() +
 # geom_label_repel(aes(label = Model)) +
  xlab("RPS (fixed)") + ylab("RPS (rolling)") + 
  ggtitle("Training period fit versus rolling updates") +
  scale_color_viridis_c(option = "C") 

#png(filename = "./png/Fig3.png", height=400, width=800)
tiff(filename = "./tif/Fig3.tif", height=5, width=10, units = "in", res=300)
plot_grid(fitvfirst+geom_point(cex=3), firstvroll+geom_point(cex=3), labels=c("A", "B"), ncol = 2)
dev.off()

#----------------- SYSTEMATIC SELECTION ------------------------#

result_table$highlight[result_table$selected==0] <- 16
result_table$highlight[result_table$selected==1] <- 5
result_table$highlight[c(1,58)] <- 16
result_table$highlight[42] <- 8
nmod <- nrow(result_table)

palette(plasma(6))
#png(filename = "./png/Fig4.png", width=600, height=500)
tiff(filename = "./tif/Fig4.tif", height = 6.5, width = 7.5, units = "in", res = 300)
par(mfrow=c(2,2), mar=c(4,4,3,4))

plot(c(1:nmod), result_table$rps, 
     main="RPS", ylab="", xlab="Model no.", cex=0.7,
     pch=result_table$highlight, col=as.factor(result_table$stage))
legend(x=30, y=0.65, legend=c("1","2","3","4 (final)"),
       title="Selected model",
       col=plasma(5)[1:4],
       pch=c(5,5,5,8),
       cex = 0.7,
       bty="n")
plot(c(1:nmod), result_table$AIC_train, 
     main="AIC (training set)", ylab="", xlab="Model no.", cex=0.7, 
     pch=result_table$highlight, col=as.factor(result_table$stage))
plot(c(1:nmod), (1-result_table$C2575), xlim=c(1,nmod), 
     main="Empirical coverage - 25-75%", ylab="", xlab="Model no.", 
     pch=result_table$highlight, cex=0.7, 
     col=as.factor(result_table$stage))
par(new = T)
plot(c(1:nmod), result_table$C2575_qwd, xlim=c(1,nmod), axes=F, xlab=NA, ylab=NA, 
     type="l", lty="dashed", col="darkgrey")
axis(side = 4)
mtext(side = 4, line = 2, 'Average interval width')
plot(c(1:nmod), (1-result_table$C1090), 
     main="Empirical coverage - 10-90%", ylab="", xlab="Model no.", cex=0.7,
     pch=result_table$highlight,col=as.factor(result_table$stage))
par(new = T)
plot(c(1:nmod), result_table$C1090_qwd, xlim=c(1,nmod), axes=F, xlab=NA, ylab=NA,
     type="l", lty="dashed", col="darkgrey")
axis(side = 4)
mtext(side = 4, line = 2, 'Average interval width')
dev.off()

# PIT histograms for selected models
selected <- c(1,6,23,33,42,52)
#png(filename ="./png/Fig5.png", height=400, width=600)
tiff(filename = "./tif/Fig5.tif", height = 6, width = 10, units = "in", res = 300)
par(mfrow=c(2,3))
i <- 1
for (m in selected){
  pit(osa.final.first[[m]], J=30, plot=list(ylim=c(0,2), main=paste0("Model no: ", selected[i])))
  i <- i+1
}
dev.off()

## ----------------------------- FINAL MODEL -------------------------------- ##

f <- 42  #index of final model
t.pred <- obsmth-48  #length of test period
nblock <- dim(cases)[2]
m.final <- model.list.all[[f]]
final.OSApred.first <- osa.final.first[[f]]
final.rps.first <- models.scores.first[[f]][,,2]
C1_mat <- coverage[[f]][[3]][,,1]
C2_mat <- coverage[[f]][[3]][,,2]
C3_mat <- coverage[[f]][[3]][,,3]

# SEASONAL WAVES
t <- 5:72
coefs <- coefficients(m.final)
log_ar <- coefs[1]+coefs[2]*sin(2*pi*t/12)+coefs[3]*cos(2*pi*t/12)
log_ne <- coefs[4]+coefs[5]*sin(2*pi*t/12)+coefs[6]*cos(2*pi*t/12)

#png(filename = "./png/Fig6.png", height=350, width=550)
tiff(filename = "./tif/Fig6.tif", height = 5, width = 8, units = "in", res = 300)
par(mfrow = c(2,1), mar = c(3, 4, 2, 2))
plot(t, exp(log_ar), type = "l", ylab = "AR component", xlab = "", xaxt = "n")
#axis(1, at=yr-1, labels=timeall[yr], las=2 , cex.axis=0.8)
par(mar = c(5, 4, 1, 2))
plot(t, exp(log_ne), type = "l", ylab = "NE component", xlab = "", xaxt = "n")
axis(1, at=yr-1, labels=timeall[yr], las=2 , cex.axis=0.8)
dev.off()

# FIT FOR HIGHEST INCIDENCE BLOCKS
df.rps <- data.frame(id=df_wide$OBJECTID, value=colMeans(final.rps.first))
avginc <- colMeans(inc)
highinc <- names(sort(avginc,decreasing = T))[1:4]
df_wide$Block[which(df_wide$OBJECTID%in%highinc)] 
#"GOPIKANDAR"   "KATHIKUND"    "BOARIJOR"     "SUNDARPAHARI"
blknames <- paste0(df_wide$Block[which(df_wide$OBJECTID%in%highinc)],
                 " (RPS = ",
                 round(df.rps$value[df.rps$id%in%highinc],2),")")
#png(filename = "./png/Fig7.png", height=600, width=800)
tiff(filename = "./tif/Fig7.tif", height = 7, width = 9, units = "in", res = 300)
plot_hhh42(m.final, units=highinc, names=blknames, ylab="No. reported cases",
           col=c("forestgreen","skyblue", "orange"))
dev.off()

#png(filename = "./png/Fig8.png", height=400, width=550)
tiff(filename = "./tif/Fig8.tif", height = 5, width = 6.5, units = "in", res = 300)
hist(df.rps$value, breaks = 50, 
     xlab = "Average Ranked Probability Score (RPS)",
     main = "") #Distribution of RPS across all blocks (n = 502)
dev.off()



# ----------------------- 3-/4-step-ahead prediction ------------------------- #

comp.pred.window <- data.table(OBJECTID=df_wide$OBJECTID,one.ahead.first=scores_mf_first,three.ahead=scores_3ahead_avg,four.ahead=scores_4ahead_avg)
comp.pred.window.melt <- melt(comp.pred.window, id.vars="OBJECTID")

p1<-ggplot(data=comp.pred.window,aes(one.ahead.first,three.ahead)) +
  geom_point() + 
  theme_bw() + 
  theme(legend.position = "none") + 
  xlab("One-month-ahead") + 
  ylab("Three-months-ahead") + 
  geom_abline(intercept = 0, slope = 1)
p2<-ggplot(data=comp.pred.window,aes(one.ahead.first,four.ahead)) +
  geom_point() + 
  theme_bw() + 
  theme(legend.position = "none") + 
  xlab("One-month-ahead") +
  ylab("Four-months-ahead") + 
  geom_abline(intercept = 0, slope = 1)

#png(filename = "./png/Fig9.png", height=300, width=600)
tiff(filename = "./tif/Fig9.tif", height=5, width=10, units = "in", res = 300)
plot_grid(p1,p2,ncol=2, labels=c("A","B"))
dev.off()


# ------------------ QUANTILE INTERVALS FOR EXAMPLE BLOCKS --------------------#

quants_1ahd <- tidy.quants(quants[[54]])
quants_3ahd <- tidy.quants(quants3_final)
quants_4ahd <- tidy.quants(quants4_final)

# Choose example blocks across the range of RPS
View(df.rps)
mean(final.rps.first) #[1] 0.4201546

# Highest RPS = 3.44, for block 2674 "PAKUR"
which(df_wide$OBJECTID=="2674")
# [1] 492
# RPS = 1.000 for block 1901 BHAGWANPUR HAT
which(df_wide$OBJECTID=="1901")
# [1] 427

# df_wide$Block[which(df_wide$OBJECTID=="1901")]
# df_wide$OBJECTID[427]

examples.id <- c(492,427)
#random.id <- sample(1:502,5)
examples <- df_wide$Block[examples.id]

i <- 1
for (blk in examples.id){
print(i)
print(blk)
print(examples[i])
df.1ahd <- make.df(final.OSApred.first[[1]],quants_1ahd,blk,48)
df.3ahd <- make.df(pred3.final[[1]],quants_3ahd,blk,50)
df.4ahd <- make.df(pred4.final[[1]],quants_4ahd,blk,51)

p1 <- quantplot(df.1ahd,"correct1","1-month-ahead",T)
p2 <- quantplot(df.3ahd,"correct1","3-month-ahead",F)
p3 <- quantplot(df.4ahd,"correct1","4-month-ahead",F)

#png(filename = paste0("./png/predwindow_1-4",examples[i],".png"), width=1000, height=300)
tiff(filename = paste0("./tif/predwindow_1-4",examples[i],".tif"), width=14, height=4, units = "in", res = 300)
print(plot_grid(p1, p2, p3, ncol=3))
dev.off()
i <- i+1
}

################################################################################
#                             SUPPLEMENTARY FIGURES                            #
################################################################################

# ------------------ DISTRICTS WITH HIGH ESTIMATED DISPERSION ----------------- #

agg.dist<-aggregate(df_wide[,5:76],by=list(df_wide$District),FUN="sum")
odd.dists<-agg.dist[which(agg.dist$Group.1%in%c("AURANGABAD","BANKA","JEHANABAD","NAWADA")),]
odd.dists.melt<-melt(odd.dists)

#png(filename = "./png/S1Fig.png", height=250,width=1000)
tiff(filename = "./tif/S1Fig.tif", height=3,width=12, units = "in", res = 300)
ggplot(odd.dists.melt,aes(as.numeric(variable),value))+xlab("Month")+ylab("No. reported cases")+geom_line()+facet_grid(~ Group.1)+theme_bw()
dev.off()


# -------- PIT histograms adding lags to a baseline, seasonal model ---------- #

# AIC_seas<-vector(length=12)
# BIC_seas<-vector(length=12)
# scores_seas<-matrix(nrow=12,ncol=4)
# calibp_lag<-matrix(nrow=12,ncol=4)
# fits.seas<-list()
# preds.seas<-list()
# for (i in 1:12){
#   ctl <- list(end = list(f = ~1, offset=population(stsobj)), 
#               ar = list(f = addSeason2formula(~1+t,S=1,period=stsobj@freq)), 
#               ne = list(f = addSeason2formula(~1+t,S=1,period=stsobj@freq), weights = neighbourhood(stsobj) == 1),
#               max_lag=i,
#               subset=13:72,
#               family = "NegBin1")
#   mod<-profile_par_lag(stsobj,control=ctl)
#   osa<-oneStepAhead_hhh4lag(mod, tp=c(48,65), type = "first", which.start = "current")
#   scores_seas[i,]<-colMeans(scores(osa,which=SCORES))
#   calib<-calibrationTest(osa,which="rps", individual=T)
#   calibp_lag[i,]<-calib[["p.value"]]
#   AIC_seas[i]<-AIC(mod)
#   BIC_seas[i]<-BIC(mod)
#   
#   fits.seas<-list.append(fits.seas,mod)
#   preds.seas<-list.append(preds.seas,osa)
# }

#png(filename = "./png/S2Fig.png", height = 500, width = 800)
tiff(filename = "./tif/S2Fig.tif", height = 6, width = 10, units = "in", res = 300)
par(mfrow = c(3,4))
for (osa in preds.seas){
  pit(osa, J = 50)
}
dev.off()


rps <- ggplot(as.data.frame(scores_seas), aes(x=1:12,y=V2)) + 
  geom_line() +
  xlab("Number of lags") +
  ylab("RPS") +
  theme_classic()
cal <- ggplot(as.data.frame(calibp_lag), aes(x=1:12,y=V2)) + 
  geom_line() +
  xlab("Number of lags") +
  ylab("Calibration p-value") +
  theme_classic()
#png(filename = "./png/S2a.png", height = 300, width = 800)
tiff(filename = "./tif/S2a.tif", height = 3, width = 8, units = "in", res = 300)
plot_grid(rps,cal, labels = c("A","B"))
dev.off()

odd.blk.id<-df.rps$id[df.rps$value>2.5]
odd.blk<-which(df_wide$OBJECTID%in%odd.blk.id)
#png("./png/S3Fig.png",height=600,width=900)
tiff(filename = "./tif/S3Fig.tif", height = 6, width = 8, units = "in", res = 300)
par(mfrow=c(2,3))
ind<-1
for (i in odd.blk){plot(c(1:72),cases[,i],type="l",
                        xlab="",
                        xaxt="n",
                        ylab="No. reported cases",
                        main=paste0(df_wide$Block[i]," (",round(df.rps$value[df.rps$id%in%odd.blk.id],2)[ind],")"))
  axis(1, at=yr-1, labels=c(2013:2019), cex.axis=0.8)
  ind<-ind+1}
dev.off()
