# Author:         Emily S Nightingale
# Institutions:   London Schoool of Hygiene and Tropical Medicine, London, UK
# Date Published: XX July 2019
################################################################################
# Produce figures in manuscript
################################################################################
################################################################################

library(ggrepel)
library(viridis)

source("C:/Users/phpuenig/Documents/VL/surveillance_modelling/1_data_run.R")
source("C:/Users/phpuenig/Documents/VL/surveillance_modelling/plot_hhh4_amended.R")

set.seed(101)
SCORES <- c("logs", "rps", "dss", "ses")

path <- "C:/Users/phpuenig/Dropbox/VL/Monthly Prediction/surveillance forecasting/final model selection/Results/"
load(paste0(path,"models_osafirst.Rdata"))
load(paste0(path,"osa_first_osafirst.Rdata"))
load(paste0(path,"osa_rolling_osafirst.Rdata"))
load(paste0(path,"scores_first_osafirst.RData"))
load(paste0(path,"scores_rolling_osafirst.RData"))
load(paste0(path,"quants_osafirst.Rdata"))
load(paste0(path,"U_osafirst.RData"))
load(paste0(path,"pred3_final_osafirst.RData"))
load(paste0(path,"pred4_final_osafirst.RData"))
load(paste0(path,"quants3_final.RData"))
load(paste0(path,"quants4_final.RData"))
load(paste0(path,"scores_mf_first.RData"))
load(paste0(path,"scores_mf_rolling.RData"))
load(paste0(path,"scores_3ahead_final.RData"))
load(paste0(path,"scores_4ahead_final.RData"))

setwd("C:/Users/phpuenig/Dropbox/VL/Monthly prediction/surveillance forecasting/final model selection/figures/paper")

#-----------------------DESCRIPTIVE-----------------------# 

timeall <- format(seq(as.Date("2013-1-1"), by="month", length.out = 73), "%Y-%m-%d")
yr <- seq(1,73,12)

par(mfrow=c(1,1))
# png(filename="timeseries_jan13dec18.png", width=700, height=500)
tiff(filename="Fig1.tif",height=4.5,width=5.2, units = "in", res = 300)
barplot(rowSums(cases), xaxt="n", xlab="", ylab="No. reported cases", space=0)
axis(1, at=yr-1, labels=timeall[yr], las=2 , cex.axis=0.8)
dev.off()

inc18 <- data.frame(id=df_wide$OBJECTID, value=colSums(cases[61:72,])*1e4/colMeans(pops[61:72,]))
tiff(filename = "Fig2.tiff", height = 4.5, width = 5.2, units = "in", res = 300)
mapplot(shapefile=shp, data=inc18, legend_title="Cases per 10,000", value = value, pal = "C") +
  ggtitle("Total block-level incidence per 10,000 in 2018", )
dev.off()

# Monthly average per 10,000 population
# png(filename="month_inc_map_13-18.png", width=600, height=500)
inc18mth <- data.frame(id=df_wide$OBJECTID, value=colMeans(cases[61:72,]*1e4/pops[61:72,]))
tiff(filename = "Fig2a.tiff", height = 4.5, width = 5.2, units = "in", res = 300)
mapplot(shapefile=shp, data=inc18mth, legend_title="Cases per 10,000", value = value, pal = "C") +
  ggtitle("Block-level, average monthly incidence per 10,000 over 2018")
dev.off()

blocktotal <- colSums(cases)
#png(filename="hist_block_total.png", width=600, height=500)
tiff(filename = "Fig2b.tiff", height = 4.5, width = 5.2, units = "in", res = 300)
hist(blocktotal, breaks=30, xlab="Total cases Jan 2013-Dec 2018", 
     main="Block-level distribution of cases over whole time period (N = 502)")
dev.off()

#----------------- RANDOM MODEL ASSESSMENT ---------------------#

load("C:/Users/phpuenig/Dropbox/VL/Monthly prediction/surveillance forecasting/exploratory/random model draw/result_table_randmods30.RData")
results <- result_table

fitvfirst <- ggplot(results,aes(AIC_train,rps.first, col=k)) +
  theme(legend.position = c(0.6,0.2)) +
  labs(colour="No. parameters") +
  geom_point() +
 # geom_label_repel(aes(label = Model)) +
  xlab("AIC") + ylab("RPS") + 
  ggtitle("Training period fit") +
  scale_color_viridis_c(option = "C") 
firstvroll <- ggplot(results,aes(rps.first,rps.roll, col=k)) +
  theme(legend.position = "none") +
  geom_point() +
 # geom_label_repel(aes(label = Model)) +
  xlab("RPS (fixed)") + ylab("RPS (rolling)") + 
  ggtitle("Training period fit versus rolling updates") +
  scale_color_viridis_c(option = "C") 

tiff(filename = "fitvspred_randmods30.png", height=5, width=10, units = "in", res=300)
plot_grid(fitvfirst+geom_point(cex=3), firstvroll+geom_point(cex=3), labels=c("A", "B"), ncol = 2)
dev.off()

tiff(filename="Fig3.tif",height=5.2,width=5.2, units = "in", res = 300)
fitvfirst + 
  #scale_color_gradient(low=lshtm[4], high=lshtm[3]) +
  #scale_color_gradient(low=lshtm[2], high=lshtm[4]) +
  geom_point(cex=3) 
dev.off()


#----------------- SYSTEMATIC SELECTION ------------------------#

load("C:/Users/phpuenig/Dropbox/VL/Monthly prediction/surveillance forecasting/final model selection/Results/selected_models_osafirst.RData")
test <- c(48,71)

selected <- c(1,6,23,33,42,52)
# png(filename ="PIT_selected_osafirst.png", height=600, width=900)
tiff(filename = "FigS2.tiff", height = 5, width = 6.5, units = "in", res = 300)
par(mfrow=c(2,3))
i <- 1
for (m in selected.models){
  osa <- oneStepAhead_hhh4lag(m,tp=test,type="first",which.start = "current")
  pit(osa, J=30, plot=list(ylim=c(0,2), main=paste0("Model no: ", selected[i])))
  i <- i+1
}
dev.off()


load("C:/Users/phpuenig/Dropbox/VL/Monthly Prediction/surveillance forecasting/final model selection/results/result_table_osafirst.Rdata")
result_table$Miscal.first <- factor(result_table$Miscal.first, levels = c("No evidence", "Borderline", "Strong"))

# Fit to the training data 
# png(filename = "fit_vs_rpsfirst.png", height=600, width=700)
tiff(filename = "Fig4.tiff", height = 5.2, width = 5.2, units = "in", res = 300)
ggplot(result_table,aes(AIC_train,rps.first,colour=as.factor(Miscal.first))) +
  geom_point(cex=3) + theme_bw() + labs(colour="Miscalibration") +
  scale_colour_manual(values = plasma(3)) +
 # geom_label_repel(aes(label = Model)) + 
  xlab("Fit to training data (AIC)") +
  ylab("Predictive power (RPS)")
dev.off()

results <- result_table
results$highlight[results$selected==0] <- 16
results$highlight[results$selected==1] <- 5
results$highlight[c(1,52)] <- 16
results$highlight[42] <- 8
#palette(c("black","dodgerblue","forestgreen","firebrick","orange"))
palette(viridis(5))

nmod <- nrow(results)

# png(filename = "pred_fit_utility_nopois_osafirst.png", width=700, height=600)
tiff(filename = "Fig5.tiff", height = 6.5, width = 7.5, units = "in", res = 300)
par(mfrow=c(2,2), mar=c(4,4,3,4))

plot(c(1:nmod), results$rps.first, 
     main="RPS", ylab="", xlab="Model no.", cex=0.7,
     pch=results$highlight, col=as.factor(results$stage))
legend(x=30, y=0.65, legend=c("1","2","3","4 (final)"),
       title="Selected model",
       col=viridis(5)[1:4],
       pch=c(5,5,5,8),
       cex = 0.7,
       bty="n")
plot(c(1:nmod), results$AIC_train, 
     main="AIC (training set)", ylab="", xlab="Model no.", cex=0.7, 
     pch=results$highlight, col=as.factor(results$stage))
plot(c(1:nmod), results$U2575, xlim=c(1,nmod), 
     main="Empirical coverage - 25-75%", ylab="", xlab="Model no.", 
     pch=results$highlight, cex=0.7, 
     col=as.factor(results$stage))
par(new = T)
plot(c(1:nmod), results$U2575_qwd, xlim=c(1,nmod), axes=F, xlab=NA, ylab=NA, 
     type="l", lty="dashed", col="darkgrey")
axis(side = 4)
mtext(side = 4, line = 2, 'Average interval width')
plot(c(1:nmod), results$U1090, 
     main="Empirical coverage - 10-90%", ylab="", xlab="Model no.", cex=0.7,
     pch=results$highlight,col=as.factor(results$stage))
par(new = T)
plot(c(1:nmod), results$U1090_qwd, xlim=c(1,nmod), axes=F, xlab=NA, ylab=NA,
     type="l", lty="dashed", col="darkgrey")
axis(side = 4)
mtext(side = 4, line = 2, 'Average interval width')
dev.off()


## ----------------------------- FINAL MODEL -------------------------------- ##

setwd("C:/Users/phpuenig/Dropbox/VL/Monthly Prediction/surveillance forecasting/
      final model selection/figures/paper")

#------------------------------------------------------------------------------#

f <- 42  #index of final model
t.pred <- obsmth-48  #length of test period
nblock <- dim(cases)[2]
m.final <- model.list.all[[f]]
final.OSApred <- osa.final.rolling[[f]]
final.OSApred.first <- osa.final.first[[f]]
final.rps.first <- models.scores.first[[f]][,,2]
final.rps.roll <- models.scores.rolling[[f]][,,2]
U1_mat <- U.osafirst[[f]][[3]][,,1]
U2_mat <- U.osafirst[[f]][[3]][,,2]
U3_mat <- U.osafirst[[f]][[3]][,,3]

# FIT FOR HIGHEST INCIDENCE BLOCKS
df.rps <- data.frame(id=df_wide$OBJECTID, value=colMeans(final.rps.first))
avginc <- colMeans(inc)
highinc <- names(sort(avginc,decreasing = T))[1:4]
df_wide$Block[which(df_wide$OBJECTID%in%highinc)] 
#"GOPIKANDAR"   "KATHIKUND"    "BOARIJOR"     "SUNDARPAHARI"
blknames <- paste0(df_wide$Block[which(df_wide$OBJECTID%in%highinc)],
                 " (RPS = ",
                 round(df.rps$value[df.rps$id%in%highinc],2),")")
# png(filename = "final_fit_highinc4.png", height=600, width=800)
palette(magma(8))
tiff(filename = "Fig6.tiff", height = 6, width = 7.5, units = "in", res = 300)
plot_hhh42(m.final, units=highinc, names=blknames, ylab="No. reported cases",
           col=c("forestgreen","skyblue", "orange"))
dev.off()


## QUANTILE INTERVALS FOR EXAMPLE BLOCKS

tidy.quants <- function(quants){
t.pred <- nrow(quants[[1]])
nblock <- ncol(quants[[1]])
qsLow_final<-array(dim=c(3,t.pred,nblock))
qsHi_final<-array(dim=c(3,t.pred,nblock))
qsLow_final[1,,] <- quants[[1]]
qsHi_final[1,,] <- quants[[6]]
qsLow_final[2,,] <- quants[[2]]
qsHi_final[2,,] <- quants[[5]]
qsLow_final[3,,] <- quants[[3]]
qsHi_final[3,,] <- quants[[4]]

qs_final<-list(qsLow_final,qsHi_final)  
return(qs_final)
}

quants_1ahd <- tidy.quants(quants_OSAfirst[[42]])
quants_3ahd <- tidy.quants(quants3_final)
quants_4ahd <- tidy.quants(quants4_final)

# Choose example blocks across the range of RPS
View(df.rps)
mean(final.rps.first) #[1] 0.4201441

# Highest RPS = 3.47, for block 2674 Pakur
which(df_wide$OBJECTID=="2674")
# [1] 492
# RPS = 0.4196 for block 1878 "PURNAHI"
which(df_wide$OBJECTID=="1878")
# [1] 404
# RPS = 1.001 for block 1857 ISHUAPUR
which(df_wide$OBJECTID=="1857")
# [1] 383
which(df_wide$OBJECTID=="1908")
# [1] 434

df_wide$Block[which(df_wide$OBJECTID=="1454")]
df_wide$OBJECTID[222]

examples.id <- c(383,404,492,434)
examples.id <- sample(1:491,5)
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

tiff(filename = paste0("./quant_plots/predwindow_1-4",examples[i],".png"), width=1000, height=300)
print(plot_grid(p1, p2, p3, ncol=3))
dev.off()
i <- i+1
}

#: Jalalpur, Bihar (Average RPS = 0.5)

#pdf(file="./figures/assessing prediction/obsvspred_werrbars_3qu.pdf", height=6,width=8)
# png(file="./figures/assessing prediction/obsvspred_jokihaat_1090.png", height=400,width=500)
ggplot(df.pred1,aes(Month,Obs,col=as.factor(correct1)), ylab="No. reported cases") +
  geom_point(pch=19) + 
  geom_point(aes(Month,Pred), col="forestgreen", pch=4, cex=1.2) +
  geom_errorbar(aes(ymax=q10,ymin=q90),col="forestgreen") +
  scale_color_manual(values=c("orange","black")) + ylab("No. reported cases") + 
  ylab("No. reported cases") +
  theme(legend.position="none") +
  ggtitle("Predicted distribution: Q10-Q90") + theme_bw()
# dev.off()
# png(file="./figures/assessing prediction/obsvspred_jokihaat_2575.png", height=400,width=500)
ggplot(df.pred1,aes(Month,Obs,col=as.factor(correct2)), ylab="No. reported cases") +
  geom_point(pch=19) + 
  geom_point(aes(Month,Pred), col="forestgreen", pch=4, cex=1.2)+
  geom_errorbar(aes(ymax=q75,ymin=q25),col="forestgreen") + 
  scale_color_manual(values=c("orange","black")) + ylab("No. reported cases") + 
  ylab("No. reported cases") +
  theme(legend.position="none") +
  ggtitle("Predicted distribution: Q25-Q75") + theme_bw()
# dev.off()
# png(file="./figures/assessing prediction/obsvspred_jokihaat_4555.png", height=400,width=500)
ggplot(df.pred1,aes(Month,Obs,col=as.factor(correct3)), ylab="No. reported cases") +
  geom_point(pch=19) + geom_point(aes(Month,Pred), col="forestgreen", pch=4, cex=1.2)+
  geom_errorbar(aes(ymax=q55,ymin=q45),col="forestgreen") + 
  scale_color_manual(values=c("orange","black")) + ylab("No. reported cases") + 
  ylab("No. reported cases") +
  theme(legend.position="none") +
  ggtitle("Predicted distribution: Q45-Q55") + theme_bw()
# dev.off()

# Example OSA predictions
png(filename = "fanplot_jalalpur_1ahd.png", height=400,width=500)
ggplot(df.pred,aes(Month,Obs), ylab="No. reported cases") +
  geom_ribbon(aes(ymin=q10,ymax=q90),lty=0,fill="gold",alpha=0.5) +
  geom_ribbon(aes(ymin=q25,ymax=q75),lty=0,fill="orange",alpha=0.5) +
  geom_ribbon(aes(ymin=q45,ymax=q55),lty=0,fill="red",alpha=0.5) +
  geom_line(aes(Month,Pred),col="white",lwd=1.2) +
  geom_point(pch=19) +
 # scale_color_manual(values=c("orange","black")) +
  ylab("No. reported cases") + 
  theme(legend.position="none") +
  ggtitle("1-month-ahead prediction for Jalalpur, Bihar (Average RPS = 0.5)")
dev.off()


# ------------------------ 3 month ahead prediction ---------------------------#


permut.test2 <- function(null.scores,new.scores){
  set.seed(seed) 
  models.scores<-list(null=null.scores,new=new.scores)
  permut.test <- permutationTest(
    models.scores$null,
    models.scores$new,
    nPermutation = 10000)
  print(permut.test)
}
permut.test2(models.scores.first[[42]][-c(1,2),,2],scores_final_3ahead)
# $`diffObs`
# [1] -0.02443085
# 
# $pVal.permut
# [1] 9.999e-05
# 
# $pVal.t
# [1] 1.955191e-19
permut.test2(models.scores.first[[42]][-c(1:3),,2],scores_final_4ahead)
# $`diffObs`
# [1] -0.02747089
# 
# $pVal.permut
# [1] 9.999e-05
# 
# $pVal.t
# [1] 9.230229e-17
permut.test2(scores_final_3ahead[-1,],scores_final_4ahead)
# $`diffObs`
# [1] -0.004419976
# 
# $pVal.permut
# [1] 0.02529747
# 
# $pVal.t
# [1] 0.02354382

comp.pred.window <- data.table(OBJECTID=df_wide$OBJECTID,one.ahead.first=scores_mf_first,one.ahead.roll=scores_mf_rolling,three.ahead=scores_3ahead_avg,four.ahead=scores_4ahead_avg)
comp.pred.window.melt <- melt(comp.pred.window, id.vars="OBJECTID")
boxplot(value ~ variable, data=comp.pred.window.melt, ylab="RPS")
p1<-ggplot(data=comp.pred.window,aes(one.ahead.first,three.ahead)) +
  geom_point() + 
  theme_bw() + theme(legend.position = "none") + 
  xlab("One-month-ahead") + ylab("Three-months-ahead") + 
  geom_abline(intercept = 0, slope = 1)
p2<-ggplot(data=comp.pred.window,aes(one.ahead.first,four.ahead)) +
  geom_point() + 
  theme_bw() + theme(legend.position = "none") + xlab("One-month-ahead") +
  ylab("Four-months-ahead") + geom_abline(intercept = 0, slope = 1)
png(filename = "ahead3-4_vs_OSAroll.png", height=500, width=1000)
plot_grid(p1,p2,ncol=2, labels=c("A","B"))
dev.off()

png(filename = "scores_ahead3-4.png", height=500, width=1000)
ggplot(data=comp.pred.window,aes(three.ahead,four.ahead)) + 
  geom_point() +
  theme_bw() + theme(legend.position = "none") +
  xlab("Three months ahead") + ylab("Four months ahead")
dev.off()


# 
# cor(scores_mf_first,scores_3ahead_avg) #0.9867457
# cor(scores_mf_first,scores_4ahead_avg) #0.9874804
# cor(scores_mf_rolling,scores_3ahead_avg) #0.9867032
# cor(scores_mf_rolling,scores_4ahead_avg) #0.9873889
# 
# mean(models.scores.first[[42]][,,2]) # 0.4201441
# mean(models.scores.rolling[[42]][,,2]) # 0.4197563
# mean(scores_3ahead_avg) # 0.4413778
# mean(scores_4ahead_avg) # 0.438052
# mean(scores_3ahead_avg-models.scores.rolling[[42]][-c(1:3),,2]) #0.03119174
# mean(scores_4ahead_avg-models.scores.rolling[[42]][-c(1:3),,2]) #0.02786594
# 


comp.scores <- data.frame(osa_first=melt(score_m1_first)[,3],
                          osa_rolling=melt(score_m1_rolling)[,3],
                          ahead3=melt(scores_3ahead)[,3],
                          ahead4=melt(rbind(matrix(NA,nrow=1,ncol=nblock),
                                            scores_4ahead))[,3])
comp.scores.block <- data.frame(osa_first=colMeans(score_m1_first),
                                osa_rolling=colMeans(score_m1_rolling),
                                ahead3=colMeans(scores_3ahead),
                                ahead4=colMeans(rbind(matrix(NA,nrow=1,ncol=nblock),
                                                      scores_4ahead)))

comp.scores.final <- data.frame(block=melt(score_mf_first)[,2],
                                osa_first=melt(score_mf_first)[,3],
                                osa_rolling=melt(score_mf_rolling)[,3],
                                ahead3=melt(scores_final_3ahead)[,3], 
                                ahead4=melt(rbind(matrix(NA,nrow=1,ncol=nblock),
                                                  scores_final_4ahead))[,3])
comp.scores.block.final <- data.frame(osa_first=colMeans(score_mf_first),
                                      osa_rolling=colMeans(score_mf_rolling),
                                      ahead3=colMeans(scores_final_3ahead),
                                      ahead4=colMeans(scores_final_4ahead))

png(filename = "scores_ahead3-4.png", height=600, width=600)
ggplot(data=comp.scores.block.final,aes(ahead3,ahead4)) + geom_point() +
  theme_bw() + theme(legend.position = "none") + 
  xlab("Three months ahead") + ylab("Four months ahead")
dev.off()


# ------------------------------------- DISTRICTS WITH HIGH ESTIMATED DISPERSION ---------------------------------------- #

highinc.ind<-which(names(cases)%in%highinc)

m1d2<-update(model.list.all[[1]], family=as.factor(df_wide$District))
agg.dist<-aggregate(df_wide[,5:76],by=list(df_wide$District),FUN="sum")
odd.dists<-agg.dist[which(agg.dist$Group.1%in%c("AURANGABAD","BANKA","JEHANABAD","NAWADA")),]
odd.dists.melt<-melt(odd.dists)
png(filename = "odd_district_disp.png", height=300,width=900)
ggplot(odd.dists.melt,aes(as.numeric(variable),value))+xlab("Month")+ylab("No. reported cases")+geom_line()+facet_grid(~ Group.1)+theme_bw()
dev.off()
# par(mfrow=c(1,3))
# for (i in odd.dists){plot(c(1:66),agg.dist[i,-1],type="l", main=agg.dist$Group.1[i], xlab="Month", ylab="No. reported cases")}
# mtext()

par(mfrow=c(2,3))
rand<-sample(c(1:36)[-odd.dists],6)
for (i in rand){plot(c(1:66),agg.dist[i,-1],type="l",main=agg.dist$Group.1[i])}

odd.blk.id<-df.rps$id[df.rps$value>2.5]
odd.blk<-which(df_wide$OBJECTID%in%odd.blk.id)
png("high_rps_blks.png",height=600,width=900)
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

# Compare patterns with random selection of blocks
par(mfrow=c(2,3))
rand<-sample(c(1:491)[-odd.blk],6)
for (i in rand){plot(c(1:66),cases[,i],type="l",main=df_wide$Block[i])}

