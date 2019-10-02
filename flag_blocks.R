# Author:         Emily S Nightingale
# Institutions:   London Schoool of Hygiene and Tropical Medicine, London, UK
# Date Published: XX July 2019
################################################################################
# 
################################################################################
################################################################################

library(tmap)
library(ggplot2)
library(reshape2)
library(surveillance)
library(hhh4addon)
library(cowplot)

source('~/VL/data_setup/1_data_run_addcovariates.R')
source('~/VL/surveillance_modelling/NStepAhead.R')
setwd("C:/Users/phpuenig/Dropbox/VL/Monthly prediction/surveillance forecasting")

# Final model with four month lags and simple seasonality in AR and NE
# subset<-5:72
# ctl <- list(end = list(f = ~ 1, offset=population(stsobj)),
#             ar = list(f = addSeason2formula(~1, S=1, period=stsobj@freq)),
#             ne=list(f= addSeason2formula(~1, S=1, period=stsobj@freq), 
#                     weights=neighbourhood(stsobj)==1),
#             max_lag = 4,
#             subset=subset,
#             family = "NegBin1")
# model<-profile_par_lag(stsobj,ctl)
#save(model, file="final_model.RData")

# Three-month-ahead path forecasts from June 2018
start <- 60
preds3 <- stepaheadN(model, start=start, type="first", n=3)
save(preds3, file="./predictions/preds3_mardec18.RData")
preds3_q10 <- stepaheadN(model, start=start, type="first", n=3, prob=0.1)
save(preds3_q10, file="./predictions/preds3_q10_mardec18.RData")


load("final_model.RData")
load("./predictions/preds3_mardec18.RData")
load("./predictions/preds3_q10_mardec18.RData")

# Compare predicted 90% quantiles to observed case counts and flag blocks for 
# which observed exceeds expected upper bound
q90 <- preds3[[3]]
q10 <- preds3_q10[[3]]
flag.hi <- (cases[(start+3):nrow(cases),] > q90 & 
            cases[(start+3):nrow(cases),] > 1)
flag.lo <- (cases[(start+3):nrow(cases),] < q10)

flag.hilo <- (flag.hi | flag.lo)
which(flag.hilo[t,]==TRUE)
all.flagged <- which(colSums(flag.hilo) > 0)
                      
# Plot observations up to fit and prediction
rownames(q10) <- timeall[(start+3):nrow(cases)]
rownames(q90) <- timeall[(start+3):nrow(cases)]

plotdata <- setNames(melt(cases[61:72,]), c("Month","Block","Obs"))
plotdata <- merge(plotdata, melt(q10), by.x=c("Month","Block"),
                  by.y=c("Var1","Var2"), all.x=T)
plotdata <- merge(plotdata, melt(q90), by.x=c("Month","Block"),
                  by.y=c("Var1","Var2"), all.x=T)
names(plotdata)[4:5] <- c("q10","q90")
plotdata$Month <- ymd(plotdata$Month)
plotdata <- mutate(plotdata, correct = (Obs <= q90 & Obs >= q10))
plotdata$correct[is.na(plotdata$correct)] <- 1

plot_flagged <- function(data, flags){
p <- ggplot(data[flags,], aes(Month, Obs, group = Block)) + 
      geom_ribbon(aes(x = Month, ymin = q10, ymax = q90), 
                  alpha = 0.5, fill = "forestgreen") + 
      geom_point(aes(shape = as.factor(correct))) +
      scale_shape_manual(values = c(4, 19)) +
      facet_wrap(vars(Block), scales = "free_y") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = -1),
            axis.title.x = element_blank(),
            legend.position = "none") +
      ylab("No. reported cases")
return(p)
}
plot_map <- function(mapdata, flag.hi, flag.lo, start){
  mapdata$flag[mapdata$region%in%names(which(flag.hi==TRUE))] <- "High"
  mapdata$flag[mapdata$region%in%names(which(flag.lo==TRUE))] <- "Low"
  map <- ggplot(mapdata, aes(long, lat, group=group, fill = flag)) + 
    geom_polygon() +
    geom_path() +
    scale_fill_manual(values = c("red","blue","white"), 
                      labels = c("Above expected", "Below expected"),
                      na.translate=F) +
    theme(legend.title=element_blank(), 
          legend.position = "bottom",
          legend.text = element_text(size = 10),
          axis.text = element_blank(), 
          axis.ticks = element_blank(),
          axis.line=element_blank()) +
    xlab("") + ylab("") +
    ggtitle(timeall[start+2+t]) +
    coord_quickmap()
  return(map)
}
# xlab(paste(min(plotdata2$Month),max(plotdata2$Month),sep=" - ")) +

t <- 1
flag.t <- plotdata$Block%in%names(which(flag.hilo[t,]==TRUE))

pdf(paste0("./figures/flagged_",timeall[start+3],"_",timeall[start+2+nrow(q90)],".pdf"), height = 7, width = 15)
for (t in 1:nrow(q90)){
  blks <- plot_flagged(plotdata, flag.hilo[t,])
  map <- plot_map(map_data(shp), flag.hi[t,], flag.lo[t,], start)
print(plot_grid(map,blks))
}
dev.off()


################################################################################
# Reporting delay
################################################################################

load("C:/Users/phpuenig/Dropbox/VL/Data/KAMIS/Clean/linelist_1318_clean_VL.Rdata") #42464 cases

vl$delay <- vl$diag.date - vl$reg_date
head(vl$delay)
hist(as.numeric(vl$delay, breaks=100))

# run beginning of data_setup to check merging of patient and diagnosis records
alldata.vl <- alldata[alldata$case_type=="VL",]

# want to check delay from first diagnosis to registration date
alldata.vl <- arrange(alldata.vl,res_patient_code_original,date_of_diagnosis)
alldata.vl <- alldata.vl[!duplicated(alldata.vl$res_patient_code_original),]
alldata.vl$delay <- ymd(alldata.vl$reg_date) - ymd(alldata.vl$date_of_diagnosis)
head(alldata.vl$delay)
hist(as.numeric(alldata.vl$delay[abs(alldata.vl$delay)<60], breaks=200))


# par(mfrow=c(1,2))
# for (t in 1:nrow(q90)){
#   clr <- rep('white', length(shp))
#   clr[flag.hi[t,]] <- 'red'
#   clr[flag.lo[t,]] <- 'blue'
#   plot(shp, col = clr, main=timeall[68+t])
#   #legend(x="topright", legend = c("Above expected", "Below expected"), 
#   #       fill = c("red","blue"), border = "white")
#   flag <- plotdata$Block%in%names(which(flag.hilo[t,]==TRUE))
#   plot_flagged(plotdata, flag)
# }

# plotdata1 <- plotdata[plotdata$Block=="1411",]
# ggplot(plotdata1, aes(Month, Obs)) + 
#   geom_point() + 
#   geom_ribbon(aes(x = 1:n_distinct(Month), ymin = q10, ymax = q90), alpha = 0.5) 
