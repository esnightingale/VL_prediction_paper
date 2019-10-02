#---------------------------------------------------------------------------- #

# Investigate potential impact of % urban population, onset to diag and diag to
# treatment as covariates.
# 
#---------------------------------------------------------------------------- #

source("C:/Users/phpuenig/Documents/VL/data_setup/1_data_run_addcovariates.R")

# Location for saving results
setwd("C:/Users/phpuenig/Dropbox/VL/Monthly Prediction/surveillance forecasting/
      fits")

seed=101
probs=list(0.25,0.75)
SCORES <- c("logs", "rps", "dss", "ses")

# Run all models on a subset in order to later compare with up to 4th order 
# lagged model.
subset<-5:72

# Define test period for evaluating predictive power
test<-c(48,71)

nb <- poly2nb(shp)
head(nb)
coord <- coordinates(shp)
shp$long <- coord[, 1]
shp$lat <- coord[, 2]

library(sf)
mapsf <- st_as_sf(shp)
mapsf$urban <- blockdata2$urban_propn
ggplot(mapsf) + geom_sf(aes(fill = urban)) +
  geom_text(aes(long, lat, label = admin3), color = "white", cex=0.6) +
  theme_bw() + guides(fill = FALSE)

hist(SIR,breaks=50)
avgSIR <- colMeans(SIR)
head(sort(avgSIR, decreasing=T))

highSIR <- names(head(sort(avgSIR, decreasing=T)))
View(blockdata[blockdata$OBJECTID%in%highSIR,])
ggplot(input2[input2$OBJECTID%in%highSIR,],
       aes(interval_start,count,group=Block,colour=Block)) +
       geom_line()

# SIR over 2018
mapsf$SIR <- colMeans(SIR[61:72,])
ggplot(mapsf) + geom_sf(aes(fill = SIR)) +
  geom_text(aes(long, lat, label = admin3), color = "white", cex=0.6) +
  theme_bw() + 
  scale_fill_gradient2(
    midpoint = 1, low = "blue", mid = "white", high = "red"
  )

# Incidence relative to target 2018

mapsf$SIR.tgt18 <- SIR.tgt18
ggplot(mapsf) + 
  geom_sf(aes(fill = SIR.tgt18)) +
  #geom_text(aes(long, lat, label = admin3), color = "grey", cex=1) +
  theme_bw() + 
  scale_fill_gradient2(
    midpoint = 1, low = "blue", mid = "white", high = "red"
  ) 
  


#---------------------------------------------------------------------------- #

ggplot(input2,aes(jitter(urban_propn),jitter(count)))+geom_point()
ggplot(input2,aes(jitter(urban_propn),jitter(incidence)))+geom_point()

# Highest incidence comes from blocks with no urban population

blkavg_urban <- group_by(input2, State, District, Block) %>%
                summarise(avg_urban_propn=mean(urban_propn))

#---------------------------------------------------------------------------- #

c1 <- list(end = list(f = ~1, offset=population(stsobj)),
            ar=list(f=addSeason2formula(~1, S=1, period=stsobj@freq)),
            ne=list(f=addSeason2formula(~1, S=1, period=stsobj@freq)),
            max_lag = 4,
            subset=subset,
            family = "NegBin1")
m1 <- profile_par_lag(stsobj,control=c1)

c2 <- list(end = list(f = ~1+urbanprop, offset=population(stsobj)),
           ar=list(f=addSeason2formula(~1, S=1, period=stsobj@freq)),
           ne=list(f=addSeason2formula(~1, S=1, period=stsobj@freq)),
           max_lag = 4,
           subset=subset,
           family = "NegBin1")
m2 <- profile_par_lag(stsobj,control=c2)

c3 <- list(end = list(f = ~1, offset=population(stsobj)),
           ar=list(f=addSeason2formula(~1+urbanprop, S=1, period=stsobj@freq)),
           ne=list(f=addSeason2formula(~1, S=1, period=stsobj@freq)),
           max_lag = 4,
           subset=subset,
           family = "NegBin1")
m3 <- profile_par_lag(stsobj,control=c3)

c4 <- list(end = list(f = ~1, offset=population(stsobj)),
           ar=list(f=addSeason2formula(~1, S=1, period=stsobj@freq)),
           ne=list(f=addSeason2formula(~1+urbanprop, S=1, period=stsobj@freq)),
           max_lag = 4,
           subset=subset,
           family = "NegBin1")
m4 <- profile_par_lag(stsobj,control=c4)

model.list <- list(m1,m2,m3,m4)
save(model.list, file="test_covs.Rdata")

select <- modelassess(model.list, test, type="first")
scores <- score_tidy(select)

#png(file="PIT_select2.png",height=700,width=1000)
par(mfrow=c(2,2))
for (m in model.list){
  osa <- oneStepAhead_hhh4lag(m,tp=test,type="first",which.start = "current")
  pit(osa,J=50)
}
#dev.off()

# ------------------------- Onset->diagnosis->treatment ------------------------ #

summary(as.numeric(vl_blk$days_fever_before_diag))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    0.00   20.00   26.00   33.38   32.00  730.00   10617

summary(as.numeric(vl_blk$diagtotrt))
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max.      NA's 
# -1689.000     0.000     0.000     5.299     2.000  1530.000      2982 


cov_byblk <-  group_by(vl_blk,State,District,Block) %>%
           summarise(daysfever=median(days_fever_before_diag,na.rm=T),
                     diagtotrt=median(diagtotrt,na.rm=T))

cov_byblk <- merge(cov_byblk,mtchdNames[,c(7:9,15)],by.x=c("State","District","Block"),by.y=c("kamis_master_state","kamis_master_dist","kamis_master_block"),all=T)
cov_byblk <- merge(cov_byblk,shp@data,by="dot.name",all=T)

df.daysfever <- data.frame(id=cov_byblk$OBJECTID,value=cov_byblk$daysfever)
mapplot(shp,df.daysfever,"Days fever before diagnosis")

df.diagtotrt <- data.frame(id=cov_byblk$OBJECTID,value=as.numeric(cov_byblk$diagtotrt))
# Remove outlying value to clarify colour scale
df.diagtotrt$value[df.diagtotrt$value>300]<-NA 
mapplot(shp,df.diagtotrt,"Diagnosis to treatment")

# -------------------------------------- By sex ----------------------------- #

sexratio <-  group_by(input_s,State,District,Block,) %>%
  summarise(daysfever=median(days_fever_before_diag,na.rm=T),
            diagtotrt=median(diagtotrt,na.rm=T))

