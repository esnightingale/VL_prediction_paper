detach(package:plyr)
library(gganimate)
library(gifski)

source("~/VL/surveillance_forecasting/1_data_run_addcovariates.R")
setwd("C:/Users/phpuenig/Dropbox/VL/Monthly prediction/surveillance forecasting/final model selection")

# stsobj2 has populations and blank rows up to end of 2025

subset<-5:156
ctl <- list(end = list(f = ~ -1 + ri(type="iid", corr="all"), offset=population(stsobj2)),
            ar = list(f = addSeason2formula(~1, S=1, period=stsobj2@freq)),
            ne=list(f= addSeason2formula(~1, S=1, period=stsobj2@freq), weights=neighbourhood(stsobj2)==1),
            max_lag = 1,
            subset=subset,
            family = "NegBin1")
model<-profile_par_lag(stsobj2,ctl)

# osa_ri <- oneStepAhead(model, tp=c(61,71), type="first", which.start="current", 
#                      keep.estimates = T)
# mean(scores(osa_ri, which="rps"))

set.seed(101)
nsims<-1000

simfc <- function(model, # model fit to object with NA rows where necessary for months to be forecasted
                    nsims, # Number of simulations to run
                    obsmth, # Months for which data observed
                    fcmth){ # Month up to which forecasts should be made
  
totmth <- obsmth + fcmth
y.start <- model$stsObj@observed[(obsmth-model$max_lag+1):obsmth,]
sim<-array(dim=c(fcmth,ncol(model$stsObj@observed),nsims))

for (n in 1:nsims){
  temp<-simulate(model, subset = (obsmth+1):totmth, y.start = y.start) 
  sim[,,n] <- temp@observed
}
rownames(sim)<-timeall[(obsmth+1):totmth]
colnames(sim)<-df_wide$OBJECTID

return(sim)
}

sim1 <- simfc(model, 1000, 72, 84)
sim2 <- simfc(model, 100, 60, 12)

#save(sim1,file="./Results/simulation/sim_2025.RData")
#save(pop_t, file="./Results/simulation/pops_2025.RData")
#save(sim2,file="./Results/simulation/sim_18.RData")

pathfc <- stepaheadN(model,start=72,type="first",n=84)
simfc <- simfc(model, nsims=1000, obsmth=72, fcmth=3)

# Plot one simulation to check 
obs_sim1<-rbind(cases,sim[,,1])
rownames(obs_sim1)<-timeall
obs_sim1<-melt(obs_sim1)
obs_sim1$pop<-melt(pop)[,3]
obs_sim1$inc<-obs_sim1$value*1e4/obs_sim1$pop
obs_sim1$Var1<-ymd(obs_sim1$Var1)
obs_sim1$obspred<-rep("obs",nrow(obs_sim1))
obs_sim1$obspred[obs_sim1$Var1>ymd("2018-12-01")]<-"pred"

png("./Results/simulation/sim1_2013_2025.png",height=700,width=1200)
ggplot(obs_sim1,aes(Var1,inc,group=Var2))+
  #scale_x_date(limits = as.Date(c('2013-01-01','2025-12-01')))+
  geom_line(col="snow4")+
  geom_hline(yintercept=1/12,col="white",lty="dashed")+
  geom_vline(xintercept=ymd("2019-01-01"), col="indianred")+
  scale_x_date(date_breaks = "1 year", date_labels = "%Y")+
  xlab("Month")+
  ylab("Incidence per 10,000")
dev.off()


