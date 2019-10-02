setwd("C:/Users/phpuenig/Dropbox/VL/Monthly prediction/surveillance forecasting/final model selection")
source("data_run_update2604.R")

# stsobj2 has populations and blank rows up to end of 2025

subset<-5:156
ctl <- list(end = list(f = ~1, offset=population(stsobj2)),
            ar = list(f = addSeason2formula(~1, S=1, period=stsobj2@freq)),
            ne=list(f= addSeason2formula(~1, S=1, period=stsobj2@freq), weights=neighbourhood(stsobj2)==1),
            max_lag = 4,
            subset=subset,
            family = "NegBin1")
model<-profile_par_lag(stsobj2,ctl)

set.seed(101)
nsims<-1000
sim<-array(dim=c(fcmth,ncol(cases2),nsims))
for (n in 1:nsims){
  temp<-simulate(model, subset = (obsmth+1):totmth, y.start = stsobj2@observed[69:72,]) 
  sim[,,n] <- temp@observed
}
rownames(sim)<-timeall[-(1:obsmth)]
colnames(sim)<-df_wide$OBJECTID
#save(sim,file="./Results/simulation/sim_2025.RData")
#save(pop_t, file="./Results/simulation/pops_2025.RData")

load("./Results/simulation/sim_2025.RData")

# Plot one simulation to check 
obs_sim1<-rbind(cases,sim[,,1])
rownames(obs_sim1)<-timeall
obs_sim1<-melt(obs_sim1)
obs_sim1$pop<-melt(pop_t)[,3]
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


# Aggregate by year


sim_agg<-array(dim=c(7,dim(sim)[2:3]))
colnames(sim_agg)<-colnames(sim)
rownames(sim_agg)<-2019:2025
for (n in 1:1000){ 
  sim_agg[7,,n]<-colSums(sim[73:84,,n])
  sim_agg[6,,n]<-colSums(sim[61:72,,n])
  sim_agg[5,,n]<-colSums(sim[49:60,,n])
  sim_agg[4,,n]<-colSums(sim[37:48,,n])
  sim_agg[3,,n]<-colSums(sim[25:36,,n])
  sim_agg[2,,n]<-colSums(sim[13:24,,n])
  sim_agg[1,,n]<-colSums(sim[1:12,,n])
} 
save(sim_agg,file="./Results/simulation/sim_agg.RData")


pop_agg<-melt(pop_t)
pop_agg$Year<-substr(pop_agg$Var1,1,4)
pop_agg<-aggregate(pop_agg$value,by=list(pop_agg$Var2,pop_agg$Year),FUN="mean")
pop_agg<-dcast(pop_agg,Group.2~Group.1,value.var="x")
rownames(pop_agg)<-pop_agg$Group.2
pop_agg<-pop_agg[,-1]
save(pop_agg,file="./Results/simulation/pop_agg.RData")

# sim_agg<-array(dim=c(7,nblock,nsims))
# for (i in 1:nsims){
#   temp<-data.frame(sim[,,i])
#   temp$Year<-substr(rownames(temp),1,4)
#   temp<-melt(temp)
#   sim_agg[,,i]<-aggregate(temp$value,by=list("variable","Year"),FUN="sum")
# }
