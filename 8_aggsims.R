load("C:/Users/phpuenig/Dropbox/VL/Data/KAMIS/Analysis data/input_jan13dec18.RData")
load("C:/Users/phpuenig/Dropbox/VL/Data/KAMIS/Analysis data/block_agg_jan13dec18.RData")

setwd("C:/Users/phpuenig/Dropbox/VL/Monthly prediction/surveillance forecasting/final model selection")
load("./Results/simulation/sim_2025.RData")


## Need to remove Aurangabad district from sims to match with previous analysis
# blockdata2 <- filter(blockdata,!(blockdata$District%in%c("GAYA","JAMUI","KAIMUR (BHABUA)","ROHTAS")))
# blockdata2[,17:88] <- apply(blockdata2[,17:88],2,FUN=function(x){x[is.na(x)]<-0; return(x)})
# sim <- sim[,which(blockdata2$District!="AURANGABAD"),]

input2 <- filter(input,!(District%in%c("AURANGABAD","GAYA","JAMUI","KAIMUR (BHABUA)","ROHTAS")))
#input2 <- filter(input,!(District%in%c("GAYA","JAMUI","KAIMUR (BHABUA)","ROHTAS")))
input2$count[is.na(input2$count)] <- 0
df_wide <- dcast(input2, OBJECTID+State+District+Block~interval_start, value.var = "count")
df_wide$dot_name <- with(df_wide,paste(State,District,Block, sep=":"))
input2$pop <- as.numeric(input2$pop)
pops_wide <- dcast(input2, OBJECTID+State+District+Block~interval_start, value.var = "pop")
pops_wide <- mutate(pops_wide, dot_name = paste(State,District,Block, sep=":")) %>%
  dplyr::select(dot_name, everything())
rownames(pops_wide) <- pops_wide$OBJECTID
growthR <- pops_wide[,75]/pops_wide[,76]
pops_wide <- cbind(pops_wide,matrix(NA,ncol=fcmth,nrow=nrow(pops_wide)))
for (i in (obsmth+6):ncol(pops_wide)){pops_wide[,i] <- pops_wide[,i-1]*growthR}
names(pops_wide)[-c(1:5)] <- timeall

# Calculate median of sims and bind with observed data
mdfc_blk<-apply(sim, c(1,2), median, names=T)
data1325<-cbind(df_wide,ceiling(t(mdfc_blk))) 

# Long format for ggplot
melt1325<-melt(data1325)
melt1325$Year<-substr(melt1325$variable,1,4)

# Reformat forecasted populations to match
meltpop1325<-melt(pops_wide)
names(meltpop1325)[6:7] <- c("interval_start","pop")
meltpop1325$Year<-substr(meltpop1325$interval_start,1,4)
meltpop1325<-meltpop1325[order(meltpop1325$dot_name,meltpop1325$interval_start),]

# Aggregate monthly pops by mean and monthly counts by sum then bind
aggpop1325<-aggregate(meltpop1325$pop,by=list(meltpop1325$dot_name,meltpop1325$Year),FUN="mean")
aggpop1325<-aggpop1325[order(aggpop1325$Group.1,aggpop1325$Group.2),]
agg1325<-setNames(aggregate(melt1325$value,by=list(melt1325$dot_name,melt1325$OBJECTID,melt1325$Year),FUN="sum"),c("dot_name","OBJECTID","Year","Cases"))
agg1325<-cbind(agg1325[order(agg1325$dot_name,agg1325$Year),],aggpop1325[,3])
names(agg1325)[5]<-"pop"
agg1325 <- mutate(agg1325, inc = Cases*1e4/pop)

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


pop_agg<-melt(pops)
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



# Map predicted values per year - mean and quantiles

melt_sim <- setNames(melt(sim_agg),c("Year","OBJECTID","Sim","Cases"))
quants <- group_by(melt_sim,Year,OBJECTID) %>%
  dplyr::summarise(q25=quantile(Cases,probs=0.25), 
                   q75=quantile(Cases, probs=0.75))
agg1325 <- merge(agg1325,quants,by=c("Year","OBJECTID"), all.x=T)
agg1325 <- mutate(agg1325, q25_inc = q25*1e4/pop, q75_inc = q75*1e4/pop, Year = as.numeric(Year)) %>%
  arrange(OBJECTID,Year)
save(agg1325,file="./Results/agg1325.RData")



pred_agg <- matrix(NA,ncol=7,nrow=502)
colnames(pred_agg) <- paste("med",c(2019:2025), sep="")
for (y in 2019:2025) {
  pred_agg[,(y-2018)] <- agg1325$inc[agg1325$Year==y]
}
pred_agg_upper <- matrix(NA,ncol=7,nrow=502)
colnames(pred_agg_upper) <- paste("upper",c(2019:2025), sep="")
for (y in 2019:2025) {
  pred_agg_upper[,(y-2018)] <- agg1325$q75_inc[agg1325$Year==y]
}
pred_agg_lower <- matrix(NA,ncol=7,nrow=502)
colnames(pred_agg_lower) <- paste("lower",c(2019:2025), sep="")
for (y in 2019:2025) {
  pred_agg_lower[,(y-2018)] <- agg1325$q25_inc[agg1325$Year==y]
}


# Aggregate by district
melt_sim <- setNames(melt(sim),c("Month","OBJECTID","Sim","Cases"))

melt_sim <- merge(df_wide[,1:4], melt_sim, by="OBJECTID") %>%
            arrange(OBJECTID, Month, Sim) 

dist_sims <- group_by(melt_sim, District, Month, Sim) %>%
            dplyr::summarise(dist_tot = sum(Cases))
dist_sims_agg <- aggregate(dist_sims$dist_tot, 
                             by=list(District=dist_sims$District,
                                     Month=dist_sims$Month),
                             FUN="median")

# Path forecast
melt_predmom <- setNames(melt(pathfc),c("Month","OBJECTID","Cases"))
dist_predmom <- merge(df_wide[,1:4], melt_predmom, by="OBJECTID") %>%
                arrange(OBJECTID, Month) %>%
                group_by(District,interval_start) %>%
                dplyr::summarise(Count = sum(count))


dist_obs <- group_by(input2[input2$period_id>60,],
                            District,interval_start) %>%
                   dplyr::summarise(Count = sum(count))
  
dist <- cbind(dist_obs,sim=dist_sims_agg$x)
ggplot(dist, aes(interval_start, Count, group=District)) +
  geom_point() +
  geom_line(aes(interval_start, sim)) +
  facet_wrap(~District, scales="free_y")

ggsave(filename = "./figures/district_10sims.png", height=5, width=7)
