source("C:/Users/phpuenig/Dropbox/VL/Monthly Prediction/surveillance forecasting/final model selection/data_run_update2604.R")

setwd("C:/Users/phpuenig/Dropbox/VL/Monthly Prediction/surveillance forecasting/final model selection/Results")
load("pred3_final_osafirst.RData")
load("pred4_final_osafirst.RData")
load("ncms_4ahd.RData")

1-ncms.4ahd[[1]]
# 0.9475432 0.8568583 0.6839309

fc4mth<-pred4.final[[1]]
rownames(fc4mth)<-timeall[52:72]
fc4mth<-melt(fc4mth)
cases_train<-melt(cases[1:48,])
obspred_4ahd<-rbind(cases_train,fc4mth)

cases_test<-melt(cases[52:72,])

hist(fc4mth$value,freq=T, col=rgb(1,0,0,0.3), breaks=30)
hist(cases_test$value, freq=T, add=T,col=rgb(1,0,0,0.3), breaks=30)

summary(as.numeric(cases_test$value))
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   0.000   0.000   0.815   1.000  32.000  
summary(fc4mth$value)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0118  0.1316  0.4147  0.8189  1.0299 11.2871  

cases_test$error<-fc4mth$value-cases_test$value


length(which(abs(cases_test$error)<3))/nrow(cases_test)
# 0.9689812
length(which(abs(cases_test$error)<2))/nrow(cases_test)
# 0.92952

ggplot(cases_test,aes(value,error))+geom_point(col="indianred")+geom_hline(yintercept = 0)

png(file="./figures/report figures/fc_validation_2018_2.png", height=600, width=800)
plot(NULL,xlim=c(0,max(cases18_2_melt$value)), ylim=c(min(cases18_2_melt$error_mean),max(cases18_2_melt$error_mean)), xlab="Observed",ylab="Forecast error", main="Error of forecasts from mean/median of 1,000 simulations")
abline(h=0,col="grey",lty="dashed")
polygon(x = c(0:14, rev(0:14)), y = c(rep(-2,15),rep(2,15)), col =  adjustcolor("snow4", alpha.f = 0.1), border = NA)
points(cases18_2_melt$value,cases18_2_melt$error_mean,col="indianred")
points(cases18_2_melt$value,cases18_2_melt$error_med,col="dodgerblue")
lines(lowess(cases18_2_melt$value,cases18_2_melt$error_med), col="dodgerblue")
lines(lowess(cases18_2_melt$value,cases18_2_melt$error_mean), col="indianred")
legend(2,-7,c("Median","Mean"),fill = c("dodgerblue","indianred"), border = "white", bty="n")
axis(2,at=c(-10,-5,-2,0,2))
dev.off()



# All 2018 data
fc18<-fc_melt_ci[fc_melt_ci$year=="2018",]
cases18<-cases[61:72,]

hist(fc18$median,freq=T, col=rgb(1,0,0,0.3), breaks=15)
hist(cases18, freq=T, add=T,col=rgb(1,0,0,0.3), breaks=15)
hist(fc18$mean, freq=T, add=T,col=rgb(1,0,0,0.3), breaks=15)

summary(as.numeric(cases18))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.0000  0.0000  0.7386  1.0000 18.0000 
summary(fc18$median)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.0000  0.0000  0.5107  1.0000 18.0000 
summary(fc18$mean)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.0000  0.1390  0.6551  0.9443 18.0000 

cases18_melt<-melt(cases18)
error_med<-fc18$median-cases18_melt$value
error_mean<-fc18$mean-cases18_melt$value
cases18_melt<-cbind(cases18_melt,error_med,error_mean)

png(file="./figures/fc_validation_2018.png", height=600, width=800)
plot(cases18_melt$value,cases18_melt$error_med,col="dodgerblue", xlab="Observed",ylab="Forecast error (1000 sims)")
points(cases18_melt$value,cases18_melt$error_mean,col="indianred")
abline(h=0,col="grey",lty="dashed")
lines(lowess(cases18_melt$value,cases18_melt$error_med), col="dodgerblue")
lines(lowess(cases18_melt$value,cases18_melt$error_mean), col="indianred")
legend(2,-7,c("Median","Mean"),fill = c("dodgerblue","indianred"))
dev.off()

length(which(abs(cases18_melt$error_med)<5))/nrow(cases18_melt)
# 0.9573999
length(which(abs(cases18_melt$error_med)<3))/nrow(cases18_melt)
# 0.8915479
length(which(abs(cases18_melt$error_mean)<5))/nrow(cases18_melt)
# 0.9614732


