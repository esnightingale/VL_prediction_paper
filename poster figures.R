library(hhh4addon)
library(rlist)
library(MASS)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(gridExtra)
# detach(package:plyr)
library(dplyr)
library(rgdal)
library(tmap)

source("C:/Users/phpuenig/Documents/VL/data_setup/1_data_run_addcovariates.R")
load("C:/Users/phpuenig/Dropbox/VL/Monthly prediction/surveillance forecasting/final model selection/Results/agg1325.RData")

setwd("C:/Users/phpuenig/Dropbox/VL/Monthly Prediction/surveillance forecasting/Paper/poster figures")

lshtm <- c("#00AEC7", "#0D5257", "#FE5000", "#00BF6F", "#FFB81C", "#000000", "#A7A8AA", "#FFFFFF")
speak <- c("#0C6B3F", "#F59432", "#232965")
pal <- cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
dist <- readOGR(dsn="C:/Users/phpuenig/Dropbox/VL/RiskMapper/example shapefiles/india", 
                layer="test_IND_district")
dist@data <- setNames(data.frame(dist, do.call(rbind, strsplit(as.character(dist$dot_name),split=":"))),c("hid_id","dot_name","continent","country","state","district"))
dist@data$vl <- 0
dist@data$vl[(dist$state=="BIHAR")|
               (dist$state=="JHARKHAND" & 
                  (dist$district %in% c("DUMKA","GODDA","PAKUR","SAHIBGANJ")))]<- 1
india <- tm_shape(dist) +
  tm_fill(col = "vl", palette = speak[c(1,2)]) + 
  # tm_borders(col=) +
  tm_layout(legend.show = F, frame = F) +
  tm_compass(type = "4star", position = c("right","top"), text.size = 1.5) +
  tm_scale_bar(breaks = c(0, 250,500), position = c("left","bottom"), text.size = 1) 

tiff(filename="India_affected.tiff", height=600, width=600)
india
dev.off()

tmap_leaflet(india) 

tm_shape(dist) + tm_fill(col="state")+tm_layout(legend.show=F)

# State totals over time
# With Year aggregated dataset agg1325 from 7_simulate_2025.R
state_agg <- setNames(data.frame(agg1325, do.call(rbind, strsplit(as.character(agg1325$dot_name),split=":"))),c(names(agg1325),"state","district","block"))
 
# DOESN'T WORK FOR NO REASON

# state_agg <- filter(state_agg, Year<2019) %>%
#              group_by(Year,state) %>%
#              dplyr::summarise(yrtot=sum(cases))

palette(pal)
state_agg <- setNames(aggregate(state_agg$Cases,by=list(state_agg$Year,state_agg$state), FUN="sum"),c("Year","State","Count"))
tiff(filename = "state_totals.tif", height=300, width=400)
ggplot(state_agg[state_agg$Year<2019,], aes(x=Year,y=Count,group=State, fill=State)) + 
  geom_bar(stat="identity") +
  ylab("No. reported cases") +
  xlab("Year") +
  labs(fill = "") +
  theme_classic() +
  theme(legend.position = c(0.8,0.9),
        text = element_text(size=14),
        axis.text = element_text(size=14),
        legend.text = element_text(size = 14)) +
  #       panel.border = element_blank(),
  #       panel.grid.major = element_blank(),
  #       panel.grid.minor = element_blank(),
  #       ) +
  scale_fill_manual(values=lshtm[c(3,1)]) 
dev.off()


blksamp<-sample(1:566,2)
ggplot(data=input2[input$block_id%in%blksamp,],
       aes(interval_start,count,group=block_id)) +
  xlab("Time") +
  ylab("No. reported cases") +
  geom_line() +
  facet_grid(~block_id)

exblks <- input[input$block_id%in%c(199,309),]
tiff(filename = "example_series.tiff", height = 400, width = 700)
ggplot(data=exblks,
       aes(interval_start,count,group=as.factor(block_id), col=as.factor(block_id))) +
  xlab("Time") +
  ylab("No. reported cases") +
  geom_line(lwd=1.1) +
  theme(legend.position = "none") +
  scale_colour_manual(values=lshtm[c(4,3)]) 
dev.off()

pal2 <- magma(6)
palette(lshtm)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# To use for fills, add

plotdata <- agg1325[agg1325$Year%in%c(2013,2016,2018),]
plotdata$Year_f <- rev(as.factor(plotdata$Year))

tiff(filename = "year1318_density.tiff", height = 300, width = 550)
ggplot(data = agg1325[agg1325$Year%in%c(2013,2016,2018),],
       aes(x=Cases,
           #fill=as.factor(Year), 
           colour = as.factor(Year))) + #fill=as.factor(Year), 
       geom_density(lwd = 1.1) + #fill = lshtm[1], col = lshtm[1]
       #scale_x_continuous(trans='log10') +
       scale_fill_manual(values=lshtm[c(3,5,4)]) +
       scale_colour_manual(values=lshtm[c(3,5,4)]) +
       theme_classic() +
       theme(text = element_text(size=14),
             axis.text = element_text(size=14),
             legend.text = element_text(size = 14)) +
       ylab("Density") +
       xlab("Total reported cases per block") +
       labs(fill="", colour = "") +
       theme(legend.position = c(0.8,0.7))
    #   scale_color_gradient(low=lshtm[4], high=lshtm[3]) 
dev.off()

ggplot(data = agg1325[agg1325$Year==2018,],
       aes(x=Cases,stat(density))) + #fill=as.factor(Year), 
       geom_histogram(bins = 40)
ggplot(data = agg1325[agg1325$Year==2016,],
       aes(x=Cases,stat(density))) + #fill=as.factor(Year), 
  geom_histogram(bins = 40)
ggplot(data = agg1325[agg1325$Year==2013,],
       aes(x=Cases,stat(density))) + #fill=as.factor(Year), 
  geom_histogram(bins = 40)

load("C:/Users/phpuenig/Dropbox/VL/Monthly Prediction/surveillance forecasting/exploratory/random model draw/result_table_randmods30.Rdata")
results <- result_table
fitvfirst <- ggplot(results,aes(AIC_train,rps.first, col=k)) +
  theme_classic() +
  theme(legend.position = c(0.8,0.3),
        text = element_text(size=14),
        axis.text = element_text(size=14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) +
  labs(colour="No. of parameters") +
  geom_point() +
 # geom_label_repel(aes(label = Model)) +
  xlab("AIC") + ylab("RPS")  

tiff(filename="AIC_RPS_randmods_poster.tif",height=400,width=400)
fitvfirst + 
  scale_color_gradient(low=lshtm[4], high=lshtm[3]) +
  #scale_color_gradient(low=lshtm[2], high=lshtm[4]) +
  geom_point(cex=3) 
dev.off()
 
png(filename="AIC_RPS_randmods_poster.png",height=400,width=400)
 fitvfirst
dev.off()



