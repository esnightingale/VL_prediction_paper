# Author:         Emily S Nightingale
# Institutions:   London Schoool of Hygiene and Tropical Medicine, London, UK
# Date Published: XX July 2019
################################################################################
# This script loads the data and sets up the necessary packages/objects for all 
# subsequent analyses. 
# It reads in block-level shapefiles and a long format input file which contains 
# pre-aggregated case counts and population estimates for each block in the 
# region of interest and each month between January 2013 and December 2018.
# Primary outputs are three sts objects for use with the functions in the 
# surveillance package:
#  a) stsobj1 - excludes all zero-count blocks and blocks which become 
#               disconnected by exclusion of zero-count blocks. This is used
#               for preliminary analyses in "a2_correlation.R" to investigate 
#               the extent of spatial and temporal correlation.
#  b) stsobj  - Primary object used in the main analysis.
#  c) stsobj2 - Includes dummy rows to allow forecasting ahead of the observed 
#               data.
################################################################################
################################################################################

# HHH4addon package is still in development and must be downloaded from github:
# install_github("jbracher/hhh4addon")

# Packages required for all analyses/output:
library(plyr)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(maptools)
library(data.table)
library(rgdal)
library(spdep)
library(pracma)
library(surveillance)
library(stringr)
library(png)
library(abind)
library(hhh4addon)
library(lattice)
library(rlist)
library(MASS)
library(ggrepel)
library(cowplot)
library(gridExtra)
library(fanplot)
library(reshape2)

setwd("C:/Users/phpuenig/Dropbox/VL/Data")

#---------------------------- DATA SET UP ------------------------------#

# Load maps
map<-readOGR(dsn = "TAHSILS_CENSUS_2001",layer = "TAHSILS_CENSUS_2001")
shp<-map[(map$state == "BIHAR" & !(map$district %in% c("GAYA","JAMUI","KAIMUR (BHABUA)","ROHTAS")))|
         (map$state == "JHARKHAND" & (map$district == "DUMKA" & !(map$admin3 %in% c("Jamtara","Kundahit","Narayanpur_D","Nala"))|
                                   (map$district %in% c("GODDA","PAKAUR","SAHIBGANJ")))),]
shp@data$OBJECTID <- rownames(shp@data)

# Define adjacency matrix
nbOrd <- nbOrder(poly2adjmat(shp), maxlag = 7)

# Load input data (month-aggregated case counts for each block)
load("./KAMIS/Analysis data/archive/input_jan13dec18_1406.RData")

# No. months to be forecasted (up to 2025)
fcmth <- 84

# Run user-defined functions 
source("C:/Users/phpuenig/Documents/VL/surveillance_modelling/2_functions.R")

#-----------------------------------------------------------------------#

# Remove four non-endemic districts and assume remaining missing numbers are 0 
input2 <- input[!(input$District %in% c("GAYA","JAMUI","KAIMUR (BHABUA)","ROHTAS")),]
input2$count[is.na(input2$count)] <- 0

# Double check no observations outside region of interest
input2<-input2[input2$State %in% c("BIHAR", "JHARKHAND"),]

# Reshape case and population data. Assign IDs to match shapefiles.
df_wide <- dcast(input2, OBJECTID+State+District+Block~interval_start, value.var = "count")
rownames(df_wide) <- df_wide$OBJECTID
nblock <- nrow(df_wide)
pops_wide <- dcast(input2, OBJECTID+State+District+Block~interval_start, value.var = "pop")
rownames(pops_wide) <- pops_wide$OBJECTID

# Define time period
obsmth <- ncol(df_wide)-4
totmth <- obsmth + fcmth
tp <- c((obsmth+1):totmth)
start <- as.POSIXct(names(df_wide[5]))
start.sts <- c(as.numeric(format(start,"%Y")), as.numeric(format(start,"%m")))
freq <- 12/(month(as.POSIXct(names(df_wide[6]))) - month(as.POSIXct(names(df_wide[5]))))
time <- format(seq(start, by="month", length.out = max(tp)))
timeall <- format(seq(start, by="month", length.out = totmth))

# Surveillance requires case and population data in matrix form
cases <- t(df_wide[,-1:-4])
pops <- t(pops_wide[,-1:-4])
nblock <- ncol(cases)
inc <- cases*1e4/pops[1:obsmth,]

#length(which(is.na(colSums(cases))))

# Add projected populations to pops
growthR <- pops[obsmth,]/pops[(obsmth-1),]
pops <- rbind(pops, matrix(NA, nrow=fcmth, ncol=ncol(pops)))
for (i in (obsmth+1):totmth){pops[i,] <- pops[i-1,]*growthR}
rownames(pops) <- time

# Using relative population sizes in surveillance makes computation easier
# Observed time period:
popfrac <- pops[1:obsmth,]/rowSums(pops[1:obsmth,]) 
colnames(popfrac) <- rownames(shp@data)
# Including time period to be forecasted:
popfrac2 <- pops/rowSums(pops) 
colnames(popfrac2) <- rownames(df_wide)
# Add dummy rows in which to write forecasted case counts
addrows <- matrix(nrow=fcmth,ncol=nblock)
cases2 <- rbind(cases,addrows)
rownames(cases2) <- timeall

# Population density is a potential covariate
popdens <- pops[1:obsmth,]/shp@data$Shape_Area
logpopdens <- log(popdens)
input2$popdens <- melt(popdens)[,3]
input2$logpopdens <- log(input2$popdens)
input2$popdens_c <- (input2$popdens-mean(input2$popdens,na.rm=T))/sqrt(var(input2$popdens))

# Define sts objects for use with surveillance functions
# Primary:
stsobj <- sts(observed = cases, start = start.sts, frequency = freq, neighbourhood = nbOrd, map = shp, population = popfrac)
# Secondary (with empty rows for predicted time points)
stsobj2 <- sts(observed = cases2, start = start.sts, frequency = freq, neighbourhood = nbOrd, map = shp, population = popfrac2)

# Third stsobj for preliminary correlation analysis (excluding blocks with no cases)
cases1 <- cases[,-which(colSums(cases) == 0)] #exclude zero count blocks
shp1 <- shp[shp$OBJECTID %in% colnames(cases1),] 
nb.r <- poly2nb(shp1)
Amat <- nb2mat(nb.r, style = "B", zero.policy = T)

# Find any blocks that are unconnected in the non-zero map
island_index <- which(sapply(nb.r, FUN = function(x)(x[1] == 0)))
island_index

# # Visualize islands
# clr <- rep('white', length(shp1))
# clr[island_index] <- 'red'
# plot(shp1, col = clr)
# # Find the nearest neighbor for each island (nearest by centroid)
# centroids <- coordinates(shp1)
# nb.nearest <- knn2nb(knearneigh(centroids, k = 1))
# nearest_index <- unlist(nb.nearest[island_index])
# # Create adjacency matrix to account for islands (needs to be reciprocal)
# Amat[cbind(island_index, nearest_index)] <- 1
# Amat[cbind(nearest_index, island_index)] <- 1

# Additionally exclude islands created by removing zero count blocks (i.e. zero cases in neighbourhood)
cases1 <- cases1[,-island_index]
popfrac1 <- popfrac[, which(colnames(cases) %in% colnames(cases1))]
shp1 <- shp1[-island_index,]
plot(shp1)
nbOrd1 <- nbOrder(poly2adjmat(shp1), maxlag = 7)

stsobj1 <- sts(observed = cases1, start = start.sts, frequency = freq, neighbourhood = nbOrd1, map = shp1, population = popfrac1)

