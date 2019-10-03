# Author:         Emily S Nightingale
# Institutions:   London Schoool of Hygiene and Tropical Medicine, London, UK
# Date Published: October 2019
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
start <- as.POSIXct(names(df_wide[5]))
start.sts <- c(as.numeric(format(start,"%Y")), as.numeric(format(start,"%m")))
freq <- 12/(month(as.POSIXct(names(df_wide[6]))) - month(as.POSIXct(names(df_wide[5]))))
time <- format(seq(start, by="month", length.out = obsmth))

# Surveillance requires case and population data in matrix form
cases <- t(df_wide[,-1:-4])
pops <- t(pops_wide[,-1:-4])
rownames(pops) <- time
nblock <- ncol(cases)
inc <- cases*1e4/pops

# Using relative population sizes in surveillance makes computation easier
# Observed time period:
popfrac <- pops/rowSums(pops) 
colnames(popfrac) <- rownames(shp@data)

# Population density is a potential covariate
popdens <- pops/shp@data$Shape_Area
logpopdens <- log(popdens)

# Define sts object for use with surveillance functions
stsobj <- sts(observed = cases, start = start.sts, frequency = freq, neighbourhood = nbOrd, map = shp, population = popfrac)
