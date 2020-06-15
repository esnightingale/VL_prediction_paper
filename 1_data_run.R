# Author:         Emily S Nightingale
# Institutions:   London Schoool of Hygiene and Tropical Medicine, London, UK
# Date Published: October 2019
################################################################################
# This script loads the data and sets up the necessary packages/objects for all 
# subsequent analyses. 
# It reads in block-level shapefiles and a long format input file which contains 
# pre-aggregated case counts and population estimates for each block in the 
# region of interest and each month between January 2013 and December 2018.
# Primary output is an sts object, for use with the functions in the 
# surveillance package.

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
library(reshape2)

data.path <- "C:/Users/phpuenig/Dropbox/VL/Data/"
script.path <- "C:/Users/phpuenig/Documents/VL/VL_prediction_paper/"

#---------------------------- DATA SET UP ------------------------------#

setwd(data.path)

# Block map of India
map<-readOGR(dsn = "TAHSILS_CENSUS_2001", layer = "TAHSILS_CENSUS_2001")
#VL region of interest
VL<-map[map$state == "BIHAR"|
           (map$state == "JHARKHAND" & (map$district == "DUMKA" & !(map$admin3 %in% c("Jamtara","Kundahit","Narayanpur_D","Nala"))|
                                          (map$district %in% c("GODDA","PAKAUR","SAHIBGANJ")))),]
# Endemic region of Bihar and Jharkhand
shp<-map[(map$state == "BIHAR" & !(map$district %in% c("GAYA","JAMUI","KAIMUR (BHABUA)","ROHTAS")))|
         (map$state == "JHARKHAND" & (map$district == "DUMKA" & !(map$admin3 %in% c("Jamtara","Kundahit","Narayanpur_D","Nala"))|
                                   (map$district %in% c("GODDA","PAKAUR","SAHIBGANJ")))),]
shp@data$OBJECTID <- rownames(shp@data)

# Define adjacency matrix
nbOrd <- nbOrder(poly2adjmat(shp), maxlag = 7)

# Load input data (month-aggregated case counts for each block)
load("input_sim.RData")
input <- input_sim
# load("./KAMIS/Analysis data/input_jan13dec18.RData")

# Run user-defined functions 
source(paste0(script.path,"2_functions.R"))

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
