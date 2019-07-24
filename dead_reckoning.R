#### Dead reckoning 

## The necessary libraries
library(TrackReconstruction)
library(BayesianAnimalTracker)
library(tidyverse)
library(dplyr)
library(SoDA)

###############################################################################
################ 1.Upload DailyDiary data and tidy it #########################
###############################################################################

Filename='Filename_split_31.txt'

demo.data=read_delim(Filename, 
                     "\t", escape_double = FALSE, trim_ws = TRUE)
# creation of DateTime
demo.data$Time=paste(sprintf("%02d",demo.data$Hours),
                     sprintf("%02d",demo.data$Minutes),
                     sprintf("%02d",demo.data$Seconds),sep=':')

demo.data$DateTime <- as.POSIXct(strptime(paste(demo.data$Date, demo.data$Time), "%d/%m/%Y %H:%M:%S"))

# change some Features names and orientations
demo.data=demo.data %>% mutate(AccSurge=Acc_x,AccSway=-Acc_y,AccHeave=-Acc_z,
                               MagSurge=Mag_x,MagSway=-Mag_y,MagHeave=-Mag_z)

###############################################################################
######################## 2.Estimation max speed ###############################
###############################################################################

# upload gps data
gps <- read_csv("13095_190422.csv", 
                col_types = cols(Date = col_date(format = "%m/%d/%Y"),Time = col_time(format = "%H:%M:%S")))

# Eliminate outliers
gps=gps %>% mutate(Latitude=replace(Latitude,Latitude<=-2000,NA)) %>% 
  mutate(Longitude=replace(Longitude,Longitude<=-2000,NA)) 

# Convert datetime data from GPS into DateTime POSIXct format
gps$DateTime <-as.POSIXct(strptime(paste(gps$Date, as.character(gps$Time)),"%Y-%m-%d %H:%M:%S"))

# Matching time of the GPS and DD:
gpsdata <- gps[gps$DateTime %in% demo.data$DateTime, ]

# Format GPS data using the GPStable function
gpsformat <- GPStable(gpsdata)

# obtain coordinates in meters
mmm=geoXY(gpsdata$Latitude, gpsdata$Longitude,unit = 1)

# speed estimation
n.points=length(gpsdata$DateTime)
t.b=difftime(gpsdata$DateTime[2:n.points],gpsdata$DateTime[1:n.points-1], units = "secs")
vel=(pathelements(mmm[,1],mmm[,2])$steps)/as.numeric(t.b) # meters per second
speed=max(vel)

###############################################################################
######################## 3.Estimation of DR path ##############################
###############################################################################

# Estimation of the betas:
betas <- Standardize(1, 1, 1, -1, -1, 1, min(demo.data$MagSurge), 
                     max(demo.data$MagSurge), min(demo.data$MagHeave),
                     max(demo.data$MagHeave), min(demo.data$MagSway), 
                     max(demo.data$MagSway), min(demo.data$AccSurge), 
                     max(demo.data$AccSurge), min(demo.data$AccHeave), 
                     max(demo.data$AccHeave), min(demo.data$AccSway), 
                     max(demo.data$AccSway))

# Estimation of the Declination and Inclination of magnetic field.
# extracted from: http://www.geomag.bgs.ac.uk/data_service/models_compass/wmm_calc.html
D_MF <- 5.403 #Main Field Declination (in degrees east)
I_MF <- -41.477 #Main Field Inclination (in degrees east)

# Convert Date and Time data from DD into DateTime POSIXct format:
Data=demo.data %>% dplyr::select(MagSurge,MagSway,MagHeave,AccSurge,AccSway,
                          AccHeave,Date,Time)
# Generation of DeadReckoning patch using the function from "TrackReconstruction" package
DRoutput <- DeadReckoning(Data, betas, decinc = c(D_MF, I_MF), 
                          Hz = 40, RmL = 2, DepthHz = 1, SpdCalc = 3, MaxSpd = speed)
# change to POSIXct
DRoutput$DateTime=as.POSIXct(strptime(DRoutput$DateTime,"%d/%m/%Y %H:%M:%S"))

###############################################################################
################## 4.GPS correction of DeadReckoning path ######################
###############################################################################

# Define starting and ending points from GPS data:
K.demo <- nrow(gpsformat)
DRstart <- min(which(DRoutput$DateTime==gpsformat$DateTime[1]))
DRend <- max(which(DRoutput$DateTime==gpsformat$DateTime[K.demo]))
# Thin the data 
DRworking <- DRoutput[c(DRstart:DRend)[c(DRstart:DRend)%%40==1], ]
# Calculate the northing in km for GPS data and for DR paths:
T.demo <- nrow(gpsformat)
GPSnorthing <- c(cumsum(gpsformat$DistanceKm[-1]*cos(gpsformat$BearingRad[-T.demo])))
GPSeasting <- c(cumsum(gpsformat$DistanceKm[-1]*sin(gpsformat$BearingRad[-T.demo])))
#Original unit of DR is in meters, so we convert it in kilometers:
DRnorthing <- (DRworking$Ydim - DRworking$Ydim[1])/1000 
DReasting <- (DRworking$Xdim - DRworking$Xdim[1])/1000 
# Calculation for northing and easting:
nlist <- as.dataList(DRnorthing, GPSnorthing, Ytime = format(gpsformat$DateTime, "%d-%b-%Y %H:%M:%S"), Xtime = format(DRworking$DateTime, "%d-%b-%Y %H:%M:%S"), s2G=0.0001, timeUnit = 40*60, betaOrder=1)
npost <- BMAnimalTrack(nlist, BMControl(print=TRUE, returnParam=TRUE))
elist <- as.dataList(DReasting, GPSeasting, Ytime = format(gpsformat$DateTime, "%d-%b-%Y %H:%M:%S"), Xtime = format(DRworking$DateTime, "%d-%b-%Y %H:%M:%S"), s2G=0.0001, timeUnit = 40*60, betaOrder = 1)
epost <- BMAnimalTrack(elist, BMControl(print = TRUE, returnParam = TRUE))

# Final results
cPathInKM <- cbind(epost$etaMar[,1], npost$etaMar[,1]) 
cPathInDeg <- KMToDeg(cPathInKM, gpsformat [1, c(3, 2)])
