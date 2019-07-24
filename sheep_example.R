# This script performs the analysis for sheep data 

## The necessary libraries
library(tidyverse)
library(circular)
library(Rcpp)
library(abc)
library(abctools)
library(HDInterval)
library(gridExtra)
library(grid)
library(cowplot)
#############################################################################
### In this script we use the Dead_Reckoning output (corrected with gps)####
################## as the "real" trajectory ################################
#############################################################################

#############################################################################
##################### 1. Calculation of the headings ########################
#############################################################################

# dead reckoning
source('dead_reckoning.R') 
data=cPathInDeg[, ]


# Select part of the trajectory to obtain the same behaviour 

n.length=length(data$Longitude)-21000
x_vals=data$Longitude[1:n.length]
y_vals=data$Latitude[1:n.length]

range.ti=range(DRworking$DateTime[1:n.length])
difftime(range.ti[2],range.ti[1])

#plot the path
plot(x_vals,y_vals,asp=1,type='l')

n.tot=length(x_vals)
adj=x_vals[2:n.tot]-x_vals[1:(n.tot-1)]
op=y_vals[2:n.tot]-y_vals[1:(n.tot-1)]

# calculation of the headings
si<-sign(op)
si[si==0]<-1  #corrects for sign(0) == 0
headings= si*(adj<0)*pi+atan(op/adj)

#########################################################################
##### 2. Calculation of turning points: From Potts et.al 2018 ###########
#########################################################################

# Set window size and the threshold angle in degrees 
window_size <- 15
thresh <- 5

# Calculate sin and cosine of headings, together with x- and y-values 
# (assuming a constant speed)
cosines <- cos(headings)
sines <- sin(headings)
x_vals <- cumsum(cosines)
y_vals <- cumsum(sines)

## Change in scale
x_vals=(x_vals-x_vals[1])/100
y_vals=(y_vals-y_vals[1])/100


sd_length <- length(cosines) - window_size
ave_cos_array <- cosines[1:sd_length]/window_size
ave_sin_array <- sines[1:sd_length]/window_size
for(counter in 1:(window_size-1))
{
  # Shift cosine and sin arrays back by 1
  cosines <- cosines[-1]
  sines <- sines[-1]
  ave_cos_array <- ave_cos_array + cosines[1:sd_length]/window_size
  ave_sin_array <- ave_sin_array + sines[1:sd_length]/window_size
}
circ_sd <- (-2)*log(sqrt((ave_sin_array)^2 + (ave_cos_array)^2))
ave_circ_sd <- sum(circ_sd) / length(circ_sd)

# Find candidate turning points by looking for spikes in SCSD
turning <- 0
chgpt_array <- c(1)
for(counter in 1:length(circ_sd))
{
  if((circ_sd[counter] > ave_circ_sd)&(turning == 0))
  {
    # Started turning.  Note turning point
    start_chgpt <- counter 
    turning <- 1
  }
  if((circ_sd[counter] < ave_circ_sd)&(turning == 1))
  {
    # Turning point is the mean of the start turning-point and the end shifted to the right by resolution/2 
    # to account for averaging being done forwards in time
    chgpt <- floor((start_chgpt + counter)/2 + window_size/2)
    chgpt_array <- c(chgpt_array,chgpt)
    turning <- 0
  }
}
chgpt_array <- c(chgpt_array, length(cosines))

# Post processing of changepoints to remove any where the switch in direction is less than thresh_angle
new_chgpt_array <- vector(mode="numeric", length=0)
prev_heading <- atan2(y_vals[chgpt_array[2]]-y_vals[chgpt_array[1]],
                      x_vals[chgpt_array[2]]-x_vals[chgpt_array[1]])
for(counter in 2:(length(chgpt_array)-1))
{
  next_heading <- atan2(y_vals[chgpt_array[counter+1]]-y_vals[chgpt_array[counter]],
                        x_vals[chgpt_array[counter+1]]-x_vals[chgpt_array[counter]])                              
  if(abs(next_heading - prev_heading) > pi)
  {
    turn_angle <- pi*2-abs(next_heading - prev_heading)
  }
  else
  {
    turn_angle <- abs(next_heading - prev_heading)
  }
  if(turn_angle*180/pi > thresh)
  {
    new_chgpt_array <- c(new_chgpt_array, chgpt_array[counter])
    prev_heading = next_heading  
  }
}


######################################################################
############# 3. Estimation of the parameters usin rstan and #########
#############   calcuation of the Summary statistics  ################
######################################################################
library(rstan)
#1-Upload packges and some functions 
source('pathelements.R')
Rcpp::sourceCpp('cppObs.cpp')

pee=pathelements(x.tpoints,y.tpoints)

# data for stan
data.stan <- list(Num=length(pee$steps), steps=pee$steps,angles=pee$turns)
# fit with stan 
fit <- stan(file="model.stan", data=data.stan,chains=3)

# outs and results
traceplot(fit,pars=c("lambda","kappa"))
outs <- rstan::extract(fit, permuted = TRUE) # return a list of arrays 

# Summaries calculated at the "observed path" every dt
stepp=pathelements(x_vals,y_vals)$steps[1]
t_vals=seq(from=1,by=stepp,length.out=length(x_vals))
dt=.5
obs.dt=cppObs(x_vals,y_vals,t_vals,dt)[-1,] # eliminate the first row (reference to (0,0))
p.obs=pathelements(obs.dt$sx,obs.dt$sy)

summaries.true=c(1/mean(p.obs$steps),
                 log(A1inv(mean.circular(cos(p.obs$turns[2:length(p.obs$turns)])))),
                 sd(p.obs$turns[2:length(p.obs$turns)]),
                 sd(p.obs$steps))


######################################################################
########################### 4. ABC Inference #########################
######################################################################

########################
#### Neural Network ####
########################

set.seed(16007)
fit.nn=abc(target=summaries.true,param=sim.params,
           sumstat = Sumsim,tol=.001,method='neuralnet')

###################
#### Loclinear ####
###################

fit.ln=abc(target=summaries.true,param=sim.params,
           sumstat = Sumsim,tol=.001,method='ridge')

##########################
#### Simple Rejection ####
##########################

fit.rec=abc(target=summaries.true,param=sim.params,
            sumstat = Sumsim,tol=.001,method='rejection')

