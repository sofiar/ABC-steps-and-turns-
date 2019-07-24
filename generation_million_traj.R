# This script generates the millon simulations need for the experiments. 

## The necessary libraries
library(circular)
library(Rcpp)
library(abc)
library(abctools)
library(HDInterval)
## Necessary functions 
source('pathelements.R')
Rcpp::sourceCpp('cppObs.cpp')
Rcpp::sourceCpp('FastCrw.cpp')

# simulate one trajectory.
# movement parameters:
t_w=5
t_k = 20
dt = 0.5 # time interval for observations
nsteps = 1500 # number of moves performed by the animal
maxt=750
true_sim <- cppFastCRW_exp(t_k,  t_w, nsteps, maxt) # simulate RW using the cpp function
oz <- cppObs(true_sim$x,true_sim$y,true_sim$t,dt) # these are the observed data
nobs=length(oz$st)

# plot trajectory
plot(true_sim$x,true_sim$y)
points(oz$sx,oz$sy,col='red')


############################## simulations ####################################

nsims <- 1e6 # total number of simulated trajectories
nsam <- 25000 # maximum number of real steps in the simulations

# sample from priors
s_w <- runif(nsims,0.1,20)
s_k <- runif(nsims,5,100)

suppressWarnings(warning("as.circular"))
suppressWarnings(warning("conversion.circularmuradians0counter"))
suppressWarnings(warning("rvonmises"))

# creation of summary matrix
Sumsim=matrix(NA,nsims,4)

for( j in 1:nsims){
  
  sim=cppFastCRW_exp(s_k[j], s_w[j], nsam,maxt) # simulate RW using the cpp function. 
  # observe the movement process at time interval dt
  sz <- cppObs(sim$x,sim$y,sim$t,dt) # these are the observed data
  sobs = length(sz$sx)
  
  ps=pathelements(sz$sx,sz$sy)
  # summary statistics
  Sumsim[j,]=c(1/mean(ps$steps),
              log(A1inv(mean.circular(cos(ps$turns)))),
               sd(ps$turns),
               sd(ps$steps))
} 
colnames(Sumsim)=c('w.hat','k.hat','mean.steps','sd.turns','sd.steps',
                   'si', 'sum.steps/T','msd')

sim.params=cbind(s_w,s_k)
colnames(sim.params)=c('w','k')
