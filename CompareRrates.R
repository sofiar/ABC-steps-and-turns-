# This script performes the Relative scale of observations analysis. 

#set the R values
Rates=c(0.06,0.2,0.5, 0.8,  1,  1.3,  2,  2.4, 3,  3.5,  4,  4.5, 4.8,5,
7,10,12,11,8,8.7,9.5)
Mean.time=dt/Rates
nn=length(Rates)

# Set parameters
Alphas=1/Mean.time 
kappas=c(10,20,40,60,70)

# Creation of error matrices
Error.lc.w=list(numeric(length (kappas)))
Error.nn.w=list(numeric(length (kappas)))
Error.rec.w=list(numeric(length (kappas)))

Error.lc.k=list(numeric(length (kappas)))
Error.nn.k=list(numeric(length (kappas)))
Error.rec.k=list(numeric(length (kappas)))

for (s in 1:length(kappas))
{
  Error.lc.w[[s]]= matrix(NA,50,nn)
  Error.nn.w[[s]]= matrix(NA,50,nn)
  Error.rec.w[[s]]= matrix(NA,50,nn)
  
  Error.lc.k[[s]]= matrix(NA,50,nn)
  Error.nn.k[[s]]= matrix(NA,50,nn)
  Error.rec.k[[s]]= matrix(NA,50,nn)
  
}

# Loop
for (i in 1:nn) # for every rate R
{
  for (jj in 1:length(kappas)) # for every kappa value
  {
    for (k in 1:50) # simulation of 50 trajectories  
    {
      tsim <- cppFastCRW_exp(kappas[jj],Alphas[i], nsam, maxt) # simulate RW using the cpp function
      ozz <- cppObs(tsim$x,tsim$y,tsim$t,dt) # compute the observation  
      
      #compute summaries for the true one
      pee=pathelements(ozz$sx,ozz$sy)
      
      pe.ws=1/mean(pee$steps)
      pe.ks=log(A1inv(mean.circular(cos(pee$turns))))
      
      sss=c(pe.ws,pe.ks,
            sd(pee$turns),
            sd(pee$steps))

      # ABC inference   
      fit.nn=abc(target=sss[quienes],param=sim.params,
                 sumstat = Sumsim[,quienes],tol=.005,method='neuralnet')
      
      fit.lc=abc(target=sss[quienes],param=sim.params,
                 sumstat = Sumsim[,quienes],tol=.005,method='ridge')
      
      # save Error
      Error.nn.w[[jj]][k,i]=sqrt((median(fit.nn$adj.values[,1])-Alphas[i])^2)
      Error.lc.w[[jj]][k,i]=sqrt((median(fit.lc$adj.values[,1])-Alphas[i])^2)
      Error.rec.w[[jj]][k,i]=sqrt((median(fit.nn$unadj.values[,1])-Alphas[i])^2)

      Error.nn.k[[jj]][k,i]=sqrt((median(fit.nn$adj.values[,2])-kappas[jj])^2)
      Error.lc.k[[jj]][k,i]=sqrt((median(fit.lc$adj.values[,2])-kappas[jj])^2)
      Error.rec.k[[jj]][k,i]=sqrt((median(fit.nn$unadj.values[,2])-kappas[jj])^2)
    }
  }
}




