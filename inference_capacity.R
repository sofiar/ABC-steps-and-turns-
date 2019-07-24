# This script performes the analysis of the inference capacity of the 
# ABC methods

# Cross validation analysis
nCv=100 # set Nrep
tols=c(0.001,.005,.01,.1) # set threshold values
a=seq(1,length(sim.params[,1]))[sim.params[,1]<25 & sim.params[,2]<70] # Limit the possible values of the true parameters
cvsamp <- sample(a, nCv) #sample parameters
cv.param.true=sim.params[cvsamp,]
cv.fits.nn=list()
cv.fits.rj=list()
cv.fits.lc=list()

for (t in 1:length(tols))
{
  cv.fits.nn[[t]]=matrix(NA,ncol=2,nrow=nCv)
  cv.fits.rj[[t]]=matrix(NA,ncol=2,nrow=nCv)
  cv.fits.lc[[t]]=matrix(NA,ncol=2,nrow=nCv)
  
  for (m in 1:nCv)  
  {
  ### NN
  fit.nn.cv=abc(target=Sumsim[cvsamp[m],],param=sim.params[-cvsamp[m],],
                sumstat = Sumsim[-cvsamp[m],],tol=tols[t],method='neuralnet')
  cv.fits.nn[[t]][m,]=c(median(fit.nn.cv$adj.values[,1]),median(fit.nn.cv$adj.values[,2]))
  
  ### RJ
  cv.fits.rj[[t]][m,]=c(median(fit.nn.cv$unadj.values[,1]),median(fit.nn.cv$unadj.values[,2]))
  
  ### LC
  fit.lc.cv=abc(target=Sumsim[cvsamp[m],],param=sim.params[-cvsamp[m],],
                sumstat = Sumsim[-cvsamp[m],],tol=tols[t],method='loclinear')
  cv.fits.lc[[t]][m,]=c(median(fit.lc.cv$adj.values[,1]),median(fit.lc.cv$adj.values[,2]))
  }
}
  
#### Errors and MV ####
Error.cv=list()
sd.cv=list()

## tol 0.001
Error.cv[[1]]=matrix(NA,ncol=2,nrow=3)
rownames(Error.cv[[1]])=c('Neural Network','Linear','Rejection')
colnames(Error.cv[[1]])=c('lambda','kappa')

sd.cv[[1]]=matrix(NA,ncol=2,nrow=3)
rownames(sd.cv[[1]])=c('Neural Network','Linear','Rejection')
colnames(sd.cv[[1]])=c('lambda','kappa')

Error.cv[[1]][1,1]=sqrt(sum((cv.fits.nn[[1]][,1]-cv.param.true[,1])^2)/(nCv))
Error.cv[[1]][1,2]=sqrt(sum((cv.fits.nn[[1]][,2]-cv.param.true[,2])^2)/(nCv))

Error.cv[[1]][2,1]=sqrt(sum((cv.fits.lc[[1]][,1]-cv.param.true[,1])^2)/(nCv))
Error.cv[[1]][2,2]=sqrt(sum((cv.fits.lc[[1]][,2]-cv.param.true[,2])^2)/(nCv))

Error.cv[[1]][3,1]=sqrt(sum((cv.fits.rj[[1]][,1]-cv.param.true[,1])^2)/(nCv))
Error.cv[[1]][3,2]=sqrt(sum((cv.fits.rj[[1]][,2]-cv.param.true[,2])^2)/(nCv))

sd.cv[[1]][1,1]= mean(abs(cv.fits.nn[[1]][,1]-cv.param.true[,1])/(cv.param.true[,1]))
sd.cv[[1]][1,2]= mean(abs(cv.fits.nn[[1]][,2]-cv.param.true[,2])/(cv.param.true[,2]))
  
sd.cv[[1]][2,1]= mean(abs(cv.fits.lc[[1]][,1]-cv.param.true[,1])/(cv.param.true[,1]))
sd.cv[[1]][2,2]= mean(abs(cv.fits.lc[[1]][,2]-cv.param.true[,2])/(cv.param.true[,2]))
  
sd.cv[[1]][3,1]= mean(abs(cv.fits.rj[[1]][,1]-cv.param.true[,1])/(cv.param.true[,1]))
sd.cv[[1]][3,2]= mean(abs(cv.fits.rj[[1]][,2]-cv.param.true[,2])/(cv.param.true[,2]))
  
## tol 0.005
Error.cv[[2]]=matrix(NA,ncol=2,nrow=3)
rownames(Error.cv[[2]])=c('Neural Network','Linear','Rejection')
colnames(Error.cv[[2]])=c('lambda','kappa')

sd.cv[[2]]=matrix(NA,ncol=2,nrow=3)
rownames(sd.cv[[2]])=c('Neural Network','Linear','Rejection')
colnames(sd.cv[[2]])=c('lambda','kappa')

Error.cv[[2]][1,1]=sqrt(sum((cv.fits.nn[[2]][,1]-cv.param.true[,1])^2)/(nCv))
Error.cv[[2]][1,2]=sqrt(sum((cv.fits.nn[[2]][,2]-cv.param.true[,2])^2)/(nCv))

Error.cv[[2]][2,1]=sqrt(sum((cv.fits.lc[[2]][,1]-cv.param.true[,1])^2)/(nCv))
Error.cv[[2]][2,2]=sqrt(sum((cv.fits.lc[[2]][,2]-cv.param.true[,2])^2)/(nCv))

Error.cv[[2]][3,1]=sqrt(sum((cv.fits.rj[[2]][,1]-cv.param.true[,1])^2)/(nCv))
Error.cv[[2]][3,2]=sqrt(sum((cv.fits.rj[[2]][,2]-cv.param.true[,2])^2)/(nCv))

sd.cv[[2]][1,1]=mean(abs(cv.fits.nn[[2]][,1]-cv.param.true[,1])/(cv.param.true[,1]))
sd.cv[[2]][1,2]=mean(abs(cv.fits.nn[[2]][,2]-cv.param.true[,2])/(cv.param.true[,2]))

sd.cv[[2]][2,1]=mean(abs(cv.fits.lc[[2]][,1]-cv.param.true[,1])/(cv.param.true[,1]))
sd.cv[[2]][2,2]=mean(abs(cv.fits.lc[[2]][,2]-cv.param.true[,2])/(cv.param.true[,2]))
  
sd.cv[[2]][3,1]=mean(abs(cv.fits.rj[[2]][,1]-cv.param.true[,1])/(cv.param.true[,1]))
sd.cv[[2]][3,2]=mean(abs(cv.fits.rj[[2]][,2]-cv.param.true[,2])/(cv.param.true[,2]))
  

## tol 0.01
Error.cv[[3]]=matrix(NA,ncol=2,nrow=3)
rownames(Error.cv[[3]])=c('Neural Network','Linear','Rejection')
colnames(Error.cv[[3]])=c('lambda','kappa')

sd.cv[[3]]=matrix(NA,ncol=2,nrow=3)
rownames(sd.cv[[3]])=c('Neural Network','Linear','Rejection')
colnames(sd.cv[[3]])=c('lambda','kappa')

Error.cv[[3]][1,1]=sqrt(sum((cv.fits.nn[[3]][,1]-cv.param.true[,1])^2)/(nCv))
Error.cv[[3]][1,2]=sqrt(sum((cv.fits.nn[[3]][,2]-cv.param.true[,2])^2)/(nCv))

Error.cv[[3]][2,1]=sqrt(sum((cv.fits.lc[[3]][,1]-cv.param.true[,1])^2)/(nCv))
Error.cv[[3]][2,2]=sqrt(sum((cv.fits.lc[[3]][,2]-cv.param.true[,2])^2)/(nCv))

Error.cv[[3]][3,1]=sqrt(sum((cv.fits.rj[[3]][,1]-cv.param.true[,1])^2)/(nCv))
Error.cv[[3]][3,2]=sqrt(sum((cv.fits.rj[[3]][,2]-cv.param.true[,2])^2)/(nCv))

sd.cv[[3]][1,1]=mean(abs(cv.fits.nn[[3]][,1]-cv.param.true[,1])/(cv.param.true[,1]))
sd.cv[[3]][1,2]=mean(abs(cv.fits.nn[[3]][,2]-cv.param.true[,2])/(cv.param.true[,2]))

sd.cv[[3]][2,1]=mean(abs(cv.fits.lc[[3]][,1]-cv.param.true[,1])/(cv.param.true[,1]))
sd.cv[[3]][2,2]=mean(abs(cv.fits.lc[[3]][,2]-cv.param.true[,2])/(cv.param.true[,2]))

sd.cv[[3]][3,1]=mean(abs(cv.fits.rj[[3]][,1]-cv.param.true[,1])/(cv.param.true[,1]))
sd.cv[[3]][3,2]=mean(abs(cv.fits.rj[[3]][,2]-cv.param.true[,2])/(cv.param.true[,2]))

## tol 0.1
Error.cv[[4]]=matrix(NA,ncol=2,nrow=3)
rownames(Error.cv[[4]])=c('Neural Network','Linear','Rejection')
colnames(Error.cv[[4]])=c('lambda','kappa')

sd.cv[[4]]=matrix(NA,ncol=2,nrow=3)
rownames(sd.cv[[4]])=c('Neural Network','Linear','Rejection')
colnames(sd.cv[[4]])=c('lambda','kappa')

Error.cv[[4]][1,1]=sqrt(sum((cv.fits.nn[[4]][,1]-cv.param.true[,1])^2)/(nCv))
Error.cv[[4]][1,2]=sqrt(sum((cv.fits.nn[[4]][,2]-cv.param.true[,2])^2)/(nCv))

Error.cv[[4]][2,1]=sqrt(sum((cv.fits.lc[[4]][,1]-cv.param.true[,1])^2)/(nCv))
Error.cv[[4]][2,2]=sqrt(sum((cv.fits.lc[[4]][,2]-cv.param.true[,2])^2)/(nCv))

Error.cv[[4]][3,1]=sqrt(sum((cv.fits.rj[[4]][,1]-cv.param.true[,1])^2)/(nCv))
Error.cv[[4]][3,2]=sqrt(sum((cv.fits.rj[[4]][,2]-cv.param.true[,2])^2)/(nCv))

sd.cv[[4]][1,1]= mean(abs(cv.fits.nn[[4]][,1]-cv.param.true[,1])/(cv.param.true[,1]))
sd.cv[[4]][1,2]= mean(abs(cv.fits.nn[[4]][,2]-cv.param.true[,2])/(cv.param.true[,2]))

sd.cv[[4]][2,1]= mean(abs(cv.fits.lc[[4]][,1]-cv.param.true[,1])/(cv.param.true[,1]))
sd.cv[[4]][2,2]= mean(abs(cv.fits.lc[[4]][,2]-cv.param.true[,2])/(cv.param.true[,2]))

sd.cv[[4]][3,1]=mean(abs(cv.fits.rj[[4]][,1]-cv.param.true[,1])/(cv.param.true[,1]))
sd.cv[[4]][3,2]=mean(abs(cv.fits.rj[[4]][,2]-cv.param.true[,2])/(cv.param.true[,2]))

names(Error.cv)=as.character(tols)
names(sd.cv)=as.character(tols)


#### Test Coverages property 

##### Neural Network ABC #####

diag1=cov.pi(param=as.data.frame(sim.params),sumstat = as.data.frame(Sumsim),
             method='neuralnet'
             ,tol=c(0.005,0.001,.01,.1),multicore = TRUE,cores=4,
             testsets = cvsamp)

##### LocLinear ABC #####

diag11=cov.pi(param=as.data.frame(sim.params),sumstat = as.data.frame(Sumsim),
              method='loclinear'
              ,tol=c(0.0001),diagnostics = c('CGR','KS'),
              testsets = cvsamp)

##### Simple Rejection ABC #####

diag3=cov.pi(param=as.data.frame(sim.params),sumstat = as.data.frame(Sumsim),
             method='rejection'
             ,tol=c(0.005,0.001,.01,.1),diagnostics = c('CGR','KS'),
             testsets = cvsamp)

