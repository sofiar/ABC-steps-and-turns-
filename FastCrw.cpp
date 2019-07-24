// This Rcpp function generates the correlated random walk path with parameters k and w
// Returns: x: x-coordinates 
//          y: y-coordinates 
//          t: times of that coordinates
#include <Rcpp.h>
using namespace Rcpp;

DataFrame cppFastCRW_exp(double k, double w, double ns, double maxx){

  RNGScope scp;
  Function rvonmises("rvonmises");
  Function rexp("rexp");
  Function runif("runif");
  Function cumsum("cumsum");
  
  double nk = k;
  double nw = w;
  double nsteps = ns;
  double mm = maxx;
  
  NumericVector vtm(nsteps);
  NumericVector vls(nsteps);
  
  NumericVector aux2(1);
  NumericVector aux(1);
  
  NumericVector tu(nsteps);
  NumericVector x(nsteps);
  NumericVector y(nsteps);
  NumericVector t(nsteps);
  NumericVector d(nsteps);
  NumericVector phi(nsteps);
  NumericVector temp(nsteps);
  
  vtm=rexp(nsteps,nw);
  phi=rvonmises(nsteps, 0, nk);
  d=cumsum(phi); 
  double temp2=0;
  int i=1;
  while((i<nsteps) && (temp2<mm)){
  t[i] = t[i-1] + vtm[i];
  x[i] = x[i-1] + cos(d[i]) * vtm[i];
  y[i] = y[i-1] + sin(d[i]) * vtm[i]; 
  temp[i]=999;
  temp2=t[i];
  i++;   
    
  }
  
  return DataFrame::create(Named("x")= x[temp>0],Named("y")= y[temp>0],
                           Named("t")= t[temp>0]);
  
}
  
  
  