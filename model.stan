data{
int Num; //number of steps  
vector[Num] steps; // step lengths
vector[Num-1] angles; // turning angles
}

parameters {
real lambda;  
real kappa;
}

model {
  steps~exponential(lambda);
  angles~von_mises(0, kappa);
}
