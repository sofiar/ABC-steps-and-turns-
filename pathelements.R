# Function to estimate the turning angle and step length 
# Returns: "steps": steps lengths of the trajectory
#          "turns": turning angles of the trajectory
#          "direction": directions of the trajectory
#          "cosine": cosine of the angle
#          "sine": sine of the angle

pathelements <- function(X,Y){
  n = length(X)
  adj = X[2:n]-X[1:n-1] 
  op = Y[2:n]-Y[1:n-1]
  step = (adj^2 + op^2)^0.5
  aca = which(step==0)
  step[aca]=0.0000001
  
  si<-sign(adj)
  si[si==0]<-1  #corrects for sign(0) == 0
  ang = si*(op<0)*pi+atan(adj/op)
  adif <- ang[2:length(ang)]-ang[1:(length(ang)-1)]
  
  ## we correct so that it is between -pi and pi
  adif[which(adif < -pi)] = adif[which(adif < -pi)] + 2*pi
  adif[which(adif > pi)] = adif[which(adif > pi)] - 2*pi  
  
  turns <- adif
  direction <- ang  
  co = adj/step # cosine
  si = op/step  # sine
  res = list(step,turns,direction,co,si)
  names(res) = c("steps","turns","direction","cosine","sine")
  return(res)
}


  


