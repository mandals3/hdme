library("lme4")
library("vtnorm");
#########  Y: n.sample * n.marker.
#########  ID: subject ID
calculate.ICC <- function(ID, Y){
  
  n.marker = ncol(Y); 
  groupFactor = factor(ID);
  ICC = rep(0,n.marker);
  for(k in 1:n.marker){
    dataModel <- data.frame(YY = Y[,k], group = groupFactor);
    LMEModel <- lmer(YY ~ 1 + (1 | group), data = dataModel)
    LMESummary <- VarCorr(LMEModel);
    sigmaB <- attr(LMESummary$group, "stddev")
    sigmaE <- attr(LMESummary, "sc")
    ICC[k] <- sigmaB^2 / (sigmaB^2 + sigmaE^2);
  }
  return(ICC);
  
}


#####  For data from repeated measurement experiment, 
#####  calculate correction factor = 1/(average ICC),
#####  and its standard deviation by bootstrap.
#####  Y: each row for one sample.
#####  First column is subejct ID, the other columns are data.

Estimate.correction.factor.and.variance <- function(ID, Y, n.boot){
  
  n.marker = ncol(Y);
  n.sample = nrow(Y);
  ICC = matrix(0, n.marker, 1+n.boot);
  
  ICC[,1] = calculate.ICC(ID, Y);
  unique.ID = unique(ID);
  n.subject = length(unique.ID);
  
  for(b in 1:n.boot){
    
    sampled.ID = sample(unique.ID, length(unique.ID), replace=TRUE);
    sample.index = rep(0,0);
    for(k in 1:n.subject) sample.index = c(sample.index, (1:n.sample)[ID==sampled.ID[k]]); 
    ICC[,b+1] = calculate.ICC(ID[sample.index], Y[sample.index,]);
    
    print(b)
  }
  
  ###### correction factor
  CF = mean(1/mean(ICC[,1]));
  sd = sd(1.0/apply(ICC[,2:(1+n.boot)], 2, mean));
  return(list(CF=CF, sd = sd));
}



###### Example
var.e = rbeta(100, 3, 3);
ID = c(t(cbind(1:50,1:50, 1:50)));  ### 50 subjects, each has 3 samples
Y = matrix(0, 150, 100);            ### 150 samples, 100 markers 

for(k in 1:100) Y[,k] = c(t(matrix(rnorm(50*3)*sqrt(var.e[k]), 50, 3) + rnorm(50)));  
Estimate.correction.factor.and.variance(ID, Y, 50)