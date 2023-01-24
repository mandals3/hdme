library("lme4")
library("mvtnorm");

calculate.ICC <- function(data){
  
  n.marker = ncol(data) - 1; 
  groupFactor = factor(data[,1]);
  Y = data[,2:ncol(data)];
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


##############  For data from repeated measurement experiment, 
##############  calculate correction factor = 1/(average ICC),
##############  and its standard error by bootstrap resampling
##############  of subjects.
##############  data: first column is subject ID, followed by marker data 
##############  n.bbot: number of bootstrap, default = 100
Estimate.correction.factor.and.variance <- function(data, n.boot = 50){
  
  ID = data[,1];
  unique.ID = unique(ID);
  n.subject = length(unique.ID);
  
  n.marker = ncol(data) - 1;
  n.sample = length(ID);
  ICC = matrix(0, n.marker, 1+n.boot);
  
  ICC[,1] = calculate.ICC(data);
  
  
  for(b in 1:n.boot){
    
    sampled.ID = sample(unique.ID, length(unique.ID), replace=TRUE);
    sample.index = rep(0,0);
    my.ID = rep(0,0);
    for(k in 1:n.subject){
      sample.index = c(sample.index, (1:n.sample)[ID==sampled.ID[k]]); 
      my.ID = c(my.ID, rep(k, sum(ID==sampled.ID[k])));
    }
    ICC[,b+1] = calculate.ICC(cbind(my.ID,data[sample.index,2:ncol(data)]));
    
    print(b)
  }
  
  ###### correction factor
  CF = mean(1/mean(ICC[,1]));
  sd = sd(1.0/apply(ICC[,2:(1+n.boot)], 2, mean));
  return(list(CF=CF, sd = sd));
}



#######  M: # of markers; 
####### 30 subjects, each 2 repeated samples, residue variance = beta(3,3)
Simulation.Table.1a <- function(M, n.boot = 50){
  
  var.e = rbeta(M, 3, 3)^2;  ##### average measurement error variance = 0.5, CF=1.5
  ID = c(t(cbind(1:30,1:30)));  ### 30 subjects, each has 2 repeated samples
  Y = matrix(0, 30*2, M);            
  for(k in 1:M) Y[,k] = c(t(matrix(rnorm(30*2)*sqrt(var.e[k]), 30, 2) + rnorm(30)));  
  result = unlist(Estimate.correction.factor.and.variance(cbind(ID, Y), n.boot));
  
  return(result);
  
}

####### 30 subjects, each 3 repeated samples, residue variance = beta(3,3)
Simulation.Table.1b <- function(M, n.boot = 50){
  
  var.e = rbeta(M, 3, 3)^2;  ##### average measurement error variance = 0.5, CF=1.5
  ID = c(t(cbind(1:30, 1:30, 1:30)));  ### 30 subjects, each has 3 samples
  Y = matrix(0, 30*3, M);            
  for(k in 1:M) Y[,k] = c(t(matrix(rnorm(30*3)*sqrt(var.e[k]), 30, 3) + rnorm(30)));  
  result = unlist(Estimate.correction.factor.and.variance(cbind(ID, Y), n.boot));
  
  return(result);
  
}

####### 30 subjects, each 2 repeated samples, residue variance = beta(5,2)
Simulation.Table.1c <- function(M, n.boot = 50){
  
  var.e = rbeta(M, 5, 2)^2;  ##### average measurement error variance = 0.71, CF=1.71
  ID = c(t(cbind(1:30, 1:30)));  ### 30 subjects, each has 2 samples
  Y = matrix(0, 30*2, M);            
  for(k in 1:M) Y[,k] = c(t(matrix(rnorm(30*2)*sqrt(var.e[k]), 30, 2) + rnorm(30)));  
  result = unlist(Estimate.correction.factor.and.variance(cbind(ID, Y), n.boot));
  
  return(result);
  
}

####### 30 subjects, each 3 repeated samples, residue variance = beta(5,2)
Simulation.Table.1d <- function(M, n.boot = 50){
  
  
  var.e = rbeta(M, 5, 2)^2;  ##### average measurement error variance = 0.71, CF=1.71
  ID = c(t(cbind(1:30,1:30, 1:30)));  ### 30 subjects, each has 3 samples
  Y = matrix(0, 30*3, M);            
  for(k in 1:M) Y[,k] = c(t(matrix(rnorm(30*3)*sqrt(var.e[k]), 30, 3) + rnorm(30)));  
  result = unlist(Estimate.correction.factor.and.variance(cbind(ID, Y), n.boot));
  
  return(result);
  
}

#######  M: # of markers; 
####### 30 subjects, each 2 repeated samples, residue variance = 1
Simulation.Table.1e <- function(M, n.boot = 50){
  
  var.e = rep(1, M);  ##### average measurement error variance = 1, CF=2
  ID = c(t(cbind(1:30,1:30)));  ### 30 subjects, each has 2 repeated samples
  Y = matrix(0, 30*2, M);            
  for(k in 1:M) Y[,k] = c(t(matrix(rnorm(30*2)*sqrt(var.e[k]), 30, 2) + rnorm(30)));  
  result = unlist(Estimate.correction.factor.and.variance(cbind(ID, Y), n.boot));
  
  return(result);
  
}

#######  M: # of markers; 
####### 30 subjects, each 2 repeated samples, residue variance = 0.25
Simulation.Table.1f <- function(M, n.boot = 50){
  
  var.e = rep(0.25, M);  ##### average measurement error variance = 0.25 CF=1.25
  ID = c(t(cbind(1:30,1:30)));  ### 30 subjects, each has 2 repeated samples
  Y = matrix(0, 30*2, M);            
  for(k in 1:M) Y[,k] = c(t(matrix(rnorm(30*2)*sqrt(var.e[k]), 30, 2) + rnorm(30)));  
  result = unlist(Estimate.correction.factor.and.variance(cbind(ID, Y), n.boot));
  
  return(result);
  
}
