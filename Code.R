library("lme4")
library("mvtnorm");


################################################################################
###########  Estimate the total contribution (or heritability) of all markers ##
###########  X: data matrix for all markers, samples by markers
###########  Y: phenotype vector
###########  return the estimated total contribution
################################################################################
Estimate.total.contribution <- function(X, Y){
  
  #### Y normalized to have mean zero and unit variance
  Y = (Y-mean(Y))/sd(Y); 
  temp = Y%*%t(Y);
  yiyj = temp[upper.tri(temp)];
  
  #### Each marker is normalized to have mean zero and unit variance;
  X = t(((t(X) - apply(X,2,mean))/apply(X,2,sd)));  
  temp = X %*% t(X)/ncol(X);
  gigj = temp[upper.tri(temp)];
  lambda = coef(lm(yiyj ~ gigj))[2];
  
  return(lambda);
}


################################################################################
###########  Estimate the overall correlation of beta coefficients for two traits
###########  from two independent studies with the common set of markers
###########  X1, X2: data matrix for all markers for two studies, samples by markers
###########  X1, X2 have the same set of markers with the same order
###########  Y1, Y2: phenotype vector for trait 1 and trait 2 from two studies
###########  return the estimated contribution1, contribution2, and the covariance.
###########  The correlation can be calculated as cov / sqrt(lambda1 * lambda2)
################################################################################
Estimate.covariance.two.studies <- function(X1, X2, Y1, Y2){
  
  #### Y normalized to have mean zero and unit variance
  Y1 = (Y1-mean(Y1))/sd(Y1);  
  #### Y normalized to have mean zero and unit variance
  Y2 = (Y2-mean(Y2))/sd(Y2);  
  yiyj12 = Y1 %*% t(Y2);    ### cross term.
  yiyj12 = c(t(yiyj12));
  
  #### Each marker is normalized to have mean zero and unit variance;
  X1 = t(((t(X1) - apply(X1,2,mean))/apply(X1,2,sd)));  
  #### Each marker is normalized to have mean zero and unit variance;
  X2 = t(((t(X1) - apply(X1,2,mean))/apply(X1,2,sd)));  
  temp = X1 %*% t(X2) / ncol(X1);
  gigj12 = c(t(temp));
  sigma = coef(lm(yiyj12 ~ gigj12))[2];
  
  lambda1 = Estimate.total.contribution(X1, Y1);
  lambda2 = Estimate.total.contribution(X2, Y2);
  
  return(list(lambda1 = lambda1, lambda2 = lambda2, cov = sigma))
  
}


################################################################################
#########  Calculate intraclass correlation coefficients (ICC) for repeated measurement experiment 
#########  Data: the first row is subject ID, the rest columns are the data from all markers
#########  The rows with the same ID are repeated measurements for the subject
#########  Return: ICC for all markers.
################################################################################

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


##############  For data from a repeated measurement experiment, 
##############  this function calculates correction factor = 1/(average ICC),
##############  and its standard error by bootstrap resampling of subjects.
##############  data: first column is subject ID, followed by marker data. 
##############  The rows with the same ID are repeated measurements for the subject
##############  n.boot: number of bootstrap samplings, default = 50
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



############Example: estimate total contribution for one trait using PCGC#######
n.marker = 5000;
n.causal.marker = 1000;
lambda = 0.5;  ##### total contribution
n.sample = 5000

###### Simulate data with 5000 markers, among which 1000 are causal, and 2000 subjects
X = matrix(rnorm(n.sample*n.marker), n.sample, n.marker);
beta = rnorm(n.causal.marker)*sqrt(lambda/n.causal.marker);
Y = X[,1:n.causal.marker] %*% beta + rnorm(n.sample)*sqrt(1-lambda);
######### estimate total contribution
Estimate.total.contribution(X, Y);

############Example: estimate overall correlation for two traits using PCGC#######
lambda1 = 0.6; ##### total contribution for trait 1
lambda2 = 0.4; ##### total contribution for trait 2
rho = 0.4;     ##### overall correlation
### Simulate data matrix for two studies
X1 = matrix(rnorm(n.sample*n.marker), n.sample, n.marker);
X2 = matrix(rnorm(n.sample*n.marker), n.sample, n.marker);
### Specify covariance matrix for betas for two traits
Sigma = matrix(c(lambda1, rho*sqrt(lambda1*lambda2), rho*sqrt(lambda1*lambda2), lambda2), 2,2)/n.causal.marker;
### Simulate effect sizes for all causal markers for two studies, assuming the same causal set
beta = rmvnorm(n=n.causal.marker,sigma = Sigma);
Y1 = X1[,1:n.causal.marker] %*% beta[,1] + rnorm(n.sample)*sqrt(1-lambda1);
Y2 = X2[,1:n.causal.marker] %*% beta[,2] + rnorm(n.sample)*sqrt(1-lambda2);
### Estimate covariance matrix of the effect sizes
Estimate.covariance.two.studies(X1, X2, Y1, Y2)
##################################################################################


############# Example: estimate correction factor and its standard error by bootstrap 
#############          using data from repeated measurement experiment. The example
#############          has 30 subjects with three repeats, for M=500 markers

M = 500; # number of markers
#########  Measure error standard deviation ~ beta(3,3) for M markers
#########  Average ICC = 0.8, and correction factor =1/0.8=1.25
var.e = rbeta(M, 3, 3)^2;  
ID = c(t(cbind(1:30,1:30, 1:30)));  ### 30 subjects, each has 3 samples
######## Simulate repeated measurements for 30 subjects
Z = matrix(0, 30*3, M);            
for(k in 1:M) Z[,k] = c(t(matrix(rnorm(30*3)*sqrt(var.e[k]), 30, 3) + rnorm(30)));  
######## Estimate correction factor and its standard error based on 50 bootstraps
unlist(Estimate.correction.factor.and.variance(cbind(ID, Z), 50));
##################################################################################