########  For a repeat measurement experiment (n.sample subjects, n.repeat repeats),
########  we calculate  


library("lme4")
library("vtnorm");
Estimate.total.contribution <- function(X, Y){
  
  Y = (Y-mean(Y))/sd(Y);  #### Y normalized to have mean zero and unit variance
  temp = Y%*%t(Y);
  yiyj = temp[upper.tri(temp)];
  
  X = t(((t(X) - apply(X,2,mean))/apply(X,2,sd)));  #### Each marker is normalized to have mean zero and unit variance;
  temp = X %*% t(X)/ncol(X);
  gigj = temp[upper.tri(temp)];
  lambda = coef(lm(yiyj ~ gigj))[2];
  
  return(lambda);
}

Estimate.covariance.two.studies <- function(X1, X2, Y1, Y2){
  
  Y1 = (Y1-mean(Y1))/sd(Y1);  #### Y normalized to have mean zero and unit variance
  Y2 = (Y2-mean(Y2))/sd(Y2);  #### Y normalized to have mean zero and unit variance
  yiyj12 = Y_1 %*% t(Y_2);    ### cross term.
  yiyj12 = c(t(yiyj12));
  
  X1 = t(((t(X1) - apply(X1,2,mean))/apply(X1,2,sd)));  #### Each marker is normalized to have mean zero and unit variance;
  X2 = t(((t(X1) - apply(X1,2,mean))/apply(X1,2,sd)));  #### Each marker is normalized to have mean zero and unit variance;
  temp = X1 %*% t(X2) / ncol(X1);
  gigj12 = c(t(temp));
  sigma = coef(lm(yiyj12 ~ gigj12))[2];
  
  lambda1 = Estimate.total.contribution(X1, Y1);
  lambda2 = Estimate.total.contribution(X2, Y2);
  
  return(list(lambda1 = lambda1, lambda2 = lambda2, cov = sigma))
  
}


n.marker = 5000;
n.causal.marker = 1000;
lambda = 0.5;
lambda1 = 0.6;
lambda2 = 0.4;
rho = 0.4;

##### lambda
X = matrix(rnorm(n.sample*n.marker), n.sample, n.marker);
beta = rnorm(n.causal.marker)*sqrt(lambda/n.causal.marker);
Y = X[,1:n.causal.marker] %*% beta + rnorm(n.sample)*sqrt(1-lambda);
Estimate.total.contribution(X, Y)

##### rho, two studies
X1 = matrix(rnorm(n.sample*n.marker), n.sample, n.marker);
X2 = matrix(rnorm(n.sample*n.marker), n.sample, n.marker);
Sigma = matrix(c(lambda1, rho*sqrt(lambda1*lambda2), rho*sqrt(lambda1*lambda2), lambda2), 2,2)/n.causal.marker;
beta = rmvnorm(n=n.causal.marker,sigma = Sigma);
Y1 = X1[,1:n.causal.marker] %*% beta[,1] + rnorm(n.sample)*sqrt(1-lambda1);
Y2 = X2[,1:n.causal.marker] %*% beta[,2] + rnorm(n.sample)*sqrt(1-lambda2);
Estimate.covariance.two.studies(X1, X2, Y1, Y2)
