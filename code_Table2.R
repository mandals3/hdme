
rm(list = ls())

####################################################################################################

#function for estimating lambda:

Estimate.total.contribution <- function(X, Y){
	Y <- scale(Y) #### Y normalized to have mean zero and unit variance
	temp <- Y %*% t(Y)
	yiyj <- temp[upper.tri(temp)]
	X <- scale(X) #### Each marker is normalized to have mean zero and unit variance
	temp <- X %*% t(X) / ncol(X)
	gigj <- temp[upper.tri(temp)]
	lambda <- coef(lm(yiyj ~ gigj))[2]
	return(lambda)
}

####################################################################################################

#parameter settings:

n.sample <- 5000
n.marker <- 10000

nCausalMarkerVec <- n.marker * c(1, 0.2, 0.1, 0.05)
nNCausalMarker <- length(nCausalMarkerVec)

lambdaVec <- c(0.1, 0.5, 0.9)
nLambda <- length(lambdaVec)

####################################################################################################

#simulations:
#(It will be time consuming if n.sample and n.marker are large. It can be modified by using parallel computing to save time.)

nSimu <- 500
resultMat1 <- resultMat2 <- resultMat3 <- matrix(NA, nSimu, nNCausalMarker * nLambda)

for(i in 1:nSimu){
	set.seed(i)
	for(j in 1:nNCausalMarker){
		n.causal.marker <- nCausalMarkerVec[j]
		for(k in 1:nLambda){
			lambda <- lambdaVec[k]
			X <- matrix(rnorm(n.sample * n.marker), n.sample, n.marker)
			beta <- rnorm(n.causal.marker) * sqrt(lambda / n.causal.marker)
			Y <- X[, 1:n.causal.marker] %*% beta + rnorm(n.sample) * sqrt(1 - lambda)
			X1 <- X + rnorm(n.sample * n.marker, 0, 0.5) #theta^2 = 0.25
			X2 <- X + rnorm(n.sample * n.marker, 0, 1) #theta^2 = 1
			s <- (j - 1) * nLambda + k
			resultMat1[i, s] <- Estimate.total.contribution(X1, Y)
			resultMat2[i, s] <- Estimate.total.contribution(X2, Y)
			resultMat3[i, s] <- Estimate.total.contribution(X, Y)
			cat(i, j, k, s, resultMat1[i, s], resultMat2[i, s], resultMat3[i, s], "\n")
		}
	}
}
#str(resultMat1)
#str(resultMat2)
#str(resultMat3)

####################################################################################################

#summary and output:

tempMeanVec <- colMeans(cbind(resultMat1, resultMat2, resultMat3))
tempSdVec <- apply(cbind(resultMat1, resultMat2, resultMat3), 2, sd)
tempWrite <- cbind(paste0(round(tempMeanVec[1:12], 3), " (", round(tempSdVec[1:12], 3), ")"), paste0(round(tempMeanVec[13:24], 3), " (", round(tempSdVec[13:24], 3), ")"))
#write.table(tempWrite, file = "result_Table2.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

####################################################################################################


