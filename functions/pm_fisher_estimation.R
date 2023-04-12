## Numerical study in the Poisson mixture model
## R function for Fisher Information matrix estimation 

pm_fisher_estimation <- function(y, n, lambda, alpha) {
    
    # y      : vector of observations
    # n      : sample size
    # lambda : vector of mean Poisson parameters for each component of the mixture
    # alpha  : vector of mixture proportions, excluding the proportion of the last mixture component
    
    K <- length(lambda)
    
    deriv1ind  <- matrix(0,2*K-1,n) 
    deriv2     <- matrix(0,2*K-1,2*K-1) 
    covderiv   <- matrix(0,2*K-1,2*K-1)
    
    ## computation of conditional expectation of the first derivatives of the complete data log-likelihood
    
    denom <- 0
    for (k in 1:(K-1)){
        denom <- denom + exp(-lambda[k])*lambda[k]^y*alpha[k]
    }
    denom <- denom + exp(-lambda[K])*lambda[K]^y*(1-sum(alpha))
    
    for (k in 1:(K-1)){
        deriv1ind[k,]   <- exp(-lambda[k])*lambda[k]^y*alpha[k]/denom * (y/lambda[k]-1)
        deriv1ind[K+k,] <- exp(-lambda[k])*lambda[k]^y*alpha[k]/denom/alpha[k] + 
            exp(-lambda[K])*lambda[K]^y*(1-sum(alpha))/denom*(-1/(1-sum(alpha)))
    }
    deriv1ind[K,] <- exp(-lambda[K])*lambda[K]^y*(1-sum(alpha))/denom * (y/lambda[K]-1)
    
    
    ## computation of conditional expectation of the second derivatives of the complete data log-likelihood
    
    
    for (k in 1:(K-1)){
        deriv2[k,k]     <- sum(exp(-lambda[k])*lambda[k]^y*alpha[k]/denom * (-y/lambda[k]^2))
        deriv2[K+k,K+k] <- sum(-exp(-lambda[k])*lambda[k]^y*alpha[k]/denom/alpha[k]^2 
                               + exp(-lambda[K])*lambda[K]^y*(1-sum(alpha))/denom*(-1/(1-sum(alpha))^2))
    }
    
    for (k in 1:(K-2)){
        for (l in (k+1):(K-1)){
            deriv2[K+k,K+l] <- sum(exp(-lambda[K])*lambda[K]^y*(1-sum(alpha))/denom*(-1/(1-sum(alpha))^2))
            deriv2[K+l,K+k] <- deriv2[K+k,K+l] 
        }
    }
    
    deriv2[K,K]<-sum(exp(-lambda[K])*lambda[K]^y*(1-sum(alpha))/denom * (-y/lambda[K]^2))
    
    ## computation of the conditional covariance matrix of the first derivatives of the complete data log-likelihood
    
    
    for (k in 1:(K-2)){
        covderiv[k,k] <- sum(exp(-lambda[k])*lambda[k]^y*alpha[k]/denom*(-1+y/lambda[k])^2) 
        covderiv[k+K,k+K] <- sum(exp(-lambda[k])*lambda[k]^y*alpha[k]/denom/alpha[k]^2
                                 +exp(-lambda[K])*lambda[K]^y*(1-sum(alpha))/denom/(1-sum(alpha))^2) 
        for (l in (k+1):(K-1)){
            covderiv[k+K,l+K] <- sum(exp(-lambda[K])*lambda[K]^y*(1-sum(alpha))/denom/(1-sum(alpha))^2) 
            covderiv[l+K,k+K] <- covderiv[k+K,l+K]
        } 
        covderiv[k,K+k] <- sum(exp(-lambda[k])*lambda[k]^y*alpha[k]/denom/alpha[k]*(-1+y/lambda[k])) 
        covderiv[k+K,k] <- covderiv[k,K+k]
        
        covderiv[K,K+k] <- sum(exp(-lambda[K])*lambda[K]^y*(1-sum(alpha))/denom*(-1+y/lambda[K])*(-1)/(1-sum(alpha))) 
        covderiv[K+k,K] <- covderiv[K,K+k]
    }
    
    covderiv[K-1,K-1] <- sum(exp(-lambda[K-1])*lambda[K-1]^y*alpha[K-1]/denom*(-1+y/lambda[K-1])^2) 
    covderiv[2*K-1,2*K-1] <- sum(exp(-lambda[K-1])*lambda[K-1]^y*alpha[K-1]/denom/alpha[K-1]^2+exp(-lambda[K])*lambda[K]^y*(1-sum(alpha))/denom/(1-sum(alpha))^2) 
    covderiv[K-1,2*K-1] <- sum(exp(-lambda[K-1])*lambda[K-1]^y*alpha[K-1]/denom/alpha[K-1]*(-1+y/lambda[K-1])) 
    covderiv[2*K-1,K-1] <- covderiv[K-1,2*K-1]
    
    covderiv[K,2*K-1] <- sum(exp(-lambda[K])*lambda[K]^y*(1-sum(alpha))/denom*(-1+y/lambda[K])*(-1)/(1-sum(alpha))) 
    covderiv[2*K-1,K] <- covderiv[K,2*K-1]
    
    covderiv[K,K] <- sum(exp(-lambda[K])*lambda[K]^y*(1-sum(alpha))/denom*(-1+y/lambda[K])^2) 
    
    
    Isco <- deriv1ind%*%t(deriv1ind)/n
    Iobs <- deriv1ind%*%t(deriv1ind)/n - deriv2/n - covderiv/n
    
    
    res <- list(Isco = Isco, Iobs = Iobs)
    
    return(res)
}

