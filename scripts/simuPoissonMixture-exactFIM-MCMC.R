## Numerical study in the Poisson mixture model
## Monte-Carlo estimation of the true Fisher information matrix based on a large sample
## Note to the reviewers : This part of the numerical experiment requires a lot
## of computing time and is therefore not executed directly. Reducing nMC allows
## us to get back to more reasonable computation times but degrades the quality
## of the estimation of the FIM by Monte-Carlo. The result of the Monte-Carlo
## approximation of the FIM obtained with nMC=10^8 is stored in the file
## PoissonMixtureTrueFIM.Rdata which is loaded in the following chunks.

nMC <- 100000000

z <- rep(0,nMC) # latent variable
y <- rep(0,nMC) # observed variable

t <- cumsum(c(alpha,1-sum(alpha)))

for (i in 1:nMC){
    u    <- runif(1)
    z[i] <- 1+sum(u>t)
    y[i] <- rpois(1,lambda[z[i]])
}

trueFIM <-  pm_fisher_estimation(y, nMC, lambda, alpha)
trueFIM$Iobs
trueFIM$Isco
trueFIM <- (trueFIM$Isco+trueFIM$Iobs)/2
save(trueFIM,file='Rfiles/PoissonMixtureTrueFIM.Rdata')
