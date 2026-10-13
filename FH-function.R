library(MASS)
library(MCMCpack)
library(statmod)


###  function for computing benchmarked estimators under quadratic loss 
Benchmark_quadratic_loss <- function(Est, Phi){
  S <- sum(ww^2/Phi)
  return( Est + ww*(C-sum(ww*Est))/(S*Phi) )
}



###  MCMCM algorithm for Fay-Herriot model with normal random effects  ###
FH_mcmc <- function(Y, X, D, mc=2000, bn=1000, constraint=F, C=NULL){
  # preparation 
  p <- dim(X)[2]
  m <- dim(X)[1]
  mat <- solve(t(X)%*%X)
  
  # initial values 
  Beta <- coef(lm(Y~X[,-1]))
  A <- 1
  
  # MCMC 
  Beta_pos <- matrix(NA, mc, p)
  A_pos <- rep(NA, mc)
  Theta_pos <- matrix(NA, mc, m)
  
  for(k in 1:mc){
    # Theta 
    Mu <- as.vector( X%*%Beta )
    Theta <- rnorm(m, (A*Y+D*Mu)/(A+D), sqrt(A*D/(A+D)))
    Theta_pos[k,] <- Theta
    if(constraint){
      Theta_pos[k,] <- Benchmark_quadratic_loss(Est=Theta, Phi=1/D)
    }
    
    # Beta
    Beta_mean <- as.vector( mat%*%t(X)%*%Theta )
    Beta_cov <- A*mat
    Beta <- mvrnorm(1, Beta_mean, Beta_cov)
    Beta_pos[k,] <- Beta
    
    # A
    Resid <- Theta - as.vector( X%*%Beta )
    A <- rinvgamma(1, 1+m/2, 1+sum(Resid^2)/2)
    A_pos[k] <- A
  }
  
  # output
  Beta_pos <- Beta_pos[-(1:bn),]
  A_pos <- A_pos[-(1:bn)]
  Theta_pos <- Theta_pos[-(1:bn),]
  return(list(Beta=Beta_pos, A=A_pos, Theta=Theta_pos))
}





###  MCMCM algorithm for Fay-Herriot model with horseshoe random effects  ###
FH_HS_mcmc <- function(Y, X, D, mc=2000, bn=1000, constraint=F, C=NULL){
  # preparation 
  p <- dim(X)[2]
  m <- dim(X)[1]
  mat <- solve(t(X)%*%X)
  
  # initial values 
  Beta <- coef(lm(Y~X[,-1]))
  A <- 1
  U <- rep(1, m)
  
  # MCMC 
  Beta_pos <- matrix(NA, mc, p)
  A_pos <- rep(NA, mc)
  U_pos <- matrix(NA, mc, m)
  Theta_pos <- matrix(NA, mc, m)
  
  for(k in 1:mc){
    # Theta 
    Mu <- as.vector( X%*%Beta )
    Theta <- rnorm(m, (U*A*Y+D*Mu)/(U*A+D), sqrt(U*A*D/(U*A+D)))
    Theta_pos[k,] <- Theta
    if(constraint){
      Theta_pos[k,] <- Benchmark_quadratic_loss(Est=Theta, Phi=1/D)
    }
    
    # U
    Resid <- Theta - Mu
    Xi <- rinvgamma(m, 1, 1+1/U)
    U <- rinvgamma(m, 1, 1/Xi + Resid^2/(2*A))
    U[U<10^(-6)] <- 10^(-6)
    U_pos[k,] <- U
    
    # Beta
    mat <- solve( t(X)%*%(X/U) )
    Beta_mean <- as.vector( mat%*%t(X)%*%(Theta/U) )
    Beta_cov <- A*mat
    Beta <- mvrnorm(1, Beta_mean, Beta_cov)
    Beta_pos[k,] <- Beta
    
    # A
    Resid <- Theta - as.vector( X%*%Beta )
    A <- rinvgamma(1, 1+m/2, 1+sum(Resid^2/U)/2)
    A_pos[k] <- A
  }
  
  # output
  Beta_pos <- Beta_pos[-(1:bn),]
  A_pos <- A_pos[-(1:bn)]
  U_pos <- U_pos[-(1:bn),]
  Theta_pos <- Theta_pos[-(1:bn),]
  return(list(Beta=Beta_pos, A=A_pos, U=U_pos, Theta=Theta_pos))
}





###  adjustment via entropic tilting (hierarchical Bayes)
HET <- function(ww, C, th_pos, maxitr=1000){
  m <- dim(th_pos)[2]
  MC <- dim(th_pos)[1]
  gam_limit <- 10
  ESS_limit <- 500
  eq_limit <- 10^(-7)
  B_pos <- apply(ww*t(th_pos), 2, sum) - C     # constraint
  gam <- 0
  for(k in 1:maxitr){
    gam_old <- gam
    den <- exp(gam_old*B_pos)*B_pos^2
    num <- exp(gam_old*B_pos)*B_pos
    if(abs(sum(num))<eq_limit){ break }
    gam <- gam_old - mean(num)/mean(den)
    d <- abs(gam - gam_old)
    if(d<10^(-5)){ break }
    # ESS
    IW <- exp(gam*B_pos)
    IW <- IW / sum(IW)
    ESS <- 1/sum(IW^2)
    if(ESS<ESS_limit){ gam <- gam_old; break }
  }

  # importance sampling
  IW <- exp(gam*B_pos)
  IW <- IW / sum(IW)
  ESS <- 1/sum(IW^2)
  sel <- sort(sample(1:MC, prob=IW, replace=T))
  th_cpos <- th_pos[sel,]
  CB <- apply(th_cpos, 2, mean)
  return(list(pos=th_cpos, CB=CB, itr=k, gam=gam, ESS=ESS))
}




###  adjustment via entropic tilting (empirical Bayes)
EET <- function(Y, X, ww, C, fit, maxitr=1000){
  m <- length(ww)
  hbeta <- apply(fit$Beta, 2, mean)
  hA <- mean(fit$A)
  if(is.null(fit$U)){ hU <- rep(1, m) 
  }else{ hU <- apply(fit$U, 2, mean) }
  PM <- apply(fit$Theta, 2, mean)
  PV <- apply(fit$Theta, 2, var)
  mPM <- PM + ww*PV/sum(ww*PV)*(C-sum(ww*PM))
  # posterior sample
  MC <- length(fit$A)
  pos <- matrix(NA, MC, m)
  for(i in 1:m){
    pos[,i] <- rnorm(MC, mPM[i], sqrt(PV[i]))
  }
  # output
  return(list(PM=mPM, PV=PV, pos=pos))
}




###  MDI: adjustment via entropic tilting and Gaussian approximation
MDI <- function(Y, X, D, C, ww, Theta_pos){
  MC <- dim(Theta_pos)[1]
  Sigma <- var(Theta_pos)
  Tht <- apply(Theta_pos, 2, mean)
  Gam <- (C-sum(ww*Tht))/as.vector(t(ww)%*%Sigma%*%ww)
  cPM <- Tht + Gam*as.vector(Sigma%*%ww)
  Theta_cpos <- mvrnorm(MC, cPM, Sigma)
  return(list(pos=Theta_cpos, CB=cPM))
}


