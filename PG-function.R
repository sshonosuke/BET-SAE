library(MASS)
library(MCMCpack)


###  MCMCM algorithm for Poisson-gamma model   ###
PG_mcmc <- function(Y, X, ni, mc=2000, bn=1000, constraint=F, C=NULL, ww=NULL){
  # preparation 
  p <- dim(X)[2]
  m <- dim(X)[1]
  
  # prior 
  B0 <- rep(0, p)
  S0 <- rep(10^4, p)
  nu0 <- eta0 <- 1
  
  # initial values 
  Z <- Y*ni
  init_fit <- glm(Z~X[,-1], offset=log(ni), family="poisson")
  Beta <- coef(init_fit)
  band <- summary(init_fit)$coefficient[,2]
  nu <- 5
  
  # MCMC 
  Beta_pos <- matrix(NA, mc, p)
  nu_pos <- rep(NA, mc)
  Lam_pos <- matrix(NA, mc, m)
  Lam_cpos1 <- Lam_cpos2 <- matrix(NA, mc, m)
  
  for(k in 1:mc){
    # Lam
    Mu <- as.vector( exp(X%*%Beta) )
    Lam <- rgamma(m, Z+nu*Mu, ni+nu)
    Lam_pos[k,] <- Lam
    if(constraint){
      Lam_cpos1[k,] <- Lam + C - sum(ww*Lam) 
      Lam_cpos2[k,] <- Lam*C/sum(ww*Lam) 
    }
    
    # Beta
    Beta_new <- Beta + band*rnorm(p)
    Mu_new <- as.vector( exp(X%*%Beta_new) )
    logdens <- sum(dgamma(Lam, nu*Mu, nu, log=T)) + sum(dnorm(Beta, B0, sqrt(S0), log=T))
    logdens_new <- sum(dgamma(Lam, nu*Mu_new, nu, log=T)) + sum(dnorm(Beta_new, B0, sqrt(S0), log=T))
    pp <- min(1, exp(logdens_new-logdens))
    if(runif(1)<pp){ 
      Beta <- Beta_new
    }
    Beta_pos[k,] <- Beta
    Mu <- as.vector( exp(X%*%Beta) ) 
    
    # nu
    nu_new <- nu + runif(1, -0.2, 0.2)
    if(nu_new<0.001){ nu_new <- 0.001 }
    logdens <- sum(dgamma(Lam, nu*Mu, nu, log=T)) + dgamma(nu, nu0, eta0, log=T)
    logdens_new <- sum(dgamma(Lam, nu_new*Mu, nu_new, log=T)) + dgamma(nu_new, nu0, eta0, log=T)
    pp <- min(1, exp(logdens_new-logdens))
    if(runif(1)<pp){ 
      nu <- nu_new
    }
    nu_pos[k] <- nu
  }
  
  # output
  Beta_pos <- Beta_pos[-(1:bn),]
  nu_pos <- nu_pos[-(1:bn)]
  Lam_pos <- Lam_pos[-(1:bn),]
  Lam_cpos1 <- Lam_cpos1[-(1:bn),]
  Lam_cpos2 <- Lam_cpos2[-(1:bn),]
  return(list(Beta=Beta_pos, nu=nu_pos, Lam=Lam_pos, cLam1=Lam_cpos1, cLam2=Lam_cpos2))
}




### adjustment via entropic tilting (hierarchical Bayes)
HET <- function(ww, C, th_pos, maxitr=1000){
  m <- dim(th_pos)[2]
  MC <- dim(th_pos)[1]
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
EET <- function(Y, X, ww, C, fit){
  m <- length(ww)
  hbeta <- apply(fit$Beta, 2, mean)
  hmu <- c( exp(X%*%hbeta) )
  hnu <- mean(fit$nu)
  
  fn <- function(gam){
    Score <- (ni*Y+hnu*hmu) / (ni+hnu+gam*ww)
    return(sum(ww*Score)-C)
  }
  
  gam <- 0
  if(fn(-100)*fn(100)<0){
    gam <- uniroot(f=fn, interval=c(-100, 100))$root
  }
  if(fn(-500)*fn(500)<0){
    gam <- uniroot(f=fn, interval=c(-500, 500))$root
  }
  if(fn(-1000)*fn(1000)<0){
    gam <- uniroot(f=fn, interval=c(-1000, 1000))$root
  }
  
  # posterior sample
  A <- ni*Y+hnu*hmu
  B <- ni+hnu+gam*ww
  MC <- length(fit$nu)
  pos <- matrix(NA, MC, m)
  for(i in 1:m){
    pos[,i] <- rgamma(MC, A[i], B[i])
  }
  # output
  return(list(PM=A/B, pos=pos, gam=gam))
}







