###   Simulation using Poisson-gamma model (scenario B)    ###
rm(list=ls())
library(progress)

source("PG-function.R")
set.seed(1)

## settings 
scenario <- 3    # scenario of true random effects (1-3)
R <- 2000        # number of replications
m <- 30          # number of areas
G <- 5
nu <- 20
X <- cbind(1, runif(m, -2, 2), rbinom(m, 1, 0.5))
ni <- rep(c(5, 10, 15, 20, 30), rep(m/G, G))
ww <- ni/sum(ni)   # weight 
p <- dim(X)[2]


## Monte Carlo 
meth <- c("HB", "HET", "EET", "PP1", "PP2")
M <- length(meth)
MSE <- CP <- AL <- matrix(NA, R, M)
dimnames(MSE)[[2]] <- dimnames(CP)[[2]] <- dimnames(AL)[[2]] <- meth

pb <- progress_bar$new(total=R)
for(r in 1:R){
  # data generation 
  mu <- exp(abs(X[,2])*X[,3])
  
  if(scenario==1){
    Lam <- rgamma(m, nu*mu, nu)
  }
  
  if(scenario==2){
    ch <- rbinom(m, 1, 0.8)
    Lam <- ch*rgamma(m, nu*mu, nu) + (1-ch)*mu
  }
  
  if(scenario==3){
    Lam <- mu*exp(0.5*rnorm(m))
  }

  C <- mean(Lam)   # benchmark target 
  Y <- rpois(m, ni*Lam)/ni
  
  # MCMC
  fit <- PG_mcmc(Y, X, ni, mc=10000, bn=2000, constraint=T, C=C, ww=ww)
  PM <- apply(fit$Lam, 2, mean)
  
  # entropic tilting (HB)
  fit_ET <- HET(ww=ww, C=C, th_pos=fit$Lam)
  cpos <- fit_ET$pos
  CB <- fit_ET$CB
  
  # entropic tilting (EB)
  epos <- EET(Y, X, ww, C, fit=fit)
  ECB <- epos$PM
  cpos_EB <- epos$pos
  
  # posterior projection (with Datta's benchmarking method)
  D1_pos <- fit$cLam1
  D2_pos <- fit$cLam2
  Datta1 <- apply(D1_pos, 2, mean)
  Datta2 <- apply(D2_pos, 2, mean)
  
  # MSE
  Est <- cbind(PM, CB, ECB, Datta1, Datta2)
  MSE[r,] <- apply((Est-Lam)^2, 2, mean)

  # CP & AL
  Pos <- list(fit$Lam, cpos, cpos_EB, D1_pos, D2_pos)
  for(model in 1:M){
    if(is.na(mean(Pos[[model]]))==F){
      CI <- apply(Pos[[model]], 2, quantile, prob=c(0.025, 0.975))
      CP[r,model] <- 100*mean(CI[1,]<Lam & CI[2,]>Lam)
      AL[r,model] <- mean(CI[2,]-CI[1,])
    }
  }

  # print
  pb$tick()
}


# save
save(list=ls(), file=paste0("PG-Result-B-s", scenario, ".RData"))

