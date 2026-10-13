###   Simulation using Fay-Herriot model (scenario B)    ###
rm(list=ls())
library(progress)

source("FH-function.R")
set.seed(1)

## settings 
scenario <- 3    # scenario of true random effects (1-3)
R <- 2000        # number of replications
m <- 30          # number of areas
G <- 5
A <- 0.3
X <- cbind(1, runif(m, -2, 2), rbinom(m, 1, 0.5))
ni <- rep(c(30, 20, 15, 10, 5), rep(m/G, G))
ww <- ni/sum(ni)   # weight 
N <- sum(ni)
D <- rep(c(0.3, 0.7, 1, 1.5, 2), rep(m/G, G))
p <- dim(X)[2]


## Monte Carlo 
meth1 <- c("HB", "HET", "EET", "MDI", "PP", "WL1", "WL2")
M1 <- length(meth1)
meth2 <- c("HB", "HET", "EET", "MDI", "PP")
M2 <- length(meth2)

MSE <- array(NA, c(R, M1, 2))
CP <- AL <- array(NA, c(R, M2, 2))
Gam <- matrix(NA, R, 2)
Data_Y <- matrix(NA, R, m)

dimnames(MSE)[[2]] <- meth1
dimnames(CP)[[2]] <- dimnames(AL)[[2]] <- meth2
dimnames(MSE)[[3]] <- dimnames(CP)[[3]] <- dimnames(AL)[[3]] <- c("Normal", "Horseshoe")

pb <- progress_bar$new(total=R)
for(r in 1:R){
  # data generation 
  mu <- (-2) + (X[,3]+1)*X[,2]^2
  
  if(scenario==1){
    Theta <- mu + sqrt(A)*rnorm(m)
    Y <- rnorm(m, Theta, sqrt(D))
  }
  
  if(scenario==2){
    Theta <- mu + rbinom(m, 1, 0.3)*sqrt(A)*rnorm(m) 
    Y <- rnorm(m, Theta, sqrt(D))
  }
  
  if(scenario==3){
    nu <- 2.5
    Theta <- mu + sqrt(A)*rt(m, nu)/sqrt(nu/(nu-2))
    Y <- rnorm(m, Theta, sqrt(D))
  }
  
  Data_Y[r,] <- Y
  
  # benchmarking target 
  C <- mean(Theta)
  
  # MCMC
  fit1 <- FH_mcmc(Y, X, D, mc=10000, bn=2000)
  fit2 <- FH_HS_mcmc(Y, X, D, mc=10000, bn=2000)
  PM1 <- apply(fit1$Theta, 2, mean)
  PM2 <- apply(fit2$Theta, 2, mean)
  
  # HB entropic tilting
  cpos1 <- HET(ww, C, th_pos=fit1$Theta)
  cpos2 <- HET(ww, C, th_pos=fit2$Theta)
  CB1 <- cpos1$CB
  CB2 <- cpos2$CB
  
  Gam[r,] <- c(cpos1$gam, cpos2$gam)
  
  
  # EB entropic tilting
  epos1 <- EET(Y, X, ww, C, fit=fit1)
  epos2 <- EET(Y, X, ww, C, fit=fit2)
  ECB1 <- epos1$PM
  ECB2 <- epos2$PM
    
  # entropic tilting with Gaussian approximation
  mpos1 <- MDI(Y, X, D, C, ww, Theta_pos=fit1$Theta)
  mpos2 <- MDI(Y, X, D, C, ww, Theta_pos=fit2$Theta)
  mCB1 <- mpos1$CB
  mCB2 <- mpos2$CB
  
  # posterior projection
  ppos1 <- FH_mcmc(Y, X, D, constraint=T, C=C)
  ppos2 <- FH_HS_mcmc(Y, X, D, constraint=T, C=C)
  pPM1 <- apply(ppos1$Theta, 2, mean)
  pPM2 <- apply(ppos2$Theta, 2, mean)
 
  # Quadratic loss 
  QL1_pw <- Benchmark_quadratic_loss(Est=PM1, Phi=1/apply(fit1$Theta, 2, var))
  QL2_pw <- Benchmark_quadratic_loss(Est=PM2, Phi=1/apply(fit2$Theta, 2, var))
  
  QL1_dw <- Benchmark_quadratic_loss(Est=PM1, Phi=1/D)
  QL2_dw <- Benchmark_quadratic_loss(Est=PM2, Phi=1/D)
  
  # MSE
  Est1 <- cbind(PM1, CB1, ECB1, mCB1, pPM1, QL1_pw, QL1_dw)
  Est2 <- cbind(PM2, CB2, ECB2, mCB2, pPM2, QL2_pw, QL2_dw)
  MSE[r,,1] <- apply((Est1-Theta)^2, 2, mean)
  MSE[r,,2] <- apply((Est2-Theta)^2, 2, mean)
  
  print(MSE[r,,1])
  print(MSE[r,,2])
  
  # coverage 
  Pos1 <- list(fit1$Theta, cpos1$pos, mpos1$pos, epos1$pos, ppos1$Theta)
  Pos2 <- list(fit2$Theta, cpos2$pos, mpos2$pos, epos2$pos, ppos2$Theta)
  for(model in 1:M2){
    CI <- apply(Pos1[[model]], 2, quantile, prob=c(0.025, 0.975))
    CP[r,model,1] <- 100*mean(CI[1,]<Theta & CI[2,]>Theta)
    AL[r,model,1] <- mean(CI[2,]-CI[1,])
    CI <- apply(Pos2[[model]], 2, quantile, prob=c(0.025, 0.975))
    CP[r,model,2] <- 100*mean(CI[1,]<Theta & CI[2,]>Theta)
    AL[r,model,2] <- mean(CI[2,]-CI[1,])
  }

  # print
  pb$tick()
}


# save
save(list=ls(), file=paste0("FH-Result-B-s", scenario, ".RData"))

