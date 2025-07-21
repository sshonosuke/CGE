###-----------------------------------------------------------###
###            Code for scalable estimation of                ###
###             crossed random effects models                 ###
###                (ordered probit model)                     ###
###-----------------------------------------------------------###
## INPUT
# Y: response vector 
# X: covariate matrix (without intercept term)
# ID: matrix of clustering factors 
# m_set: numbers of factors in each clustering
# G_set: numbers of groups for multi-way grouping 
# maxitr: maximum number of iterations for the optimization algorithm 
# pen: tuning parameters for the penalty terms of the crossed random effects
# tol; tolerance rate for convergence assessment 
# print: If `T`, the iteration number is printed during the process

## OUTPUT
# beta: estimates of regression coefficients (re-estimated version)
# beta_vcov: variance-covariance matrix of the estimates of regression coefficients
# beta_sd: standard error of the estimates of regression coefficients 
# RE: estimates of random effects (re-estimated version)
# gg: estimates of grouping parameters 
# itr: number of iterations until convergence
# beta_init: initial estimates of regression coefficients (probit model without random effects)
# threshold: estimates of thresholds 

library(ordinal)

GE_CCD_ordinal <- function(Y, X, ID, m_set, G_set=NULL, maxitr=500, pen=0.01, tol=10^(-3), print=F){
  ## preparation
  if(sd(X[,1])==0){  # omit intercept (if included)
    X <- X[,-1]
  }
  N <- length(Y)
  K <- max(Y)    # class size
  J <- dim(ID)[2]
  p <- dim(X)[2]
  mat <- solve(t(X)%*%X)
  if(is.null(G_set)){
    G_set <- round(sqrt(m_set))
  }
  
  ## initial grouping
  fit_init <- clm(factor(Y)~., data=data.frame(Y=Y, X=X), link="probit")
  hBeta <- fit_init$beta
  hC <- fit_init$alpha
  hV_unq <- list()
  for(j in 1:J){
    hV_unq[[j]] <- seq(-0.5, 0.5, length=G_set[j])
  }
  
  hg <- list()
  reg <- as.vector(X%*%hBeta)
  for(j in 1:J){
    val <- matrix(NA, m_set[j], G_set[j])
    for(g in 1:G_set[j]){
      eta <- reg + hV_unq[[j]][g]   # linear predictor
      LL <- pnorm(c(hC, Inf)[Y] - eta) - pnorm(c(-Inf, hC)[Y] - eta)
      val[,g] <- as.vector(tapply(log(LL), ID[,j], sum))
    }
    hg[[j]] <- apply(val, 1, which.max)
  }
  
  hV <- list()
  for(j in 1:J){
    hV[[j]] <- hV_unq[[j]][hg[[j]]]
  }
  
  ## iteration 
  for(r in 1:maxitr){
    ## Beta
    hBeta_old <- hBeta    # current value 
    offset <- 0
    for(j in 1:J){
      offset <- offset + (hV[[j]])[ID[,j]]
    }
    eta <- as.vector(X%*%hBeta) + offset
    # gradient 
    d1 <- c(hC, Inf)[Y] - eta
    d1[d1>10^5] <- 10^5
    d0 <- c(-Inf, hC)[Y] - eta
    d0[d0<(-10^5)] <- (-10^5)
    numer <- dnorm(d1) - dnorm(d0)
    PP <- pnorm(d1) - pnorm(d0)   # P(Y_i=k)
    grad <- (-1)*apply((numer/PP)*X, 2, sum)
    # Hessian 
    multi <- (-1)*(numer^2/PP^2 + (1/PP)*(dnorm(d1)*d1-dnorm(d0)*d0)) 
    Hes <- t(X)%*%(multi*X)
    # update of Beta
    hBeta <- hBeta - solve(Hes)%*%grad
    
    ## Grouped effects 
    reg <- as.vector(X%*%hBeta) 
    for(j in 1:J){
      # offset
      offset <- reg
      for(k in 1:J){
        if(k!=j){ offset <- offset + (hV[[k]])[ID[,k]] }
      }
      # group-wise regression (intercept term only model)
      for(g in 1:G_set[j]){
        sub <- (hg[[j]][ID[,j]]==g)
        if(sum(sub)>0){
          eta <- offset + hV_unq[[j]][g]
          # penalty terms for crossed effects
          n_g <- table(factor(hg[[j]], 1:G_set[j]))
          pen_center <- (-1)*(sum(hV_unq[[j]]*n_g) - n_g*hV_unq[[j]]) / m_set[j]
          for(k in 1:J){
            if(k!=j){ pen_center <- pen_center + mean((hV[[k]])[ID[,k]]) }
          }
          pen_center <- m_set[j]*pen_center/n_g
          pen_multi <- pen*N*n_g^2/m_set[j]^2
          # gradient 
          d1 <- c(hC, Inf)[Y[sub]] - eta[sub]
          d1[d1>10^5] <- 10^5
          d0 <- c(-Inf, hC)[Y[sub]] - eta[sub]
          d0[d0<(-10^5)] <- (-10^5)
          numer <- dnorm(d1) - dnorm(d0)
          PP <- pnorm(d1) - pnorm(d0)   # P(Y_i=k)
          grad <- (-1)*sum(numer/PP) - (hV_unq[[j]][g]-pen_center[g])*pen_multi[g]
          # Hessian 
          Hes <- (-1)*sum(numer^2/PP^2 + (1/PP)*(dnorm(d1)*d1-dnorm(d0)*d0)) - pen_multi[g]
          # update 
          hV_unq[[j]][g] <- hV_unq[[j]][g] - grad/Hes
        }
      }
      hV_unq[[j]] <- hV_unq[[j]] - mean(hV_unq[[j]])
      
      # grouping 
      val <- matrix(NA, m_set[j], G_set[j])
      for(g in 1:G_set[j]){
        mu <- offset + hV_unq[[j]][g]  # linear predictor
        LL <- pnorm(c(hC, Inf)[Y] - mu) - pnorm(c(-Inf, hC)[Y] - mu)
        val[,g] <- as.vector(tapply(log(LL), ID[,j], sum))
      }
      hg[[j]] <- apply(val, 1, which.max)
      hV[[j]] <- (hV_unq[[j]])[hg[[j]]] 
    }
    
    # convergence 
    dd <- sum(abs(hBeta-hBeta_old))/sum(abs(hBeta_old)+0.001)*100
    if(dd<tol){ break() }
    if(print & r>1){ print( paste0(r, "th iteration: difference=", dd) ) }
  }
  
  # post-hoc weighting
  hV_fuzzy <- list()
  for(j in 1:J){
    val <- matrix(NA, m_set[j], G_set[j])
    offset <- as.vector(X%*%hBeta) 
    for(k in 1:J){
      if(k!=j){ offset <- offset + (hV[[k]])[ID[,k]] }
    }
    for(g in 1:G_set[j]){
      mu <- offset + hV_unq[[j]][g]
      LL <- pnorm(c(hC, Inf)[Y] - mu) - pnorm(c(-Inf, hC)[Y] - mu)
      val[,g] <- as.vector(tapply(log(LL), ID[,j], sum))
    }
    val <- val - apply(val, 1, max)
    Pi <- exp(val)/apply(exp(val), 1, sum)
    hV_fuzzy[[j]] <- apply(t(Pi)*hV_unq[[j]], 2, sum)
  }
  
  ## re-estimate of Beta
  offset <- 0
  for(j in 1:J){
    offset <- offset + (hV_fuzzy[[j]])[ID[,j]]
  }
  # iteration 
  hBeta_revised <- hBeta
  for(r2 in 1:maxitr){
    hBeta_old <- hBeta_revised    # current value 
    eta <- as.vector(X%*%hBeta_revised) + offset
    # gradient 
    d1 <- c(hC, Inf)[Y] - eta
    d1[d1>10^5] <- 10^5
    d0 <- c(-Inf, hC)[Y] - eta
    d0[d0<(-10^5)] <- (-10^5)
    numer <- dnorm(d1) - dnorm(d0)
    PP <- pnorm(d1) - pnorm(d0)   # P(Y_i=k)
    grad <- (-1)*apply((numer/PP)*X, 2, sum)
    # Hessian 
    multi <- (-1)*(numer^2/PP^2 + (1/PP)*(dnorm(d1)*d1-dnorm(d0)*d0)) 
    Hes <- t(X)%*%(multi*X)
    # update of Beta
    hBeta_revised <- hBeta_revised - solve(Hes)%*%grad
    # convergence check
    dd <- sum(abs(hBeta_revised-hBeta_old))/sum(abs(hBeta_old)+0.001)*100
    if(dd<tol){ break() }
  }
  
  # variance-covariance matrix
  beta_vcov <- (-1)*solve(Hes)
  beta_sd <- sqrt( diag(beta_vcov) )
  
  # Result 
  Result <- list(beta=hBeta_revised, beta_vcov=beta_vcov, beta_sd=beta_sd, RE=hV_fuzzy, 
                 gg=hg, itr=r, beta_init=fit_init$beta, threshold=fit_init$alpha)
  return(Result)
}


