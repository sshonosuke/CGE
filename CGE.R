###-----------------------------------------------------------###
###            Code for scalable estimation of                ###
###             crossed random effects models                 ###
###             (logistic and Poisson models)                 ###
###-----------------------------------------------------------###
## INPUT
# Y: response vector 
# X: covariate matrix (without intercept term)
# ID: matrix of clustering factors 
# m_set: numbers of factors in each clustering
# G_set: numbers of groups for multi-way grouping 
# family: distribution of response ("binomial" or "poisson")
# maxitr: maximum number of iterations for the optimization algorithm 
# pen: tuning parameters for the penalty terms of the crossed random effects
# tol; tolerance rate for convergence assessment 
# print: If `T`, the iteration number is printed during the process

## OUTPUT
# beta: estimates of regression coefficients (re-estimated version)
# beta_sd: standard error of the estimates of regression coefficients 
# Int: intercept estimate 
# RE: estimates of random effects (re-estimated version)
# gg: estimates of grouping parameters 
# itr: number of iterations until convergence

GE_CCD <- function(Y, X, ID, m_set, G_set=NULL, family="binomial", 
                   maxitr=100, pen=0.01, tol=10^(-3), print=F){
  ## preparation
  if(sd(X[,1])==0){  # omit intercept (if included)
    X <- X[,-1]
  }
  J <- dim(ID)[2]
  p <- dim(X)[2]
  mat <- solve(t(X)%*%X)
  if(is.null(G_set)){
    G_set <- round(sqrt(m_set))
  }
  N <- length(Y)    # total sample size
  
  ## initial grouping
  fit_init <- coef(glm(Y~X, family=family))
  hBeta <- fit_init[-1]
  Int <- fit_init[1]
  hV_unq <- list()
  for(j in 1:J){
    hV_unq[[j]] <- seq(-0.5, 0.5, length=G_set[j])
  }

  hg <- list()
  for(j in 1:J){
    val <- matrix(NA, m_set[j], G_set[j])
    for(g in 1:G_set[j]){
      mu <- Int + as.vector(X%*%hBeta) + hV_unq[[j]][g]
      if(family=="binomial"){ 
        prob <- 1/(1+exp(-mu))
        LL <- dbinom(Y, 1, prob, log=T) 
      }
      if(family=="poisson"){ 
        lam <- exp(mu)
        LL <- dpois(Y, lam, log=T) 
      }
      val[,g] <- as.vector(tapply(LL, ID[,j], sum))
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
    hBeta_old <- hBeta
    offset <- Int
    for(j in 1:J){
      offset <- offset + (hV[[j]])[ID[,j]]
    }
    hBeta <- coef(glm(Y~X-1, family=family, offset=offset))
    reg <- as.vector(X%*%hBeta) 
    
    ## Grouped effects 
    for(j in 1:J){
      # offset and probability
      offset <- Int + reg
      for(k in 1:J){
        if(k!=j){ offset <- offset + (hV[[k]])[ID[,k]] }
      }
      Eta <- Int + reg
      for(k in 1:J){ Eta <- Eta + (hV[[k]])[ID[,k]] }
      
      # penalty terms for crossed effects
      n_g <- table(factor(hg[[j]], 1:G_set[j]))
      pen_center <- (-1)*(sum(hV_unq[[j]]*n_g) - n_g*hV_unq[[j]]) / m_set[j]
      for(k in 1:J){
        if(k!=j){ pen_center <- pen_center + mean((hV[[k]])[ID[,k]]) }
      }
      pen_center <- m_set[j]*pen_center/n_g
      pen_multi <- pen*N*n_g^2/m_set[j]^2
      
      # iteratively weighted least squares 
      if(family=="binomial"){
        prob <- 1/(1+exp(-Eta))
        ww <- prob*(1-prob)
        zz <- Eta + (Y-prob)/ww
      }
      if(family=="poisson"){
        ww <- exp(Eta)
        zz <- Eta + (Y-ww)/ww
      }
      numer <- ww*(zz-offset)
      for(g in 1:G_set[j]){
        sub <- (hg[[j]][ID[,j]]==g)
        if(sum(sub)>0){
          hV_unq[[j]][g] <- (sum(numer[sub]) + pen_multi[g]*pen_center[g]) / (sum(ww[sub]) + pen_multi[g])
        }
      }
      
      # grouping 
      val <- matrix(NA, m_set[j], G_set[j])
      for(g in 1:G_set[j]){
        mu <- offset + hV_unq[[j]][g]
        if(family=="binomial"){ 
          prob <- 1/(1+exp(-mu))
          LL <- dbinom(Y, 1, prob, log=T)
        }
        if(family=="poisson"){ 
          lam <- exp(mu)
          LL <- dpois(Y, lam, log=T) 
        }
        val[,g] <- as.vector(tapply(LL, ID[,j], sum))
      }
      hg[[j]] <- apply(val, 1, which.max)
      hV[[j]] <- (hV_unq[[j]])[hg[[j]]] 
    }
    
    # adjustment 
    for(j in 1:J){
      n_g <- table(factor(hg[[j]], 1:G_set[j]))
      mm <- sum(n_g*hV_unq[[j]]) / m_set[j]
      Int <- Int + mm
      hV_unq[[j]] <- hV_unq[[j]] - mm/m_set[j]
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
    offset <- Int + as.vector(X%*%hBeta) 
    for(k in 1:J){
      if(k!=j){ offset <- offset + (hV[[k]])[ID[,k]] }
    }
    for(g in 1:G_set[j]){
      mu <- offset + hV_unq[[j]][g]
      if(family=="binomial"){ 
        prob <- 1/(1+exp(-mu))
        LL <- dbinom(Y, 1, prob, log=T)
      }
      if(family=="poisson"){ 
        lam <- exp(mu)
        LL <- dpois(Y, lam, log=T) 
      }
      val[,g] <- as.vector(tapply(LL, ID[,j], sum))
    }
    val <- val - apply(val, 1, max)
    Pi <- exp(val)/apply(exp(val), 1, sum)
    hV_fuzzy[[j]] <- apply(t(Pi)*hV_unq[[j]], 2, sum)
  }
  
  # re-estimate of Beta
  offset <- Int
  for(j in 1:J){
    offset <- offset + (hV_fuzzy[[j]])[ID[,j]]
  }
  fit_revised <- glm(Y~X-1, family=family, offset=offset)
  hBeta_revised <- coef(fit_revised)
  hBeta_sd <- summary(fit_revised)$coefficients[,2]
  
  # Result 
  Result <- list(beta=hBeta_revised, beta_sd=hBeta_sd, Int=Int, 
                 RE=hV_fuzzy, gg=hg, itr=r)
  return(Result)
}


