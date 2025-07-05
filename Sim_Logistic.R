###-----------------------------------------------------------###
###              Code for one-shot simulation                 ###
###              under two-way binary outcome                 ###
###-----------------------------------------------------------###
rm(list=ls())
set.seed(123)

library(lme4)
library(glmmTMB)
source("clubbed_backfitting_schall.R")   # available at https://github.com/G28Sw/scalable-logistic-regression-with-crossed-random-effects
source("CGE.R")   # proposed method 

# settings 
scenario <- 1    # 1 or 2
N <- 5000       # sample size

J <- 2    # number of effects
p <- 5    # number of covariates
n1 <- n2 <- 2*round(sqrt(N))
m_set <- c(n1, n2) 
ID <- cbind(sample(1:n1, N, replace=T), sample(1:n2, N, replace=T))

# true coefficients
Beta <- c(1, -1, 0.5, rep(0, p-2))


# random effects
if(scenario==1){
  V1 <- 0.5*rnorm(n1)
  V2 <- rnorm(n2)
}
if(scenario==2){
  V1 <- rgamma(n1, 1, 1) - 1
  V2 <- 1 - rgamma(n2, 1, 1)
}
V <- list(V1, V2)

# covariate 
X <- cbind(1, matrix(rnorm(N*p), N, p))

# response 
Mu <- as.vector(X%*%Beta) + V1[ID[,1]] + V2[ID[,2]] 
Prob <- 1/(1+exp(-Mu))
Y <- rbinom(N, 1, Prob)


# Grouped effects model
CPT <- c()
tt <- proc.time()[1]
fit_GE <- GE_CCD(Y=Y, X=X, ID=ID, m_set=m_set, tol=10^(-3), maxitr=200, print=F)
CPT[1] <- proc.time()[1] - tt
fit_GE$itr

# Mixed model (Laplace approximation)
tt <- proc.time()[1]
fit_GLMM1 <- glmer(Y~X[,-1]+(1|ID[,1])+(1|ID[,2]), family="binomial", nAGQ=1)
CPT[2] <- proc.time()[1] - tt

# Mixed model (PQL)
tt <- proc.time()[1]
fit_GLMM2 <- glmer(Y~X[,-1]+(1|ID[,1])+(1|ID[,2]), family="binomial", nAGQ=0)
CPT[3] <- proc.time()[1] - tt

hV_GLMM1 <- hV_GLMM2 <- list()
for(j in 1:J){
  hV_GLMM1[[j]] <- coef(fit_GLMM1)[[j]][,1] - summary(fit_GLMM1)$coefficients[1,1]
  hV_GLMM2[[j]] <- coef(fit_GLMM2)[[j]][,1] - summary(fit_GLMM2)$coefficients[1,1]
}

# glmmTMB
data <- data.frame(Y=Y, X=X[,-1], ID1=ID[,1], ID2=ID[,2])
names(data)[2:(p+1)] <- paste0("X", 2:(p+1))
tt <- proc.time()[1]
fit_GLMM_tmb <- glmmTMB(Y~X2+X3+X4+X5+X6+(1|ID1)+(1|ID2), family=binomial, data=data)
CPT[4] <- proc.time()[1] - tt

hV_TMB <- list()
for(j in 1:J){
  hV_TMB[[j]] <- as.numeric(ranef(fit_GLMM_tmb)$cond[[j]][,1])
}

# Backfitting method by Ghosh et al. 
tt <- proc.time()[1]
fit_BF <- backfitting_outer_loop(X, Y, ID[,1], ID[,2])
CPT[5] <- proc.time()[1] - tt


# computation time
names(CPT) <- c("CGE", "MLE1", "MLE0", "TMB", "BF")
CPT


# Beta
Beta_est <- matrix(NA, p, 5)
dimnames(Beta_est)[[2]] <- c("CGE", "MLE1", "MLE0", "TMB", "BF")
Beta_est[,1] <- fit_GE$beta
Beta_est[,2] <- summary(fit_GLMM1)$coefficients[-1,1]
Beta_est[,3] <- summary(fit_GLMM2)$coefficients[-1,1]
Beta_est[,4] <- unname(fixef(fit_GLMM_tmb)$cond)[-1]
Beta_est[,5] <- fit_BF$beta[-1]

# squared error loss
100*apply((Beta_est - Beta[-1])^2, 2, mean)


# 95% Confidence interval
CI <- list()
CI[[1]] <- cbind(fit_GE$beta-1.96*fit_GE$beta_sd, fit_GE$beta+1.96*fit_GE$beta_sd)
sd_GLMM1 <- summary(fit_GLMM1)$coefficients[-1,2]
CI[[2]] <- cbind(Beta_est[,2]-1.96*sd_GLMM1, Beta_est[,2]+1.96*sd_GLMM1)
sd_GLMM2 <- summary(fit_GLMM2)$coefficients[-1,2]
CI[[3]] <- cbind(Beta_est[,3]-1.96*sd_GLMM2, Beta_est[,3]+1.96*sd_GLMM2)
sd_TMB <- summary(fit_GLMM_tmb)$coefficients$cond[,"Std. Error"][-1]
CI[[4]] <- cbind(Beta_est[,4]-1.96*sd_TMB, Beta_est[,4]+1.96*sd_TMB)

for(l in 1:4){
  dimnames(CI[[l]])[[1]] <- paste0("X", 1:p)
  dimnames(CI[[l]])[[2]] <- c("lower", "upper")
}
names(CI) <- c("CGE", "MLE1", "MLE0", "TMB")
CI


