###-----------------------------------------------------------------###
###          Code for application of ordered probit model           ###
###      with two-way random effects to 100K MovieLens data         ###
###-----------------------------------------------------------------###

rm(list=ls())
source("GE-CCD-ordinal.R")

## Dataset 
load("Data.RData")   

Y <- rating_matrix[,3]
X_mat <- (rating_matrix[,-(1:4)])
age_original <- as.numeric(rating_matrix[,4])
age <- scale( age_original )
knot <- quantile(age, prob=(1:9)/10)
age_spline <- cbind(age, age^2)
for(l in 1:9){
  age_spline <- cbind(age_spline, (age-knot[l])^2*(age>=knot[l]))
}
X <- cbind(age_spline, X_mat[,-c(19,32)])
p <- dim(X)[2]

uID <- rating_matrix[,1]
mID <- rating_matrix[,2]
ID <- cbind(uID, mID)
m_set <- apply(ID, 2, max)


## Grouped effects model 
CPT <- c()

tt <- proc.time()[1]
fit_GE1 <- GE_CCD_ordinal(Y=Y, X=X, ID=ID, m_set=m_set, G_set=c(20, 20), maxitr=500, print=T)
CPT[1] <- proc.time()[1] - tt

tt <- proc.time()[1]
fit_GE2 <- GE_CCD_ordinal(Y=Y, X=X, ID=ID, m_set=m_set, tol=10^(-3), maxitr=500, print=T)
CPT[2] <- proc.time()[1] - tt

## standard probit regression 
tt <- proc.time()[1]
fit_init <- clm(factor(Y)~., data=data.frame(Y=Y, X=X), link="probit")
CPT[3] <- proc.time()[1] - tt


## computation time 
names(CPT) <- c("CGE1", "CGE2", "Probit")


## Summary 
V <- fit_GE1$beta_vcov
xx <- seq(min(age), max(age), length=1000)
xx_spline <- cbind(xx, xx^2)
for(l in 1:9){
  xx_spline <- cbind(xx_spline, (xx-knot[l])^2*(xx>=knot[l]))
}
qq <- dim(age_spline)[2]
SD_spline <- sqrt(diag(xx_spline%*%V[1:qq,1:qq]%*%t(xx_spline)))

fn0 <- xx_spline%*%fit_init$beta[1:qq]
fn1 <- xx_spline%*%fit_GE1$beta[1:qq]
fn2 <- xx_spline%*%fit_GE2$beta[1:qq]

RE1 <- fit_GE1$RE[[1]]
RE2 <- fit_GE1$RE[[2]]


## packages for ggplot
library(ggplot2)
library(dplyr) 
library(gridExtra)

## Figure (coefficients)
pdf("App-reg.pdf", height=8, width=16, pointsize=16)
xx_plot <- sd(age_original)*xx + mean(age_original)

df1 <- data.frame(
  Age=rep(xx_plot, 4),
  Effect=c(fn0, fn1, fn1-1.96*SD_spline, fn1+1.96*SD_spline),
  PlotModel=rep(c("WoE", "CGE", "CGE-l", "CGE-u"), each=length(xx_plot)),
  LegendModel=rep(c("WoE", "CGE", "CGE-CI", "CGE-CI"), each=length(xx_plot))
)

p1 <- ggplot(df1, aes(x=Age, y=Effect, color=LegendModel, linetype=LegendModel, group=PlotModel)) +
  geom_line(size=0.8) +
  scale_color_manual(values=c("WoE"="red", "CGE"="black", "CGE-CI"="black")) +
  scale_linetype_manual(values=c("WoE"="solid", "CGE"="solid", "CGE-CI"="dotted")) +
  theme_minimal(base_size=22) +
  labs(x="Age", y="Effect", color=NULL, linetype=NULL) +
  theme(legend.position="bottom")

df2 <- data.frame(Regression=fit_init$beta[-(1:qq)], CGE=fit_GE1$beta[-(1:qq)],
  ymin=fit_GE1$beta[-(1:qq)]-1.96*fit_GE1$beta_sd[-(1:qq)],
  ymax=fit_GE1$beta[-(1:qq)]+1.96*fit_GE1$beta_sd[-(1:qq)],
  Model="Regression vs CGE1")
p2 <- ggplot(df2, aes(x=Regression, y=CGE)) +
  geom_errorbar(aes(ymin=ymin, ymax=ymax), width=0.05, color="blue") +
  geom_point(aes(color=Model, shape=Model), size=4) +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  scale_color_manual(values="blue") +
  scale_shape_manual(values=18) +
  theme_minimal(base_size=22) +
  labs(x="Estimate (WoE)", y="Estimate (CGE)") +
  theme(legend.position="none")

grid.arrange(p1, p2, nrow=1)
dev.off()



## Figure (random effects)
pdf("App-RE.pdf", height=6, width=12)
df <- data.frame(value=c(RE1, RE2), effect=rep(c("User effect", "Movie effect"), c(length(RE1), length(RE2))))
ggplot(df, aes(x=value)) +
  geom_histogram(bins=30, fill="steelblue", color="white") +
  facet_wrap(~effect, scales="free", nrow=1) +
  theme_minimal(base_size=17) +
  labs(x="", y="Frequency")
dev.off()
