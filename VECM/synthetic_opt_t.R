rm(list = ls())
set.seed(4321)

library(tsDyn)
library(mnormt)
library(urca)
library(vars)
library(tseries)
library(expm)
library(Rdonlp2)
library(timeSeries)
library(QRM)
library(mvtnorm)

setwd("~/Desktop/final project/codes")
source("models/VECM/synthetic_t_parametertest.R")

GetPortExp <- function(weight, assets.mean) {
  port.mean = matrix(weight, nrow=1) %*% matrix(assets.mean, ncol=1)
  return(port.mean)
}


GetPortVar <- function(weight, var.matrix) {
  port.var = matrix(weight, nrow=1) %*% var.matrix %*% matrix(weight, ncol=1)
  return(port.var)
}


Optimize <- function(asset.mean, asset.var, gamma) {
  ObjectiveFunction <- function(weight) {
    return(-GetPortExp(weight, asset.mean) + gamma * GetPortVar(weight, asset.var))
  }
  
  num.assets <- length(asset.mean)
  
  # lower bound and upper bound
  par.l = rep(0, num.assets)
  par.u = rep(0.3, num.assets)
  
  # constraints sum of weight = 1
  A = matrix(rep(1, num.assets), nrow=1)
  lin.l = 1
  lin.u = 1
  
  # initial weight
  init.weight <- rep(1/num.assets, num.assets)
  
  # do optimization
  ret = donlp2(init.weight, ObjectiveFunction, par.u=par.u, par.l=par.l, A, lin.l=lin.l, lin.u=lin.u)
  opt.weight = ret$par
  
  return(opt.weight)
}

GetFinalReturn <- function(init.price, final.price, weight) {
  final.return <- matrix(final.price-init.price, nrow=1) %*% matrix(weight, ncol=1)
  return(c(final.return))
}

#======== Data Frame ============== (all.data + new.data)
data.group1 <- GenerateSampleData(pi[[1]], df, num.total.days)
data.group2 <- GenerateSampleData(pi[[2]],df, num.total.days)
data.group3 <- GenerateSampleData(pi[[3]],df, num.total.days)

# estimated parameters (1000 for each)
fit.result.group1<-FitModel(data.group1, num.train.days)
fit.result.group2<-FitModel(data.group2, num.train.days)
fit.result.group3<-FitModel(data.group3, num.train.days)
pi.hat.group1 <- fit.result.group1[[1]]
pi.hat.group2 <- fit.result.group2[[1]]
pi.hat.group3 <- fit.result.group3[[1]]
df.hat.group1 <- fit.result.group1[[2]]
df.hat.group2 <- fit.result.group2[[2]]
df.hat.group3 <- fit.result.group3[[2]]
sigma.hat.group1 <- fit.result.group1[[3]]
sigma.hat.group2 <- fit.result.group2[[3]]
sigma.hat.group3 <- fit.result.group3[[3]]

#========= Group Exp and Var ============
EstimateMean <- function(pi,num.test.days, init.price) {
  identity <- diag(3)
  mean.est <- (identity + pi) %^% (num.test.days) %*% init.price - init.price
  return(mean.est)
}

EstimateVar <- function(pi, R, df, num.test.days) {
  identity = diag(3)
  var.est <- matrix(0, nrow=3, ncol=3)
  for (k in 0:(num.test.days - 1)) {
    var.est = var.est + t((identity + pi) %^% (k)) %*% R %*% ((identity + pi) %^% (k))
  }
  var.est = var.est * df / (df-2)
  return(var.est)
}

#========= Group expectation =======
est.mean <- list()
for(i in 1: 1000){
  est.mean.group1 <- EstimateMean(pi.hat.group1[[i]] ,5,init.price)
  est.mean.group2 <- EstimateMean(pi.hat.group2[[i]] ,5,init.price)
  est.meam.group3 <- EstimateMean(pi.hat.group3[[i]] ,5,init.price)
  est.mean[[i]] <- c(est.mean.group1,est.mean.group2,est.meam.group3)
}

#========== Group Variance ==========
est.var <- list()
for(i in 1:1000){
  est.var.group1 <- as.matrix(EstimateVar(pi.hat.group1[[i]],sigma.hat.group1[[i]], df.hat.group1[[i]], 5))
  est.var.group2 <- as.matrix(EstimateVar(pi.hat.group2[[i]],sigma.hat.group2[[i]],df.hat.group2[[i]], 5))
  est.var.group3 <- as.matrix(EstimateVar(pi.hat.group3[[i]],sigma.hat.group3[[i]], df.hat.group3[[i]],5))
  est.var[[i]] <- matrix(0, ncol = 9, nrow = 9)
  est.var[[i]][1:3,1:3] <- est.var.group1
  est.var[[i]][4:6,4:6] <- est.var.group2
  est.var[[i]][7:9,7:9] <- est.var.group3
}

#=========portfolio opt weight=============
GetWeight<- function(asset.mean, asset.var,gamma){
  opt.weight <- list()
  for(i in 1:1000){
    opt.weight[[i]]<- Optimize(asset.mean[[i]],asset.var[[i]],gamma)
  }
  return(opt.weight)
}

#========initial and final prices for final rewards===========
initial.price <- matrix(0, ncol = 9, nrow = 1000)
final.price <- matrix(0, ncol = 9, nrow = 1000)
for(i in 1:1000){
  initial.price.group1 <- data.group1[[1]][[i]][100,]
  initial.price.group2 <- data.group2[[1]][[i]][100,]
  initial.price.group3 <- data.group3[[1]][[i]][100,]
  initial.price[i,] <- c(initial.price.group1,initial.price.group2,initial.price.group3)
  final.price.group1 <- data.group1[[1]][[i]][105,]
  final.price.group2 <- data.group2[[1]][[i]][105,]
  final.price.group3 <- data.group3[[1]][[i]][105,]
  final.price[i,] <- c(final.price.group1,final.price.group2,final.price.group3)
}

#========Get final Return =============
finalStep <- function(gamma){
  opt.weight <- GetWeight(est.mean, est.var, gamma)
  final.return <- matrix(0, ncol = 1000, nrow = 8)
  for(i in 1: 1000){
    return <- GetFinalReturn(initial.price[i,], final.price[i,],opt.weight[[i]])
    final.return[,i] <- return
  }
  return(final.return)
}

gamma <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.1, 1.0, 10)

final.density<- list()
for(i in 1:8){
  return<-finalStep(gamma[i])
  final.density[[i]] = density(return)
}

# Plot density graphs
plot(final.density[[1]], type = 'l', col = '1', main = 'Estimated parameters', xlab = 'Expected Reward',ylim = c(0,0.4))
for (i in 2:8) {
  lines(final.density[[i]], col =i)
  legend("bottomright",legend = c("0.0001", "0.0005", "0.001", "0.005", "0.01", "0.1", "1.0", "10"), 
         lty=c(1,1,1,1,1,1,1,1),col=c('1','2','3','4','5','6','7','8'), bty = 'n', cex = 0.7)
}

#####################################
# True parameters
#####################################
#========= Group expectation =======
true.mean.group1 <- EstimateMean(pi[[1]],5,init.price)
true.mean.group2 <- EstimateMean(pi[[2]],5,init.price)
true.meam.group3 <- EstimateMean(pi[[3]],5,init.price)
true.mean <- c(true.mean.group1,true.mean.group2,true.meam.group3)
#========== Group Variance ==========
true.var.group1 <- EstimateVar(pi[[1]],R, df, 5)
true.var.group2 <- EstimateVar(pi[[2]],R,df, 5)
true.var.group3 <- EstimateVar(pi[[3]],R,df, 5)
true.var<- matrix(0, ncol = 9, nrow = 9)
true.var[1:3,1:3] <- true.var.group1
true.var[4:6,4:6] <- true.var.group2
true.var[7:9,7:9] <- true.var.group3

#=========portfolio opt weight=============
GetWeight.true<- function(asset.mean, asset.var,gamma){
  opt.weight<- Optimize(asset.mean,asset.var,gamma)
  return(opt.weight)
}

#========Get final Return =============
finalStep.true <- function(gamma){
  opt.weight <- GetWeight.true(true.mean, true.var, gamma)
  final.return <- matrix(0, ncol = 1000, nrow = 8)
  for(i in 1: 1000){
    return <- GetFinalReturn(initial.price[i,], final.price[i,],opt.weight)
    final.return[,i] <- return
  }
  return(final.return)
}

final.density.true<- list()
for(i in 1:8){
  return.true<-finalStep.true(gamma[i])
  final.density.true[[i]] = density(return.true)
}

# Plot density graphs
plot(final.density.true[[1]], type = 'l', col = '1', main = 'True parameters', xlab = 'Expected Reward',ylim = c(0,0.5))
for (i in 2:8) {
  lines(final.density.true[[i]], col =i)
  legend("bottomright",legend = c("0.0001", "0.0005", "0.001", "0.005", "0.01", "0.1", "1.0", "10"), 
         lty=c(1,1,1,1,1,1,1,1),col=c('1','2','3','4','5','6','7','8'), bty = 'n', cex = 0.7)
}


