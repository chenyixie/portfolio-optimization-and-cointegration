library(tsDyn)
library(urca)
library(vars)
library(tseries)
library(expm)
library(Rdonlp2)
library(PerformanceAnalytics)
library(QRM)

# clear environment and set random seed
rm(list = ls())
set.seed(4321)

# parameters
num.asset = 3     
num.group = 3    
num.train.days = 100
num.test.days  = 1
num.total.days =num.train.days + num.test.days

#============================================================================
# 1. Generate synthetic data
#============================================================================
# Generate data for one group 
GenerateGroupData <- function(pi, R, num.days, init.price) {
  num.asset <- length(init.price)
  all.data <- matrix(0, ncol = num.asset, nrow = num.days)
  all.data[1,] <- init.price
  for (i in 2:num.days){
    all.data[i, ] <- pi %*% matrix(all.data[(i-1),], ncol = 1) + matrix(all.data[(i-1), ],ncol = 1) + matrix(rmnorm(n = 1, mu = 0, Sigma = R),ncol =1)
  }
  return(all.data)
}

GenerateSampleData <- function(pi, num.days){
  all.data <- list()
  for( i in 1:1000){
    all.data[[i]] <- GenerateGroupData(pi, R, num.days, init.price)
  }
  return(all.data)
}

#generate data for 3 groups
GenerateData <- function(num.group){
  all.data <- rep(0, num.total.days*num.asset*num.group)
  dim(all.data) <- c(num.total.days, num.asset, num.group)
  for(i in 1: num.group){
    all.data[,,i] <- GenerateGroupData(pi[[i]],R,num.total.days, init.price)
  }
  return(all.data)
}

#============================================================================
# 2. fit VECM model
#============================================================================

FitModel<- function(all.data, num.days){
  train.data<- list()
  vecm.fit <- list()
  s <- list()
  R.hat <- list()
  alpha.hat <- list()
  beta.hat <- list()
  pi.hat <- list()
  for(i in 1:1000){
    train.data[[i]] <- all.data[[i]][1:num.days,]
    vecm.fit[[i]] <- VECM(train.data[[i]], r=1, lag=0, estim="ML", include="none")
    s[[i]] <- summary(vecm.fit[[i]])
    R.hat[[i]] <- cov(s[[i]]$residuals)
    alpha.hat[[i]] <- s[[i]]$coefficients
    beta.hat[[i]] <- s[[i]]$model.specific$beta
    pi.hat[[i]] <- alpha.hat[[i]] %*% t(beta.hat[[i]])
  }
  return(list(pi.hat, R.hat))
}

#============================================================================
# 3. Experiemnt 
#============================================================================
pi<- list()
pi[[1]] <- matrix(c(-0.5, -0.6, 0.6), nrow = num.asset) %*% matrix(c(0.5, 0.7, -0.9), ncol = num.asset)
pi[[2]] <- matrix(c(-0.8, -0.7, 0.9), nrow = num.asset) %*% matrix(c(0.7, 0.4, -0.92), ncol = num.asset)
pi[[3]] <- matrix(c(-0.45,-0.7, 0.5), nrow = num.asset) %*% matrix(c(0.4, 0.6, -0.7), ncol = num.asset)

r <- matrix(rnorm(num.asset* num.asset, 0, 0.5), ncol=num.asset, nrow=num.asset)
R <- r %*% t(r)

init.price <- abs(rnorm(3, 10, 2))

GetMSE <- function(vector1, vector2) {
  vector1 <- as.vector(vector1)
  vector2 <- as.vector(vector2)
  
  return(mean((vector1 - vector2) ^ 2))
}

#==========Plot Pricing Trends==============
group1.single <- GenerateGroupData(pi[[1]], R, num.train.days, init.price)
plot(group1.single[,1], type = 'l', main = "Cointegrating Assets", xlab = 'time', ylab = 'Price', ylim = c(0,16))
lines(group1.single[,2], col = 2)
lines(group1.single[,3], col = 4)
legend("bottomright",legend = c("Asset 1", "Asset 2", "Asset 3"), 
       lty=c(1,1,1),col=c('1','2','4'), bty = 'n', cex = 0.7)

#=========Plot for histograms ==============
data.group1<- GenerateSampleData(pi[[1]], num.train.days)
fit.result <- FitModel(data.group1, num.train.days)
pi.hat <- fit.result[[1]]
sigma.hat <- fit.result[[2]]
pi.MSE <- c()
sigma.MSE <- c()
for(i in 1:1000){
  pi.MSE[i] <- GetMSE(pi.hat[[i]], pi[[1]])
  sigma.MSE[i] <- GetMSE(sigma.hat[[i]], R)
}
hist(pi.MSE, xlab = 'MSE of pi',main = 'T = 100')
hist(sigma.MSE,  xlab = 'MSE of Sigma',main = 'T = 100')

