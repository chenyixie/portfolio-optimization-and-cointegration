library(tsDyn)
library(mnormt)
library(urca)
library(vars)
library(tseries)
library(expm)
library(mvtnorm)
library(QRM)

# clear environment and set random seed
rm(list = ls())
set.seed(4321)

# parameters
num.asset = 3     
num.group = 3    
num.train.days = 100
num.test.days  = 5  
num.total.days =num.train.days + num.test.days

#============================================================================
# 1. Generate synthetic data
#============================================================================
# Generate data for one group 
GenerateGroupData <- function(pi, R, df, num.days, init.price) {
  num.asset <- length(init.price)
  all.data <- matrix(0, ncol = num.asset, nrow = num.days)
  data.new    <- matrix(0, ncol = num.asset, nrow = num.days)
  all.data[1,] <- init.price
  for (i in 2:num.days){
    all.data[i, ] <- pi %*% all.data[(i-1),] + all.data[(i-1), ] + matrix(rmvt(n=1, sigma=R, df=df),ncol = 1)
  }
  data.new <- all.data
  data.new[1,] <- all.data[1,]
  for(i in 2:dim(all.data)[1]){
    data.new[i,] <- all.data[i,] - all.data[(i-1),] - pi %*% all.data[(i-1),]
  }
  return(list(all.data,data.new))
}

GenerateSampleData <- function(pi,df,num.days){
  all.data <- list()
  new.data <- list()
  for( i in 1:1000){
    all.data[[i]] <- GenerateGroupData(pi, R,df, num.days, init.price)[[1]]
    new.data[[i]] <- GenerateGroupData(pi, R,df, num.days, init.price)[[2]]
  }
  return(list(all.data,new.data))
}

#============================================================================
# 2. fit VECM model
#============================================================================

FitModel<- function(all.data,num.days){
  train.data<- list()
  train.data.new <- list()
  pi.hat <- list()
  df.hat <- list()
  R.hat <- list()
  vecm.fit <- list()
  s <- list()
  alpha.hat <- list()
  beta.hat <- list()
  pi.hat <- list()
  t.fit <- list()
  for(i in 1:1000){
    train.data[[i]] <- all.data[[1]][[i]][1:num.days,]
    train.data.new[[i]] <- all.data[[2]][[i]][1: num.days,]
    vecm.fit[[i]] <- VECM(train.data[[i]], r=1, lag=0, estim="ML", include = "none")
    s[[i]] <- summary(vecm.fit[[i]])
    alpha.hat[[i]] <- s[[i]]$coefficients
    beta.hat[[i]] <- s[[i]]$model.specific$beta
    pi.hat[[i]] <- alpha.hat[[i]] %*% t(beta.hat[[i]])
    t.fit[[i]] <- fit.mst(train.data.new[[i]], method="BFGS")
    df.hat[[i]] <- t.fit[[i]]$df
    R.hat[[i]] <- t.fit[[i]]$Sigma
  }
  return(list(pi.hat, df.hat, R.hat))
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
df = 3

GetMSE <- function(vector1, vector2) {
  vector1 <- as.vector(vector1)
  vector2 <- as.vector(vector2)
  
  return(mean((vector1 - vector2) ^ 2))
}

#=========Plot for histograms ==============
data.group1<- GenerateSampleData(pi[[1]],df, num.train.days)
fit.result <- FitModel(data.group1, num.train.days)
pi.hat <- fit.result[[1]]
df.hat <- fit.result[[2]]
sigma.hat <- fit.result[[3]]
pi.MSE <- c()
sigma.MSE <- c()
df.MSE <- c()
for(i in 1:1000){
  pi.MSE[i] <- GetMSE(pi.hat[[i]], pi[[1]])
  sigma.MSE[i] <- GetMSE(sigma.hat[[i]], R)
  df.MSE[i] <- GetMSE(df.hat[[i]], df)
}
hist(pi.MSE, xlab = 'MSE of pi',main = 'T = 100')
hist(sigma.MSE,  xlab = 'MSE of Sigma',main = 'T = 100')
hist(df.MSE,  xlab = 'MSE of df',main = 'T = 100')







