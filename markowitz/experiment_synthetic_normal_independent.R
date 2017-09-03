# Experiments for synthetic data
# 2017.8.10

rm(list=ls())
set.seed(4321)
setwd("~/Desktop/final project/codes")
source("shared_functions.R")
source("models/markowitz/normal_independent.R")

mu <- rnorm(9, 0, 0.01)
sigma <- abs(rnorm(9, 0, 0.2)) * diag(9)

init.price <- abs(rnorm(9, 5, 5))
num.days <- 500

all.mu.mse <- rep(0, 1000)
all.sigma.mse <- rep(0, 1000)

for (i in 1:1000) {
  all.data <- GenerateData(mu, sigma, init.price, num.days)
  fit.result <- FitModel(all.data)
  mu.hat <- fit.result$mu
  sigma.hat <- fit.result$sigma
  
  mu.mse <- GetMSE(mu, mu.hat)
  sigma.mse <- GetMSE(diag(sigma), diag(sigma.hat))
  
  all.mu.mse[i] <- mu.mse
  all.sigma.mse[i] <- sigma.mse
}

hist(all.mu.mse, main = 'T = 500', xlab = 'MSE of mu')
hist(all.sigma.mse, main = 'T = 500', xlab = 'MSE of sigma')

