# Experiments for synthetic data
# 2017.8.11

rm(list=ls())
library(clusterGeneration)
set.seed(4321)
setwd("~/Desktop/final project/codes")
source("shared_functions.R")
source("models/markowitz/t_dependent.R")

mu <- rnorm(9, 0, 0.01)
sigma <- genPositiveDefMat(9)$Sigma
degree <- 3

init.price <- abs(rnorm(9, 5, 5))
num.days <- 1000
num.samples <- 1000

all.mu.mse <- rep(0, num.samples)
all.sigma.mse <- rep(0, num.samples)
all.degree.mse <- rep(0, num.samples)

for (i in 1:num.samples) {
  all.data <- GenerateData(degree, mu, sigma, init.price, num.days)
  fit.result <- FitModel(all.data)
  mu.hat <- fit.result$mu
  sigma.hat <- fit.result$sigma
  degree.hat <- fit.result$degree
  
  mu.mse <- GetMSE(mu, mu.hat)
  sigma.mse <- GetMSE(sigma, sigma.hat)
  degree.mse <- GetMSE(degree, degree.hat)
  
  all.mu.mse[i] <- mu.mse
  all.sigma.mse[i] <- sigma.mse
  all.degree.mse[i] <- degree.mse
  
  cat(sprintf("\r%d", i))
}
cat("\n")

hist(all.mu.mse, main = sprintf('T = %d', num.days), xlab = 'MSE of mu')
hist(all.sigma.mse, main = sprintf('T = %d', num.days), xlab = 'MSE of sigma')
hist(all.degree.mse, main = sprintf('T = %d', num.days), xlab = 'MSE of df')