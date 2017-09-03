# Markowitz with Student t Distr.
# 2017.8.4

library(timeSeries)
library(QRM)
library(mvtnorm)

# ===== Functions for fitting model =====

FitModel <- function(all.data) {
  num.train.days <- dim(all.data)[1]
  num.assets <- dim(all.data)[2]
  dx <- diff(all.data)
  
  fit.result <- fit.mst(dx, method="BFGS")
  degree.hat <- fit.result$df
  mu.hat <- fit.result$mu
  sigma.hat <- fit.result$Sigma
  
  return(list(degree=degree.hat, mu=mu.hat, sigma=sigma.hat))
}


# ===== Functions for estimation =====

EstimateMean <- function(fit.result, num.test.days) {
  mu.hat <- fit.result$mu
  mean.est <- mu.hat * num.test.days   # 9 * 1
  return(mean.est)
}


EstimateVar <- function(fit.result, num.test.days) {
  sigma.hat <- fit.result$sigma
  degree.hat <- fit.result$degree
  
  if (degree.hat < 2) {
    degree.hat = 2.1
    }
  
  var.est <- sigma.hat * num.test.days * degree.hat / (degree.hat - 2)
  var.est.diag <- diag(var.est)
  
  var.est <- var.est * 0
  for (i in 1:dim(var.est)[1]) {
    var.est[i, i] <- var.est.diag[i]
  }
  
  return(var.est)    # 9*9
}


# ===== Functions for generating synthetic data ====

GenerateData <- function(degree, mu, sigma, init.price, num.days) {
  num.assets <- length(init.price)
  all.data <- matrix(0, nrow=num.days, ncol=num.assets)
  all.data[1, ] <- init.price
  
  # dx ~ mu + t(0, sigma, degree)
  dx <- mu + rmvt(n=(num.days-1), sigma=sigma, df=degree)
  for (day in 2:num.days) {
    all.data[day, ] <- all.data[day-1, ] + dx[day-1, ]
  }
  
  return(all.data)
}