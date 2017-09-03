# Markowitz with Student t Distr.
# 2017.8.4

library(timeSeries)

# ===== Functions for fitting model =====

FitOneAsset <- function(data.one.asset) {
  dx <- diff(data.one.asset)
  return(list(mu=mean(dx), sigma=var(dx)))
}


FitModel <- function(all.data) {
  num.assets = dim(all.data)[2]
  mu.hat <- rep(0, num.assets)
  sigma.hat <- matrix(0, nrow=num.assets, ncol=num.assets)
  
  for (i in 1:num.assets) {
    fit.result <- FitOneAsset(all.data[, i])
    mu.hat[i] <- fit.result$mu
    sigma.hat[i, i] <- fit.result$sigma
  }
  
  return(list(mu=mu.hat, sigma=sigma.hat))
}

EstimateMean <- function(fit.result, num.test.days) {
  mu.hat <- fit.result$mu
  mean.est <- mu.hat * num.test.days   # 9 * 1
  return(mean.est)
}

EstimateVar <- function(fit.result, num.test.days) {
  sigma.hat <- fit.result$sigma
  var.est <- sigma.hat * num.test.days
  return(var.est)    # 9*9
}


# ===== Functions for generate synthetic data =====

GenerateData <- function(mu, sigma, init.price, num.days) {
  num.assets <- length(init.price)
  all.data <- matrix(0, nrow=num.days, ncol=num.assets)
  all.data[1, ] <- init.price
  
  # dx ~ N(mu, sigma)
  dx <- mvrnorm(n=(num.days-1), mu=mu, Sigma=sigma)
  for (day in 2:num.days) {
    all.data[day, ] <- all.data[day-1, ] + dx[day-1, ]
  }
  
  return(all.data)
}
