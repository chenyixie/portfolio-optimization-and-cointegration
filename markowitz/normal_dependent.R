# Markowitz with Student t Distr.
# 2017.8.4

library(timeSeries)

# ===== Functions for fitting model =====

FitModel <- function(all.data) {
  num.days <- dim(all.data)[1]
  num.assets <- dim(all.data)[2]
  
  dy <- diff(as.matrix(all.data))   # 99 * 9
  t.dy <- t(dy)                     # 9 * 99
  mu.hat <- colMeans(dy)            # 9 * 1
  sigma.hat = cov(dy)

  return(list(mu=mu.hat, sigma=sigma.hat))
}


# ===== Functions for estimation =====

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


# ===== Functions for generating synthetic data ====

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
