library(zoo)
library(R.matlab)
library(tsDyn)
library(vars)
library(expm)
library(Rdonlp2)
library(timeSeries)
library(QRM)
library(mvtnorm)

rm(list=ls())

setwd("~/Desktop/final project/codes")

# ======== functions for fit and estimization =========

SliceData <- function(data, start, offset) {
  return(data[start:(start+offset-1), ])
}

FitVecm <- function(data) {
  fit         <- VECM(data, r=1, lag=0, estim="ML", include="none")
  fit.summary <- summary(fit)
  fit.alpha   <- fit.summary$coefficients
  fit.beta    <- fit.summary$model.specific$beta
  fit.pi      <- fit.alpha %*% t(fit.beta)
  data.new    <- matrix(0, ncol = dim(data)[2], nrow = dim(data)[1])
  data.new <- data
  data.new[1,] <- data[1,]
  for(i in 2:dim(data)[1]){
    data.new[i,] <- data[i,] - data[(i-1),] - as.matrix(data[(i-1),]) %*% fit.pi
  }
  fit.t <- fit.mst(data.new, method = "BFGS")
  fit.R <- fit.t$Sigma
  fit.df<- fit.t$df
  fit.result  <- list(fit.pi, fit.R, fit.df)
  
  return(fit.result)
}

EstimateMean <- function(pi, num.test.days, init.price) {
  identity <- diag(3)
  mean.est <- (identity + pi) %^% (num.test.days) %*% init.price -init.price
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


GetRealMean <- function(test.init.price, test.final.price) {
  return(test.final.price - test.init.price)
}


# ======== functions for optimization =========

GetPortExp <- function(weight, asset.exp) {
  port.exp = matrix(weight, nrow=1) %*% matrix(asset.exp, ncol=1)
  return(port.exp)
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
  # Get final return
  final.return <- matrix((final.price - init.price), nrow=1) %*% matrix(weight, ncol=1)
  return(c(final.return))
}

# ===== functions for experiment =====

RunOneWindow <- function(all.data, gamma, num.train.days, num.test.days) {
  num.total.days <- num.train.days + num.test.days
  
  all.est.mean <- rep(0, 9)
  all.est.var <- matrix(0, nrow=9, ncol=9)
  for (g in 1:3) {
    data <- all.data[[g]]
    train.data <- data[1:num.train.days, ]
    
    fit.result <- FitVecm(train.data)
    pi.hat <- fit.result[[1]]
    R.hat <- fit.result[[2]]
    df.hat <- fit.result[[3]]
    
    test.init.price <- c(data[num.train.days, 1], data[num.train.days, 2], data[num.train.days, 3])
    est.mean <- EstimateMean(pi.hat, num.test.days, test.init.price)
    est.var <- EstimateVar(pi.hat, R.hat,df.hat, num.test.days)
    est.var <- matrix(est.var, ncol = 3)
    
    all.est.mean[(g * 3 - 2):(g * 3)] <- est.mean
    all.est.var[(g * 3 - 2):(g * 3), (g * 3 - 2):(g * 3)] <- est.var
  }
  
  opt.weight <- Optimize(all.est.mean, all.est.var, gamma)
  
  all.test.init.price <- rep(0, 9)
  all.test.final.price <- rep(0, 9)
  for (g in 1:3) {
    data <- all.data[[g]]
    all.test.init.price[(g * 3 - 2):(g * 3)] <- c(data[num.train.days, 1], data[num.train.days, 2], data[num.train.days, 3])
    all.test.final.price[(g * 3 - 2):(g * 3)] <- c(data[num.total.days, 1], data[num.total.days, 2], data[num.total.days, 3])
  }
  
  opt.final.return <- GetFinalReturn(all.test.init.price, all.test.final.price, opt.weight)
  return(opt.final.return)
}


RunAllWindows <- function(all.data, gamma, num.train.days, num.test.days) {
  num.total.days <- num.train.days + num.test.days
  
  all.final.return <- c()
  start <- 1
  num.windows <- 1
  while (start + num.total.days <= dim(all.data[[1]])[1]) {
    all.data.slice <- list()
    for (g in 1:3) {
      all.data.slice[[g]] <- SliceData(all.data[[g]], start, num.total.days)
    }
    run.result <- RunOneWindow(all.data.slice, gamma, num.train.days, num.test.days)
    final.return <- run.result
    all.final.return[num.windows] <- final.return
    
    start <- start + num.test.days
    num.windows <- num.windows + 1
  }
  
  return(all.final.return)
}

# ===== Functions for fetch data of different frequency =====

GetWeeklyData <- function(all.data) {
  data.it <- all.data[[1]]
  data.cd <- all.data[[2]]
  data.ind <- all.data[[3]]
  
  num.total.days <- dim(data.it)[1]
  week.index <- seq(5, num.total.days, 5)
  data.it.weekly <- data.it[week.index, ]
  data.cd.weekly <- data.cd[week.index, ]
  data.ind.weekly <- data.ind[week.index, ]
  
  all.data.weekly <- list(data.it.weekly, data.cd.weekly, data.ind.weekly)
  return(all.data.weekly)
}

GetMonthlyData <- function(all.data) {
  data.it <- all.data[[1]]
  data.cd <- all.data[[2]]
  data.ind <- all.data[[3]]
  
  num.total.days <- dim(data.it)[1]
  
  monthly.index <- seq(20, num.total.days, 20)
  data.it.monthly <- data.it[monthly.index, ]
  data.cd.monthly <- data.cd[monthly.index, ]
  data.ind.monthly <- data.ind[monthly.index, ]
  
  all.data.monthly <- list(data.it.monthly, data.cd.monthly, data.ind.monthly)
  return(all.data.monthly)
}


# ===== Do experiment =====

load("data_for_exp_real.RData")
data.it <- data.frame(app.price, intel.price, msft.price)
data.cd <- data.frame(lb.price, pvh.price, tjx.price)
data.ind <- data.frame(cat.price, eaton.price, north.price)
all.data <- list(data.it, data.cd, data.ind)

test <- function(all.data) {
  all.data.weekly <- GetWeeklyData(all.data)
  all.data.monthly <- GetMonthlyData(all.data)
  
  all.gamma <- c(0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.01, 0.1, 1, 5, 10)
  all.results <- list()
  all.results.week <-list()
  all.results.month <-list()
  
  num.train.days <- 100
  num.test.days <- 5
  
  cat("gamma.daily, mean, total mean, sd, sharpe, maxdd\n")
  for (i in 1:length(all.gamma)) {
    gamma <- all.gamma[[i]]
    result <- RunAllWindows(all.data, gamma, num.train.days, num.test.days)
    all.results[[i]] <- result
    all.final.return <- result
    return.mean <- mean(all.final.return)
    return.sum <- sum(all.final.return)
    return.sd <- sd(all.final.return)
    all.cum.returns <- cumsum(all.final.return)
    sharpe.ratio <- sharpe(all.cum.returns)
    max.drawdown <- maxdrawdown(all.cum.returns)
    dd.value <- max.drawdown$maxdrawdown
    dd.price <- all.cum.returns[max.drawdown$from]
    dd.ratio <- dd.value / dd.price
    cat(sprintf("%f,%.5f, %.5f, %.5f, %.5f, %.5f\n",gamma, return.mean, return.sum, return.sd, sharpe.ratio, dd.ratio))
  }
  
  cat("gamma.week, mean, total mean, sd, sharpe, maxdd\n")
  for (i in 1:length(all.gamma)) {
    gamma <- all.gamma[[i]]
    result.week <- RunAllWindows(all.data.weekly, gamma, num.train.days, num.test.days)
    all.final.return.week <- result.week
    all.results.week[[i]] <- result.week
    return.mean.week <- mean(all.final.return.week)
    return.sum.week <- sum(all.final.return.week)
    return.sd.week <- sd(all.final.return.week)
    all.cum.returns.week <- cumsum(all.final.return.week)
    sharpe.ratio.week <- sharpe(all.cum.returns.week)
    max.drawdown.week <- maxdrawdown(all.cum.returns.week)
    dd.value.week <- max.drawdown.week$maxdrawdown
    dd.price.week <- all.cum.returns.week[max.drawdown.week$from]
    dd.ratio.week <- dd.value.week / dd.price.week
    cat(sprintf("%f,%.5f, %.5f, %.5f, %.5f, %.5f\n",gamma, return.mean.week,return.sum.week, return.sd.week, sharpe.ratio.week, dd.ratio.week))
  }
  
  cat("gamma.month, mean, total mean, sd, sharpe, maxdd\n")
  for (i in 1:length(all.gamma)) {
    gamma <- all.gamma[[i]]
    result.month <- RunAllWindows(all.data.monthly, gamma, num.train.days, num.test.days)
    all.results.month[[i]] <- result.month
    all.final.return.month <- result.month
    return.mean.month <- mean(all.final.return.month)
    return.sum.month <- sum(all.final.return.month)
    return.sd.month <- sd(all.final.return.month)
    all.cum.returns.month <- cumsum(all.final.return.month)
    sharpe.ratio.month <- sharpe(all.cum.returns.month)
    max.drawdown.month <- maxdrawdown(all.cum.returns.month)
    dd.value.month <- max.drawdown.month$maxdrawdown
    dd.price.month <- all.cum.returns.month[max.drawdown.month$from]
    dd.ratio.month <- dd.value.month / dd.price.month
    cat(sprintf("%f,%.5f, %.5f, %.5f, %.5f, %.5f\n",gamma, return.mean.month,return.sum.month, return.sd.month, sharpe.ratio.month, dd.ratio.month))
  }
}

test(all.data)