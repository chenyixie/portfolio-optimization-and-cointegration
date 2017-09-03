# 30/07/2017
# Reward Function with Transaction Cost 
library(zoo)
library(R.matlab)
library(tsDyn)
library(vars)
library(expm)
library(Rdonlp2)


rm(list=ls())


# ======== functions for fit and estimization =========

SliceData <- function(data, start, offset) {
  return(data[start:(start+offset-1), ])
}


GetRank <- function(data) {
  var.est <- VAR(data, p=1, type="const", lag.max=12*(100/100)^0.25, ic="AIC")
  lag.length <- max(12*(100/100)^0.25, var.est$p)
  cointest <- ca.jo(data, K= lag.length, type = "trace", ecdet = "const", spec = "transitory")
  teststat <- cointest@teststat
  cval <- cointest@cval
  
  rank <- 3
  for (r in 0:2) {
    if (teststat[3-r] <= cval[6-r]) {
      rank <- r
      break
    }
  }
  return (rank)
}


FitVecm <- function(data, rank) {
  fit         <- VECM(data, r=rank, lag=0, estim="ML", include="none")
  fit.summary <- summary(fit)
  fit.R       <- cov(fit.summary$residuals)
  fit.alpha   <- fit.summary$coefficients
  fit.beta    <- fit.summary$model.specific$beta
  fit.pi      <- fit.alpha %*% t(fit.beta)
  fit.result  <- list(fit.pi, fit.R, fit.summary$residuals)
  
  return(fit.result)
}


EstimateMean <- function(pi, num.test.days, init.price) {
  identity <- diag(3)
  mean.est <- (identity + pi) %^% (num.test.days) %*% init.price - init.price
  return(mean.est)
}


EstimateVar <- function(pi, R, num.test.days) {
  identity = diag(3)
  var.est <- matrix(0, nrow=3, ncol=3)
  
  for (k in 0:(num.test.days - 1)) {
    var.est = var.est + t((identity + pi) %^% (k)) %*% R %*% ((identity + pi) %^% (k))
  }
  
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
  final.return <- final.return - abs(final.return) * 0.005
  return(c(final.return))
}


GetRsquare <- function(train.data, res) {
  # SST
  sst = c()
  total.sst =0
  for(i in 1:3){
    sst[i] = sum((train.data[, i]-mean(train.data[, i]))^2)
    total.sst = sum(sst[i])
  }
  
  # SSE
  sse = c()
  total.sse = 0
  for (i in 1:3){
    sse[i] = sum(res[,i]^2)
    total.sse = sum(sse[i])
  }
  
  # R square
  r.square = 1 - total.sse / total.sst
  return(r.square)
}


# ===== functions for experiment =====

RunOneWindow <- function(all.data, gamma) {
  all.est.mean <- rep(0, 9)
  all.est.var <- matrix(0, nrow=9, ncol=9)
  for (g in 1:3) {
    data <- all.data[[g]]
    train.data <- data[1:100, ]
    
    rank <- GetRank(train.data)
    if (rank == 3) {
      rank <- 2
    }
    if (rank == 0) {
      rank <- 1
    }
    
    fit.result <- FitVecm(train.data, rank)
    pi.hat <- fit.result[[1]]
    R.hat <- fit.result[[2]]
    res <- fit.result[[3]]
    rsquare <- GetRsquare(train.data, res)
    
    test.init.price <- c(data[100, 1], data[100, 2], data[100, 3])
    est.mean <- EstimateMean(pi.hat, 20, test.init.price)
    est.var <- EstimateVar(pi.hat, R.hat, 20)
    
    all.est.mean[(g * 3 - 2):(g * 3)] <- est.mean
    all.est.var[(g * 3 - 2):(g * 3), (g * 3 - 2):(g * 3)] <- est.var
  }
  
  opt.weight <- Optimize(all.est.mean, all.est.var, gamma)
  
  PrintWeight <- function(weight) {
    for (w in weight) {
      cat(sprintf('%.3f ', w))
    }
    cat("\n")
  }
  
  all.test.init.price <- rep(0, 9)
  all.test.final.price <- rep(0, 9)
  for (g in 1:3) {
    data <- all.data[[g]]
    all.test.init.price[(g * 3 - 2):(g * 3)] <- c(data[100, 1], data[100, 2], data[100, 3])
    all.test.final.price[(g * 3 - 2):(g * 3)] <- c(data[120, 1], data[120, 2], data[120, 3])
  }
  
  opt.final.return <- GetFinalReturn(all.test.init.price, all.test.final.price, opt.weight)
  return(opt.final.return)
}


RunAllWindows <- function(all.data, gamma) {
  all.final.return <- c()
  start <- 1
  num.windows <- 1
  while (start + 120 <= 6000) {
    all.data.slice <- list()
    for (g in 1:3) {
      all.data.slice[[g]] <- SliceData(all.data[[g]], start, 120)
    }
    run.result <- RunOneWindow(all.data.slice, gamma)
    final.return <- run.result
    all.final.return[num.windows] <- final.return
    
    #cat(sprintf("%d~%d, %.4f\n", start, start+119, final.return))
    
    start <- start + 20
    num.windows <- num.windows + 1
  }
  
  return(all.final.return)
}


# ===== Do experiment =====

load("data_for_exp_real.RData")

data.it <- data.frame(app.price, intel.price, msft.price)
data.cd <- data.frame(lb.price, pvh.price, tjx.price)
data.ind <- data.frame(cat.price, eaton.price, north.price)
all.data <- list(data.it, data.cd, data.ind)

#all.gamma <- c(0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 10, 50, 100)
all.gamma <- c(0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.01, 0.1, 1, 5, 10)
all.results <- list()
all.means <- rep(0, length(all.gamma))
all.vars <- rep(0, length(all.gamma))

cat("gamma, mean, var, shapre, maxdd\n")
for (i in 1:length(all.gamma)) {
  gamma <- all.gamma[[i]]
  result <- RunAllWindows(all.data, gamma)
  all.results[[i]] <- result
  
  all.final.returns <- result
  return.mean <- mean(all.final.returns)
  return.var <- var(all.final.returns)
  all.cum.returns <- cumsum(all.final.returns)
  sharpe.ratio <- sharpe(all.cum.returns)
  max.drawdown <- maxdrawdown(all.cum.returns)
  dd.value <- max.drawdown$maxdrawdown
  dd.price <- all.cum.returns[max.drawdown$from]
  dd.ratio <- dd.value / dd.price
  cat(sprintf("%f, %.5f, %.5f, %.5f, %.5f\n", gamma, return.mean, return.var, sharpe.ratio, dd.ratio))
  all.means[i] <- return.mean
  all.vars[i] <- return.var
  
  plot(all.final.returns, type = 'l', main = 'Final Returns', xlab = 'number of experiment', ylab = 'Final Rewards($)')
  plot(all.cum.returns, type='l', main = 'Cumulative Return',xlab = 'number of experiment', ylab = 'Final Rewards($)')
}

final.return.gamma2 <-RunAllWindows(all.data, 0.1)
all.cum.returns <- cumsum(final.return.gamma2)
plot(all.cum.returns, type = 'l')
plot(final.return.gamma2, type = 'l',ylab = 'final return')
hist(final.return.gamma2)

plot(all.means, type='l')
plot(all.vars, type='l')