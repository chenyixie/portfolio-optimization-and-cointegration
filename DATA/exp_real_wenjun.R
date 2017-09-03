# Exp real (wenjun) - 2017.7.27
#
# 1. Get 100 days from Group1(IT):	app.price, intel.price, msft.price
# 2. Get rank of this data slice -> rank
# 3. Fit VECM use the data slice and rank -> pi.hat, R.hat
# 4. Use pi.hat and R.hat to estimate mean and variance of future 20 days
# 5. Optimize

library(zoo)
library(R.matlab)
library(tsDyn)
library(vars)
library(expm)
library(Rdonlp2)


rm(list=ls())


# ======== functions =========

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

  # lower bound and upper bound
  par.l = rep(0, 3)
  par.u = rep(0.5, 3)

  # constraints sum of weight = 1
  A = matrix(rep(1, 3), nrow=1)
  lin.l = 1
  lin.u = 1

  # initial weight
  init.weight <- rep(1/3, 3)

  # do optimization
  ret = donlp2(init.weight, ObjectiveFunction, par.u=par.u, par.l=par.l, A, lin.l=lin.l, lin.u=lin.u)
  opt.weight = ret$par

  return(opt.weight)
}


GetFinalReturn <- function(init.price, final.price, weight) {
  # Get final return
  final.return <- matrix((final.price / init.price - 1), nrow=1) %*% matrix(weight, ncol=1)
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


RunOneWindow <- function(data, gamma) {
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
  mean.est <- EstimateMean(pi.hat, 20, test.init.price)
  var.est <- EstimateVar(pi.hat, R.hat, 20)

  opt.weight <- Optimize(mean.est, var.est, gamma)

  test.final.price <- c(data[120, 1], data[120, 2], data[120, 3])
  opt.final.return <- GetFinalReturn(test.init.price, test.final.price, opt.weight)

  return(list(opt.final.return, rsquare))
}


RunAllWindows <- function(data, gamma) {
  all.final.return <- c()
  all.rsquare <- c()
  start <- 1
  num.windows <- 1
  while (start + 120 <= 6000) {
    data.slice <- SliceData(data, start, 120)
    run.result <- RunOneWindow(data.slice, gamma=0.1)
    final.return <- run.result[[1]]
    rsquare <- run.result[[2]]
    all.final.return[num.windows] <- final.return
    all.rsquare[num.windows] <- rsquare

    cat(sprintf("%d~%d, %.4f, %.4f\n", start, start+119, final.return, rsquare))

    start <- start + 20
    num.windows <- num.windows + 1
  }

  return(list(all.final.return, all.rsquare))
}


# ===== Do experiment =====

load("data_for_exp_real.RData")

data.it <- data.frame(app.price, intel.price, msft.price)
data.cd <- data.frame(lb.price, pvh.price, tjx.price)
data.ind <- data.frame(cat.price, eaton.price, north.price)

result.it <- RunAllWindows(data.it, 0.1)
all.final.returns.it <- result.it[[1]]
all.rsquare.it <- result.it[[2]]
all.cum.returns.it <- cumsum(all.final.returns.it)
plot(all.cum.returns.it, type='l')
plot(all.rsquare.it, type='l')
hist(all.rsquare.it)

result.cd <- RunAllWindows(data.cd, 0.1)
all.final.returns.cd <- result.cd[[1]]
all.rsquare.cd <- result.cd[[2]]
all.cum.returns.cd <- cumsum(all.final.returns.cd)
plot(all.cum.returns.cd, type='l')
plot(all.rsquare.cd, type='l')
hist(all.rsquare.cd)

result.ind <- RunAllWindows(data.ind, 0.1)
all.final.returns.ind <- result.ind[[1]]
all.rsquare.ind <- result.ind[[2]]
all.cum.returns.ind <- cumsum(all.final.returns.ind)
plot(all.cum.returns.ind, type='l')
plot(all.rsquare.ind, type='l')
hist(all.rsquare.ind)

plot((all.cum.returns.it + all.cum.returns.cd + all.cum.returns.ind) / 3, type='l')
