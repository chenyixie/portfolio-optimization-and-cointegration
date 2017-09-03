library(zoo)
library(R.matlab)
library(tsDyn)
library(vars)
library(expm)
library(Rdonlp2)
library(mnormt)


rm(list=ls())


# ======== functions for fit and estimization =========

SliceData <- function(data, start, offset) {
  return(data[start:(start+offset-1), ])
}


FitModel <- function(all.data) {
  num.days <- dim(all.data)[1]
  num.assets <- dim(all.data)[2]
  dy <- diff(as.matrix(all.data))# 99 * 9
  # t.dy <- t(dy)   # 9*99
  mu.hat <- colMeans(dy)           # 9 * 1
  #sigma.hat = (t.dy-mu.hat) %*% t(t.dy-mu.hat)/99
  sigma.hat<- matrix(0, ncol = num.assets, nrow = num.assets)
  for(i in 1: num.assets){
    sigma.hat[i,i] <- var(dy[,i])
  }
  # sigma.hat = cov(dy)
  epsilon.hat <- rmnorm(num.days - 1, 0, sigma.hat) # 99 * 9
  dy.hat <- mu.hat + epsilon.hat    # 99 * 9
  residual <- dy - dy.hat    # 99*9
  
  return(list(mu.hat, sigma.hat, residual))
}

EstimateMean <- function(mu.hat, num.test.days) {
  mean.est <- mu.hat * num.test.days   # 9 * 1
  return(mean.est)
}


EstimateVar <- function(sigma.hat, num.test.days) {
  var.est <- sigma.hat * num.test.days
  return(var.est)    # 9*9
}


GetRealMean <- function(test.init.price, test.final.price) {
  return(test.final.price - test.init.price)
}


# ======== functions for optimization =========

GetPortExp <- function(weight, assets.mean) {
  port.mean = matrix(weight, nrow=1) %*% matrix(assets.mean, ncol=1)
  return(port.mean)
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
  final.return <- matrix(final.price-init.price, nrow=1) %*% matrix(weight, ncol=1)
  return(c(final.return))
}


GetRsquare <- function(train.data, res) {
  # train.data: data of 9 assets (matrix) (100 * 9)
  # res: 99 * 9
  
  dy <- diff(as.matrix(train.data))
  sst = c()
  total.sst =0
  for(i in 1:9){
    sst = sum((dy[,i] - mean(dy[,i])) ^ 2)
    total.sst = sum(sst[i])
  }
  sse = c()
  total.sse = 0
  for(i in 1:9){
    sse[i] = sum(res[,i]^2)
    total.sse = sum(sse[i])
  }
  # dy <- diff(as.matrix(train.data))    # 99*9
  # total.sst <- sum((t(t(dy) - colMeans(dy)))^2)  # sum(99 * 9)
  # total.sse = sum(res ^ 2)    # sum(99*9)
  r.square = 1 - total.sse / total.sst
  return(r.square)
}



# ===== functions for experiment =====

RunOneWindow <- function(all.data, gamma, num.train.days, num.test.days) {
  num.assets <- dim(all.data)[2]
  train.data <- all.data[1:num.train.days, ]
  fit.result <- FitModel(train.data)
  # ---------------------------------
  
  mu.hat <- fit.result[[1]]
  sigma.hat <- fit.result[[2]]
  res <- fit.result[[3]]
  rsquare <- GetRsquare(train.data, res)
  est.mean <- EstimateMean(mu.hat, num.test.days)
  est.var <- EstimateVar(sigma.hat, num.test.days)  # matrix 9*9
  
  # ----------------------------------
  opt.weight <- Optimize(est.mean, est.var, gamma)
  
  all.test.init.price <- rep(0, num.assets)
  all.test.final.price <- rep(0, num.assets)
  for (i in 1:num.assets) {
    data <- all.data[, i]
    all.test.init.price[i] <- data[num.train.days]
    all.test.final.price[i] <- data[num.train.days + num.test.days]
  }
  
  opt.final.return <- GetFinalReturn(all.test.init.price, all.test.final.price, opt.weight)
  return(list(opt.final.return, res))
}


RunAllWindows <- function(all.data, gamma, num.train.days, num.test.days) {
  num.one.window.days = num.train.days + num.test.days
  
  all.final.return <- c()
  all.rsquare <- c()
  start <- 1
  num.windows <- 1
  while (start + num.one.window.days <= dim(all.data)[1]) {
    all.data.slice <- SliceData(all.data, start, num.one.window.days)
    
    run.result <- RunOneWindow(all.data.slice, gamma, num.train.days, num.test.days)
    final.return <- run.result[[1]]
    rsquare <- run.result[[2]]  # a single value
    
    all.final.return[num.windows] <- final.return
    all.rsquare[num.windows] <- rsquare
    
    #cat(sprintf("%d~%d, %.4f\n", start, start+119, final.return))
    
    start <- start + num.test.days
    num.windows <- num.windows + 1
  }
  
  return(list(all.final.return, all.rsquare))
}

#result.it <- RunAllWindows(data.it, 0.1)
#all.rsquare.it <- result.it[[2]]
#hist(all.rsquare.it, main = '', xlab = 'R-Squared (IT)' )

RunAllGammas <- function(all.data, all.gamma) {
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
    
    plot(all.cum.returns, type='l')
  }
  
  plot(all.means, type='l')
  plot(all.vars, type='l')
  
  return(all.results)
}


TestOneWindow <- function(all.data, num.train.days, num.test.days) {
  num.total.days <- num.train.days + num.test.days
  all.data.slice <- SliceData(all.data, 1, num.total.days)
  
  for (i in 1:dim(all.data.slice)[2]) {
    data.slice <- all.data.slice[, i]
    
    fit.result <- FitModel(data.slice)
    a.hat <- fit.result[[1]]
    epsilon.hat <- fit.result[[2]]
    
    real.mean <- data.slice[num.total.days] - data.slice[num.train.days]
    esti.mean <- EstimateMean(a.hat, num.test.days)
    
    plot(data.slice, type='l')
    cat(sprintf("a.hat = %f\n", a.hat))
    cat(sprintf("epsilon.hat = %f\n", epsilon.hat))
    cat(sprintf("real.mean = %f\n", real.mean))
    cat(sprintf("esti.mean = %f\n", esti.mean))
    cat("\n")
  }
}



# ===== Functions for fetch data of different frequency =====

GetWeeklyData <- function(all.data) {
  num.total.days <- dim(all.data)[1]
  
  week.index <- seq(5, num.total.days, 5)
  all.data.weekly <- all.data[week.index, ]
  
  return(all.data.weekly)
}

GetMonthlyData <- function(all.data) {
  num.total.days <- dim(all.data)[1]
  
  month.index <- seq(20, num.total.days, 20)
  all.data.monthly <- all.data[month.index, ]
  
  return(all.data.monthly)
}

# ===== Do experiment =====
load("data_for_exp_real.RData")

all.data <- data.frame(app.price, intel.price, msft.price,
                       lb.price, pvh.price, tjx.price,
                       cat.price, eaton.price, north.price)

test6 <- function(all.data) {
  all.data.weekly <- GetWeeklyData(all.data)
  all.data.monthly <- GetMonthlyData(all.data)
  
  all.gamma <- c(0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.01, 0.1, 1, 5, 10)
  all.results <- list()
  all.results.week <-list()
  all.results.month <-list()
  
  # gamma <- 0.001
  num.train.days <- 100
  num.test.days <- 5
  
  cat("gamma.daily, mean, total mean, sd, sharpe, maxdd\n")
  for (i in 1:length(all.gamma)) {
    gamma <- all.gamma[[i]]
    result <- RunAllWindows(all.data, gamma, num.train.days, num.test.days)
    all.results[[i]] <- result
    all.final.return <- result[[1]]
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
    all.results.week[[i]] <- result.week
    all.final.return.week <- result.week[[1]]
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
    all.final.return.month <- result.month[[1]]
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

test6(all.data)












