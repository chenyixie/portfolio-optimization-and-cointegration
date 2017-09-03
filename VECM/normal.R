library(zoo)
library(R.matlab)
library(tsDyn)
library(vars)
library(expm)
library(Rdonlp2)


rm(list=ls())

# ===== Functions for fitting model =====

GetRank <- function(data) {
  num.train.days = dim(data)[1]
  
  var.est <- VAR(data, p=1, type="const", lag.max=12*(num.train.days/100)^0.25, ic="AIC")
  lag.length <- max(12*(num.train.days/100)^0.25, var.est$p)
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

# ===== Functions for estimation =====

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


GetFinalReturn <- function(init.price, final.price, weight) {
  # Get final return
  final.return <- matrix((final.price - init.price), nrow=1) %*% matrix(weight, ncol=1)
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

RunOneWindow <- function(all.data, gamma, num.train.days, num.test.days) {
  num.total.days <- num.train.days + num.test.days
  
  all.est.mean <- rep(0, 9)
  all.est.var <- matrix(0, nrow=9, ncol=9)
  for (g in 1:3) {
    data <- all.data[[g]]
    train.data <- data[1:num.train.days, ]
    
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
    
    test.init.price <- c(data[num.train.days, 1], data[num.train.days, 2], data[num.train.days, 3])
    est.mean <- EstimateMean(pi.hat, num.test.days, test.init.price)
    est.var <- EstimateVar(pi.hat, R.hat, num.test.days)
    
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
    all.test.init.price[(g * 3 - 2):(g * 3)] <- c(data[num.train.days, 1], data[num.train.days, 2], data[num.train.days, 3])
    all.test.final.price[(g * 3 - 2):(g * 3)] <- c(data[num.total.days, 1], data[num.total.days, 2], data[num.total.days, 3])
  }
  
  opt.final.return <- GetFinalReturn(all.test.init.price, all.test.final.price, opt.weight)
  return(list(opt.final.return,res))
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
    
    #cat(sprintf("%d~%d, %.4f\n", start, start+119, final.return))
    
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

test1 <- function(all.data) {
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
}

test2 <- function(all.data, num.train.days, num.test.days) {
  gamma <- 0.1
  num.total.days <- num.train.days + num.test.days
  
  all.data.slice <- list()
  for (g in 1:3) {
    all.data.slice[[g]] <- SliceData(all.data[[g]], 1, num.total.days)
  }
  
  result <- RunOneWindow(all.data.slice, gamma, num.train.days, num.test.days)
  final.return <- result
  
  print(final.return)
}

test3 <- function(all.data, num.train.days) {
  all.num.train.days <- c(100, 200, 300, 400, 500)
  num.test.days <- 20
  gamma <- 0.1
  
  
  for (num.train.days in all.num.train.days) {
    result <- RunAllWindows(all.data, gamma, num.train.days, num.test.days)
    all.final.return <- result
    return.mean <- mean(all.final.return)
    
    cat(sprintf("%d+%d: mean = %.5f, var = %.5f, sharpe = %.5f, maxdd = %.5f\n", num.train.days, num.test.days, return.mean))
  }
}

test4 <- function(all.data) {
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

test4(all.data)
# qq plot for one window for gamma= 0.1 (daily freq)
#x <- RunOneWindow(all.data, 0.1, 100, 5)
#qqnorm(x[[2]])
#qqline(x[[2]])
