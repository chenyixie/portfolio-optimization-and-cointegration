# Shared functions for all models
# 2017.8.5

# ===== Functions for preparing data =====

GetAllData <- function() {
  load("data_for_exp_real.RData")
  all.data <- cbind(app.price, intel.price, msft.price,
                    lb.price, pvh.price, tjx.price,
                    cat.price, eaton.price, north.price)
  return(all.data)    # matrix 6000 * 9
}


SliceData <- function(data, start, offset) {
  return(data[start:(start+offset-1), ])
}


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


# ===== Functions for optimization =====

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


# ===== Functions for getting final return =====

GetFinalReturn <- function(init.price, final.price, weight) {
  # Get final return
  final.return <- matrix(final.price-init.price, nrow=1) %*% matrix(weight, ncol=1)
  return(c(final.return))
}


# ===== Functions for rolling windows =====

RunOneWindow <- function(all.data, gamma, num.train.days, num.test.days) {
  num.assets <- dim(all.data)[2]
  train.data <- all.data[1:num.train.days, ]
  fit.result <- FitModel(train.data)

  est.mean <- EstimateMean(fit.result, num.test.days)
  est.var <- EstimateVar(fit.result, num.test.days)  # matrix 9*9
  
  opt.weight <- Optimize(est.mean, est.var, gamma)
  
  all.test.init.price <- rep(0, num.assets)
  all.test.final.price <- rep(0, num.assets)
  for (i in 1:num.assets) {
    data <- all.data[, i]
    all.test.init.price[i] <- data[num.train.days]
    all.test.final.price[i] <- data[num.train.days + num.test.days]
  }
  
  opt.final.return <- GetFinalReturn(all.test.init.price, all.test.final.price, opt.weight)
  return(list(opt.final.return))
}


RunAllWindows <- function(all.data, gamma, num.train.days, num.test.days) {
  num.one.window.days = num.train.days + num.test.days
  
  all.final.return <- c()
  start <- 1
  num.windows <- 1
  while (start + num.one.window.days <= dim(all.data)[1]) {    
    cat(sprintf("\r%d: %d~%d", num.windows, start, start + num.one.window.days - 1))
    
    all.data.slice <- SliceData(all.data, start, num.one.window.days)
    
    run.result <- RunOneWindow(all.data.slice, gamma, num.train.days, num.test.days)
    final.return <- run.result[[1]]
    all.final.return[num.windows] <- final.return

    start <- start + num.test.days
    num.windows <- num.windows + 1
  }
  cat("\r")
  
  return(list(all.final.return))
}


RunAllGammas <- function(all.data, all.gamma, num.train.days, num.test.days) {
  cat("gamma, mean, total, sd, sharpe, maxdd\n")
  
  for (i in 1:length(all.gamma)) {
    result <- RunAllWindows(all.data, all.gamma[i], num.train.days, num.test.days)
    all.final.returns <- result[[1]]
    
    performance <- GetPerformance(all.final.returns)
    return.mean <- performance[[1]]
    total.return <- performance[[2]]
    return.sd <- performance[[3]]
    sharpe.ratio <- performance[[4]]
    dd.ratio <- performance[[5]]
    
    cat(sprintf("%f, %.5f, %.5f, %.5f, %.5f, %.5f\n", all.gamma[i], return.mean, total.return, return.sd, sharpe.ratio, dd.ratio))
  }
}


# ===== Functions for performance =====

GetPerformance <- function(all.final.returns) {
  return.mean <- mean(all.final.returns)
  return.sd <- sd(all.final.returns)
  total.return <- sum(all.final.returns)
  
  all.cum.returns <- cumsum(all.final.returns)
  sharpe.ratio <- sharpe(all.cum.returns)
  
  max.drawdown <- maxdrawdown(all.cum.returns)
  dd.value <- max.drawdown$maxdrawdown
  dd.price <- all.cum.returns[max.drawdown$from]
  dd.ratio <- dd.value / dd.price
  
  return(list(return.mean, total.return, return.sd, sharpe.ratio, dd.ratio))
}


GetMSE <- function(vector1, vector2) {
  vector1 <- as.vector(vector1)
  vector2 <- as.vector(vector2)
  
  return(mean((vector1 - vector2) ^ 2))
}


GetMAPE <- function(real_vector, forcast_vector) {
  real_vector <- as.vector(real_vector)
  forcast_vector <- as.vector(forcast_vector)
  
  n <- length(real_vector)
  mape <- 100 / n * sum(abs((real_vector - forcast_vector) / real_vector))
  
  return(mape)
}

