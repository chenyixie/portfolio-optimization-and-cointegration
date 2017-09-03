# Exp real for 3 groups using simple model - 2017.7.28
#
# 1. Get 100 days from Group1
# 2. Fit data by formula:
#
#        d(y_t) = a + e_t
#        a.hat = mean(d(y_t))
#        mean(d(y_t+s - y_t)) = a * s
#        var(d(y_t+s - y_t)) = sigma^2 * s
#
# 3. Use a.hat as the estimated mean of future 20 days (var is 1)????
# 4. Do step 1~3 for Group2 and Group3
# 5. Optimize

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


FitModel <- function(data.one.asset) {
  # dx_t = a + e_t
  # dx_t.hat = a.hat - e_t.hat 
  # res = dx_t - dx_t.hat
  
  # data.one.asset : 100 values
  # dx: 99
  # epsilon.hat: 99 random values
  # a.hat: 1
  # dx.hat: a.hat + epsilon.hat 99
  # residual: dx - dx.hat 99
  
  dx <- diff(data.one.asset)
  a.hat <- mean(dx)
  sigma.hat.square <- sum((dx - a.hat)^2)/length(dx)
  sigma.hat <-sqrt(sigma.hat.square)
  epsilon.hat <- rnorm(99, 0, sigma.hat)
  
  dx.hat <- a.hat + epsilon.hat
  residual <- dx - dx.hat
  
  #plot(dx, type='l')
  #lines(rep(a.hat, 99), col='green')
  #lines(dx.hat, col='red')
  
  #plot(residual)
  #abline(c(0,0), col = 2)
  
  return(list(a.hat, sigma.hat, residual))
}


EstimateMean <- function(a, num.test.days) {
  mean.est <- a * num.test.days
  return(mean.est)
}


EstimateVar <- function(sigma.hat, num.test.days) {
  var.est <- sigma.hat * num.test.days
  return(var.est)
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
  # train.data: data of one asset (vector)
  sst = sum((train.data - mean(train.data)) ^ 2)
  sse = sum(res^2)
  r.square = 1 - sse / sst
  return(r.square)
}



# ===== functions for experiment =====

RunOneWindow <- function(all.data, gamma) {
  num.assets <- dim(all.data)[2]
  num.train.days <- 100
  num.test.days <- 20

  all.est.mean <- rep(0, num.assets)
  all.est.var <- matrix(0, nrow=num.assets, ncol=num.assets)
  all.rsquare <- rep(0, num.assets)
  
  for (i in 1:num.assets) {
    data <- all.data[, i]
    train.data <- data[1:num.train.days]

    fit.result <- FitModel(train.data)
    a.hat <- fit.result[[1]]
    sigma.hat <- fit.result[[2]]
    res <- fit.result[[3]]
    rsquare <- GetRsquare(train.data, res)

    
    est.mean <- EstimateMean(a.hat, num.test.days)
    est.var <- EstimateVar(sigma.hat, num.test.days)
    
    all.est.mean[i] <- est.mean
    all.est.var[i, i] <- est.var
    all.rsquare[i] <- rsquare
  }

  opt.weight <- Optimize(all.est.mean, all.est.var, gamma)

  all.test.init.price <- rep(0, num.assets)
  all.test.final.price <- rep(0, num.assets)
  for (i in 1:num.assets) {
    data <- all.data[, i]
    all.test.init.price[i] <- data[num.train.days]
    all.test.final.price[i] <- data[num.train.days + num.test.days]
  }

  opt.final.return <- GetFinalReturn(all.test.init.price, all.test.final.price, opt.weight)
  return(list(opt.final.return, all.rsquare))
}


RunAllWindows <- function(all.data, gamma) {
  num.train.days = 100
  num.test.days = 20
  num.one.window.days = num.train.days + num.test.days

  all.final.return <- c()
  all.rsquare <- c()
  start <- 1
  num.windows <- 1
  while (start + num.one.window.days <= dim(all.data)[1]) {
    all.data.slice <- SliceData(all.data, start, num.one.window.days)
    
    run.result <- RunOneWindow(all.data.slice, gamma)
    final.return <- run.result[[1]]
    rsquares <- run.result[[2]]  # a vector
    
    all.final.return[num.windows] <- final.return
    all.rsquare <- rbind(all.rsquare, rsquares)

    #cat(sprintf("%d~%d, %.4f\n", start, start+119, final.return))

    start <- start + 20
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


TestOneWindow <- function(all.data) {
  all.data.slice <- SliceData(all.data, 1, 120)

  for (i in 1:dim(all.data.slice)[2]) {
    data.slice <- all.data.slice[, i]

    fit.result <- FitModel(data.slice)
    a.hat <- fit.result[[1]]
    epsilon.hat <- fit.result[[2]]

    real.mean <- data.slice[120] - data.slice[100]
    esti.mean <- EstimateMean(a.hat, 20)

    plot(data.slice, type='l')
    cat(sprintf("a.hat = %f\n", a.hat))
    cat(sprintf("epsilon.hat = %f\n", epsilon.hat))
    cat(sprintf("real.mean = %f\n", real.mean))
    cat(sprintf("esti.mean = %f\n", esti.mean))
    cat("\n")
  }
}


# ===== Do experiment =====

load("data_for_exp_real.RData")

all.data <- data.frame(app.price, intel.price, msft.price,
                       lb.price, pvh.price, tjx.price,
                       cat.price, eaton.price, north.price)


#all.gamma <- c(0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.01, 0.1, 1, 5, 10)
# all.gamma <- c(0.01, 0.03, 0.05, 0.07, 0.1, 0.2, 0.3, 0.4, 0.5)
# all.results <- RunAllGammas(all.data, all.gamma)

results <- RunAllWindows(all.data, 0.1)
all.final.return <- results[[1]]
all.rsqure <- results[[2]]

par(mfrow = c(3,3))
asset.names <- c('APP(IT)', 'INTEL(IT)', 'MSFT(IT)',
                 'LB(CD)', 'PVH(CD)', 'TJX(CD)',
                 'CAT(IND)', 'EATON(IND)', 'NORTH(IND)')
for( i in 1:9){
  hist(all.rsqure[,i], main=asset.names[i], xlab = 'r-squared')
}
                      
                    