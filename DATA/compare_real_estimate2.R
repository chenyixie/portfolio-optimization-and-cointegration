library(zoo)
library(tsDyn)
library(mnormt)
library(urca)
library(vars)
library(tseries)
library(expm)
library(Rdonlp2)
library(PerformanceAnalytics)

# clear environment and set random seed
rm(list = ls())
set.seed(1234)

# parameters
n = 3     # number of assets in each group
N = 3     # number of groups
train_days = 200   # number of train data
test_days  = 20     # number of test data
pre_days   = 200    # number previous unused data
L1 = train_days + test_days + pre_days    # number of days, 200 + 200 + 20
L2 = train_days + test_days

#============================================================================
# 1. Generate synthetic data
#============================================================================

# initial price of VECM
initial_price = rnorm(n, 10, 3)
# initial_price = initial_price + abs(min(initial_price)) + 5

# create B for each group
# B is the product of two matrix with dimension (n * r) and rank r (here r=1)
B = list()
B[[1]] = matrix(c(-0.9, -0.6, 0.7), nrow = n) %*% matrix(c(0.65, 0.7, -0.9), ncol = n)
B[[2]] = matrix(c(-0.8, -0.7, 0.9), nrow = n) %*% matrix(c(0.7, 0.4, -0.03), ncol = n)
B[[3]] = matrix(c(-0.7,-0.65, 0.5), nrow = n) %*% matrix(c(0.6, 0.5, -0.7), ncol = n)

generate_data = function(epsilon) {
  
  # function that simulate VECM data
  simvecm = function(x0, B){
    # Parameters:
    #   x0: initial price
    #   B:  n*n matrix
    #
    # Returns:
    #   x: n * L matrix (n assets, L days)
    #
    #      Let VECM has no intercept and lag = 0, the equation is:
    #
    #          x[t] = x[t-1] + B * x[t-1] + epsilon
    #
    #      where B can be written as the product of two matrix with
    #      dimension (n ???? r) and rank r, and epsilon ~ N(0, R). (what's R?)
    
    x = matrix(rep(0, n*L1), nrow = n, ncol = L1)
    
    # set first column as initial price (price of first day)
    x[,1] = x0
    
    for(i in 2:L1){
      # epsilon: n*1 vector
      epsilon[, i] = matrix(rmnorm(n=1, varcov=r %*% t(r)/n), nrow=n, ncol=1)
      x[, i] = B %*% x[, i-1] + x[, i-1] + epsilon[, i]
    }
    
    return(x)
  }
  
  
  # x0[, i]: initial price of i_th group
  x0 = matrix(rep(initial_price, N), nrow = n, ncol = N)
  
  # simAns[, , i]: simulated data of i_th group
  simAns = rep(0, n * L1 * N)
  dim(simAns) = c(n, L1, N)
  
  #a = runif(n,-1,1);  b = runif(n,-1,1);c = runif(n,-1,1); d = runif(n,-1,1); e = runif(n,-1,1); f = runif(n,-1,1)
  #B[[1]] = matrix(a, nrow = n) %*% matrix(b, ncol = n)
  #B[[2]] = matrix(c, nrow = n) %*% matrix(d, ncol = n)
  #B[[3]] = matrix(e, nrow = n) %*% matrix(f, ncol = n)
  
  # generate data for each group
  for (i in 1:N) {
    # r: a random n*n matrix to compute the varcov of epsilon
    simAns[, , i] = simvecm(x0[, i], B[[i]])
  }
  
  # price must be larger than zero
  # simAns = simAns + abs(min(simAns)) + 5
  
  # remove previous unused days
  simAns = simAns[, -(1:pre_days), ]
  
  return(simAns)
  
  
  # plot graph for group 2
  plot(log(simAns[2, ,1]), type = 'l', main = "Cointegrating Assets", ylab = "Price", xlab = "days")
  lines(log(simAns[1, ,1]), col = 'red')
  lines(log(simAns[3, ,1]), col = 'blue')
}

#num_experiment = 1
#r = list()
#simAns = list()
#for (k in 1:num_experiment) {
#  r[[k]] = matrix(rnorm(n*n, 0, 0.5), ncol=n, nrow=n)
#  simAns[[k]] = generate_data(r[[k]])
#}

#simAns = simAns[[1]]
#r = r[[1]]

experiment <- function (gam, simAns) {
  
  #============================================================================
  # 2. fit VECM model
  #============================================================================
  
  trainData = simAns[, 1:train_days, ]
  testData  = simAns[, (train_days+1):L2, ]
  
  alpha    = matrix(rep(0, n*N), nrow=n)
  beta     = matrix(rep(0, n*N), nrow=n)
  vecm.fit = list()
  R        = list()
  Pi = list()
  
  for (i in 1:N){
    vecm.fit[[i]] = VECM(t(trainData[, , i]), r=1, lag=0, estim="ML", include="none")
    s             = summary(vecm.fit[[i]])
    R[[i]]        = cov(s$residuals)
    alpha[, i]    = s$coefficients
    beta[, i]     = s$model.specific$beta
    Pi[[i]] = alpha[, i] %*% t(beta[, i])
  }
  
  #============================================================================
  # 4. Expectation for groups
  #============================================================================
  
  identity = diag(n)
  extra = list()
  
  E = matrix(rep(0, n * N), nrow = n, ncol = N)
  #for (i in 1:N){
  #  for (k in 0:(test_days-1)){
  #    E[, i] = E[,i] + Pi[[i]] %*% (identity + Pi[[i]]) %^% (k) %*% simAns[,train_days,i]
  #  }
  #  extra[[i]] = as.vector(E[,i])
  # }  
  for (i in 1:N){
    E[, i] = (identity + Pi[[i]]) %^% (test_days) %*% simAns[,train_days,i] - simAns[,train_days,i]
    extra[[i]] = as.vector(E[,i])
  }
  
  # delta(X_t) = pi X_(t-1) + e_t
  # delta(X_(t+1)) = pi X_t + e_(t+1)
  # delta(X_(t+1))- delta(X_t) = pi delta(X_t) + e_(t+1) - e_t
  # delta(X_(t+1)) = (I + pi) delta(X_t)
  
  # pi \sum (I + pi)^i X_0 = (I + pi)^t X_t - X_0
  # E[X_t] = (I + pi)^t X_0
  SingleE = c(extra[[1]], extra[[2]], extra[[3]])
  SingleE2 = c(mean(testData[1,,1]), mean(testData[2,,1]), mean(testData[3,,1]),
               mean(testData[1,,2]), mean(testData[2,,2]), mean(testData[3,,2]),
               mean(testData[1,,3]), mean(testData[2,,3]), mean(testData[3,,3]))
  
  
  SingleE3 = c(testData[1,20,1]-trainData[1,200,1],testData[2,20,1]-trainData[2,200,1],testData[3,20,1]-trainData[3,200,1],
               testData[1,20,2]-trainData[1,200,2],testData[2,20,2]-trainData[2,200,2],testData[3,20,2]-trainData[3,200,2],
               testData[1,20,3]-trainData[1,200,3],testData[2,20,3]-trainData[2,200,3],testData[3,20,3]-trainData[3,200,3]
  )
  
  #============================================================================
  # 5. Expectation for portfolio
  #============================================================================
  
  # function to calculate portfolio expectation
  exp_port <- function(x, SingleE){
    # Parameters:
    #   x: weight vector
    #   SingleE: expectation of each asset
    # Returns:
    #   fun: portfolio expectation
    fun = 0
    
    for (i in 1:(n * N)){
      fun = fun + x[i] * SingleE[i]
    }
    
    return(fun)
  }
  
  
  #============================================================================
  # 6. Variance for portfolio
  #============================================================================
  
  # function to calculate portfolio variance
  var_port <- function(x, varf){
    # Parameters:
    #   x: weight vector
    #   varf: variance of what?
    # Returns:
    #   var: variance of portfolio
    var = t(as.matrix(x)) %*% varf %*% as.matrix(x)
    return(var)
  }
  
  
  f = list()
  identity = diag(n)
  for (g in 1:N) {
    f[[g]] = matrix(0, nrow = n, ncol = n)
    # f[[g]] = f[[g]] + R[[g]]
    for (k in 0:(test_days - 1)) {
      f[[g]] = f[[g]] + t((identity + Pi[[g]]) %^% (k)) %*% R[[g]] %*% ((identity + Pi[[g]]) %^% (k))
    }
  }
  
  varf = matrix(rep(0),ncol = (n*N), nrow = (n*N))
  
  varf[1:n,1:n] = f[[1]]
  varf[(n+1):(2*n),(n+1):(2*n)] = f[[2]]
  varf[(2*n+1):(3*n),(2*n+1):(3*n)] = f[[3]]
  
  # real variance for testdata var(X_220 - X_200) = var(X_220) = ?[3*3]
  # f2: var(X_201, X_202, ..., X_220)
  f2 = list()
  f2[[1]] = cov(t(testData[,,1]))
  f2[[2]] = cov(t(testData[,,2]))
  f2[[3]] = cov(t(testData[,,3]))
  
  varf2 = matrix(rep(0),ncol = (n*N), nrow = (n*N))
  varf2[1:n,1:n] = f2[[1]]
  varf2[(n+1):(2*n),(n+1):(2*n)] = f2[[2]]
  varf2[(2*n+1):(3*n),(2*n+1):(3*n)] = f2[[3]]
  
  
  
  #============================================================================
  # 7. Portfolio optimization
  #============================================================================
  
  ## 7.1 objective function
  
  fn = function(x){
    # Parameters:
    #   x: weight vector
    # Returns:
    #   f: target value
    f = -exp_port(x, SingleE3) + gam * var_port(x, varf2)
    return(f)
  }
  
  ## 7.2 low bound and upper bound
  
  
  num_assets = n * N  # number of assets
  init_cash = 10000
  
  # initial price of all assets
  new.price = cbind(t(simAns[,train_days,1]),t(simAns[,train_days,2]),t(simAns[,train_days,3]))
  
  # lower bound and upper bound
  par.l = rep(0, num_assets)
  par.u = c(rep(0.3 * init_cash, num_assets) / abs(new.price))
  
  ## 7.3 constraints sum of weight = 1: init_price * init_amount = init_cash
  A = matrix(new.price, 1, byrow = TRUE)
  lin.l = c(init_cash);
  lin.u = c(init_cash)
  
  p0 =  rep(init_cash, num_assets)/(num_assets * new.price)
  
  ret = donlp2(p0, fn, par.u=par.u, par.l=par.l, A, lin.l=lin.l, lin.u=lin.u)
  p = ret$par
  
  print_p <- function(p) {
    cat('               p: ')
    for (i in 1:num_assets) {
      cat(sprintf('%.2f ', p[i]))
    }
    cat('\n')
  }
  
  print_p(p)
  
  #============================================================================
  # 8. Final return of portfolio in test days
  #============================================================================
  
  # price of last test day
  Tprice = cbind(t(simAns[, L2, 1]),t(simAns[, L2, 2]),t(simAns[, L2, 3]))
  Final_Return = matrix(c(Tprice-new.price), ncol=n*N) %*% matrix(c(p), nrow=n*N)
  
  # get_final_return = function(Tprice, p, init_cash) {
  #return_of_asset = abs(c(Tprice) / c(new.price) - 1)
  #final_return = sum(c(p) * c(new.price) * (return_of_asset + 1)) / init_cash - 1
  #return(final_return)
  #  return(as.matrix(Tprice) %*% as.matrix(c(p)) / init_cash - 1)
  #  }
  
  # allocate initial cash to each asset and get the amount
  # Final_Return = get_final_return(Tprice, p, init_cash)
  #  average_return = as.matrix(Tprice) %*% as.matrix(c(p0)) / init_cash - 1
  
  #  cat(sprintf("optimized return: %6.3f\n", Final_Return))
  #  cat(sprintf("optimized object: %6.3f\n", -fn(p)))
  #  cat(sprintf("  average return: %6.3f\n", average_return[1, 1]))
  #  cat(sprintf("  average object: %6.3f\n", -fn(c(p0))))
  #  cat('\n')
  
  return(list(Final_Return, R, Pi, SingleE, SingleE3))
}


num_experiment = 10

# prepare experiment data: 1000 r for 1000 simAns
epsilon = list()
simAns = list()
# epsilon = matrix(rep(0,L1*n), nrow = 3)
r = matrix(rnorm(n*n, 0, 0.5), ncol=n, nrow=n)
for (k in 1:num_experiment) {
  epsilon[[k]] = matrix(rep(0,L1*n), nrow = 3)
  for(i in 1: L2){
  epsilon[[k]][,i] = matrix(rmnorm(n=1, varcov=r %*% t(r)/n), nrow=n, ncol=1)
  }
  simAns[[k]] = generate_data(epsilon[[k]])
}


# real R and real PI
real.R <- list()
real.Pi <- list()
for (e in 1:num_experiment) {
  real.R[[e]] <- r %*% t(r)/n
  real.Pi[[e]] <- B
}


# prepare experiment data: gamma
gam = c(0.0001, 0.0003, 0.0005, 0.001, 0.005, 0.01, 0.1, 1.0)
experiment_result = list()
d = list()
esti.R = list()
esti.Pi = list()
real.Exp = list()
esti.Exp = list()

for (i in 1:8) {
  experiment_result[[i]] = rep(0, num_experiment)
  
  cat(sprintf('gamma: %f\n\n', gam[i]))
  
  for (k in 1:num_experiment) {
    cat(sprintf('experiment %d\n\n', k))
    
    result_list = experiment(gam[i], simAns[[k]])
    experiment_result[[i]][k] = result_list[[1]]
    
    if (i == 1) {
      esti.R[[k]] = result_list[[2]]
      esti.Pi[[k]] = result_list[[3]]
      esti.Exp[[k]] = result_list[[4]]
      real.Exp[[k]] = result_list[[5]]
    }
  }
  
  d[[i]] = density(experiment_result[[i]]/10000)
  plot(d[[i]], type = 'l')
}


# plot all d
plot(d[[1]], type = 'l', ylim=c(0, 6),  xlim = c(-0.5,1), main = "Density of Final Return", xlab = "yield")
for (i in 2:8) {
  lines(d[[i]], col = i)
}
legend("bottomright", c( "0.0001", "0.0003", "0.0005", "0.001", "0.005", "0.01", "0.1", "1.0"),lty = c(1,1,1),col=c(1,2,3,4,5,6,7,8), cex = 0.7, bty = "n") 

# compare real PI and estimated PI (each group has one PI)
all_mse = matrix(0, nrow=num_experiment, ncol=N)   # each experiment has N mse (because of N PI)
for (k in 1:num_experiment) {
  real_pi_k = real.Pi[[k]]    # real pi of all groups of experiment k
  esti_pi_k = esti.Pi[[k]]    # estimated pi of all groups of experiment k
  
  for (g in 1:N) {
    real_pi_g = real_pi_k[[g]]    # real pi of group g of experiment k
    esti_pi_g = esti_pi_k[[g]]    # estimated pi of group g of experiment k
    # compare real_pi_g and esti_pi_g
    mse_pi_g = sum((real_pi_g - esti_pi_g) * (real_pi_g - esti_pi_g)) / (N*N)
    all_mse[k, g] = mse_pi_g
  }
}

# plot all mse
hist(c(all_mse), freq = FALSE , main = "MSE of Estimated Parameter Pi", xlab = "MSE", xlim = c(0,0.1))

# compare real R and estimated R (each group has one R)
all_mse_R = matrix(0, nrow=num_experiment, ncol=N)   # each experiment has N mse (because of N PI)
for (k in 1:num_experiment) {
  real_r_k = real.R[[k]]    # real pi of all groups of experiment k
  esti_r_k = esti.R[[k]]    # estimated pi of all groups of experiment k
  
  for (g in 1:N) {
    real_r_g = real_r_k[[g]]    # real pi of group g of experiment k
    esti_r_g = esti_r_k[[g]]    # estimated pi of group g of experiment k
    # compare real_pi_g and esti_pi_g
    mse_r_g = sum((real_r_g - esti_r_g) * (real_r_g - esti_r_g)) / (N*N)
    all_mse_R[k, g] = mse_r_g
  }
}

# plot all mse
hist(c(all_mse_R), freq = FALSE , main = "MSE of Estimated Residual", xlab = "MSE")


# compare SingleE
all.mse.exp = rep(0, num_experiment)
for (k in 1:num_experiment) {
  real.exp.k = real.Exp[[k]]   # vector SingleE3 of experiment k
  esti.exp.k = esti.Exp[[k]]   # vector SingleE of experiment k
  
  # compare
  mse.k = mean((real.exp.k - esti.exp.k) * (real.exp.k - esti.exp.k))
  all.mse.exp[k] = mse.k
}

hist(all.mse.exp, freq = FALSE , main = "MSE of Estimated Expectation", xlab = "MSE")

# true expectation + white nose
all.mse.wn = rep(0, num_experiment)
for (k in 1:num_experiment) {
  real.exp.k = real.Exp[[k]]   # vector SingleE3 of experiment k
  real.wn.exp.k = rnorm(n*N, 0, 1) + real.exp.k
  mse.k = mean((real.exp.k - real.wn.exp.k) * (real.exp.k - real.wn.exp.k))
  all.mse.wn[k] = mse.k
}

hist(all.mse.wn, freq = FALSE , main = "MSE of Real Expectation with Noise", xlab = "MSE", breaks = 10)

# Monte Carlo Mean
#for (k in 1:num_experiment){
#  real.exp.k = real.Exp[[k]] #vector SingleE3 of experiment k

# }

