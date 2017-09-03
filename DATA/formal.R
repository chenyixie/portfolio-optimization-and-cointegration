library(tsDyn)
library(mnormt)
library(urca)
library(vars)
library(tseries)
library(expm)
library(Rdonlp2)
library(PerformanceAnalytics)

#Generate simple VECM 
set.seed(1234)
rm(list=ls())
simvecm = function(x0, n, B=matrix(0,nrow=n,ncol=n)){
  x = matrix(rep(0, n*L), nrow = n, ncol = L)
  x[,1]=x0
  r = matrix(rnorm(n*n,0,5),ncol=n,nrow=n)
  for(i in 2:L){
    #x[,i] = B %*% x[,i-1] +x[,i-1] + matrix(rnorm(n,mean=0,sd=v), nrow = n, ncol = 1)
    x[,i] = B %*% x[,i-1] +x[,i-1] + matrix(rmnorm(n=n,varcov=r %*% t(r)),nrow=n,ncol=1)#mvrnorm(n = 3, mu = 0, Sigma = 0.5)
  }
  return(x)
}


#Generate synthetic data
n = 3
L = 1000
N = 4
x0 = matrix(rep(rnorm(n,0,1),N),nrow = n,ncol = N)
simAns = rep(0,n*L*N)
dim(simAns) = c(n,L,N)
for (i in 1:N){
  B = matrix(c(-0.8,1.05,0.9),nrow=n) %*% matrix(c(0.7,-0.4,-1.02),ncol=n)
  simAns[,,i] = simvecm(x0[,i], n, B)
}

# plot graph for one group 
plot(simAns[1,,3],type = 'l', main = "Cointegrating Assets", ylab = "Price", xlab = "time")
lines(simAns[2,,3],col = 'red')
lines(simAns[3,,3],col = 'blue')

#plot graph for all data
for(i in N)
{
  Sys.sleep(1)
  plot(simAns[1,,i],type = 'l')
  lines(simAns[2,,i],col = 'red')
  lines(simAns[3,,i],col = 'blue')
}

# fit VECM model
beta = matrix(rep(0,n*N),nrow=n)
alpha = matrix(rep(0,n*N),nrow=n)
R = list()
for (i in 1:N){
  vecm.fit<- VECM(t(simAns[,,i]), r=1, lag=0, estim = "2OLS",include="none")
  s=summary(vecm.fit)
  R[[i]] = cov(s$residuals)
  beta[,i] = s$model.specific$beta
  alpha[,i] = s$coefficients
}
beta.Mean = rowMeans(beta)
alpha.Mean = rowMeans(alpha)
Bhat<-alpha.Mean%*%t(beta.Mean)  
Bhat
B

# Expectation for groups
Max_FN=50
value=rep(0,Max_FN)
FN = 10
# portexp = rep(0, Max_FN)
# portvar = rep(0, Max_FN)
identity = diag(n)
#FN = 10 #t+10
extra = list()
E = matrix(rep(0, n * N), nrow = n, ncol = N)
Exp = c()
SingleE = matrix(rep(0,Max_FN*12), ncol = 12)
for (FN in 1:Max_FN){
for (i in 1:N){
  E[, i] = (identity + alpha[, i] %*% t(beta[, i])) %^% (FN) %*% simAns[,L,i]
  extra[[i]] = as.vector(E[,i])
  }
E
SingleE[FN,] = c(extra[[1]], extra[[2]], extra[[3]], extra[[4]])
SingleE #1*12
}
plot(SingleE[,1], type = 'l')
plot(SingleE[,2], type = 'l')

# Expectation for portfolio
exp_port <- function(x){
  fun = 0
  for (i in 1: (n*N)){
    fun = fun + x[i]*SingleE[FN,i]
  }
  return(fun)
}
exp_port(rep(1/12,12))
# Variance for Groups
f = list()
f1 = list()
for (FN in 1: Max_FN){
for (j in 1:N) {
  f[[j]] = 
  f1[[j]] = matrix(0, nrow = n, ncol = n)
  for (i in 0:(FN - 1)) {
    f[[j]] = f[[j]] + t((identity + alpha[, j] %*% t(beta[, j])) %^% (i)) %*% R[[j]] %*% ((identity +
                                                                                             alpha[, j] %*% t(beta[, j])) %^% (i))
  }
}
}
f = list()
for (j in 1:N) {
  f[[j]] = list()
  for (FN in 1:Max_FN) {
    for (i in 0:(FN - 1)) {
      f[[j]][[FN]] = f[[j]][[FN]] + t((identity + alpha[, j] %*% t(beta[, j])) %^% (i)) %*% R[[j]] %*% ((identity +
                                                                                                           alpha[, j] %*% t(beta[, j])) %^% (i))
    }
  }
}

# Variance for portfolio 
varf = matrix(rep(0),ncol = (n*N), nrow = (n*N))

for (FN in 1:Max_FN){

varf[1:n,1:n] = f[[1]]
varf[(n+1):(2*n),(n+1):(2*n)] = f[[2]]
varf[(2*n+1):(3*n),(2*n+1):(3*n)] = f[[3]]
varf[(3*n+1):(4*n),(3*n+1):(4*n)] = f[[4]]
}



#######################
all_varf = list()
for (FN in 1:Max_FN) {
  f = list()
  for (j in 1:N) {
    f[[j]] = matrix(0, nrow = n, ncol = n)
    for (i in 0:(FN - 1)) {
      f[[j]] = f[[j]] + t((identity + alpha[, j] %*% t(beta[, j])) %^% (i)) %*% R[[j]] %*% ((identity +
                                                                                               alpha[, j] %*% t(beta[, j])) %^% (i))
    }
  }
  
  varf = matrix(rep(0),ncol = (n*N), nrow = (n*N))
  
  varf[1:n,1:n] = f[[1]]
  varf[(n+1):(2*n),(n+1):(2*n)] = f[[2]]
  varf[(2*n+1):(3*n),(2*n+1):(3*n)] = f[[3]]
  varf[(3*n+1):(4*n),(3*n+1):(4*n)] = f[[4]]
  
  all_varf[[FN]] = varf
}

var_port <- function(x){
  var = rep(0,50) 
  for (FN in 1:Max_FN){
  var[FN] = t(as.matrix(x)) %*% all_varf[[FN]] %*% as.matrix(x)
  }
  return(var)
}

plot(var_port(rep(1/12,12)), type = 'l')

# Portfolio Optimization
for (FN in 1:Max_FN){
final_opt<- function(gam){
  p = rep(1/12,12)    #initial
  par.l = rep(-0.3,12); par.u = rep(0.3,12)  ## lb and ub
  # objective function
  fn = function(x){
    f<- -exp_port(x)+gam*var_port(x)
    return(f)
  }
  
  # constraints sum of weight = 1
  A = matrix(rep(1,12),ncol =12)
  lin.l = 1; lin.u = 1
  
  # result
  ret = donlp2(p, fn, par.u=par.u, par.l=par.l, A, lin.l=lin.l, lin.u=lin.u)
  ret$par
  
  return(ret$par)
}
}
# value[FN] = t(as.matrix(final_opt(1))) %*% sim.price[1000,]
# portexp[FN] = exp_port(final_opt(1)) 
# portvar[FN] = var_port(final_opt(1))

value
plot(portexp, type = 'l')
plot(value, type = 'l')

#initial wealth 
sim.price = cbind(t(simAns[,,1]),t(simAns[,,2]),t(simAns[,,3]),t(simAns[,,4]))
iw = t(as.matrix(final_opt(1))) %*% sim.price[1000,]






