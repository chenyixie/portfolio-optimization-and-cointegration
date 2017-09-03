#Portfolio optimization for synthetic data _ markowitz_t_dep

# clear environment and set random seed
rm(list = ls())
set.seed(4321)
setwd("~/Desktop/final project/codes")
source("shared_functions.R")
source("models/markowitz/t_dependent.R")

#Generate Data (1000 samples, 105 days)
mu <- rnorm(9, 0, 0.01)
sigma <- genPositiveDefMat(9)$Sigma
degree <- 2


num.test.days <- 5
num.experiment <- 1000
gamma = c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.1, 1.0, 10)

init.price <- abs(rnorm(9, 5, 5))
num.days <- 105
all.data = list()
train.data = list()
test.data = list()
for(i in 1:num.experiment){
  all.data[[i]] <- GenerateData(degree,mu, sigma, init.price, num.days)
  train.data[[i]] <- all.data[[i]][1:100,]
  test.data[[i]]<- all.data[[i]][101:105,]
}

# ===== Optimizaton for estimated parameters =====
#Optimize 
fit.result = list()
asset.mean = list()
asset.var = list()
opt.weight = list()

GetWeight <- function(gamma){
  fit.result = list()
  asset.mean = list()
  asset.var = list()
  opt.weight = list()
  for(i in 1:num.experiment){
    fit.result[[i]] <- FitModel(train.data[[i]])
    asset.mean[[i]]<- EstimateMean(fit.result[[i]], num.test.days)
    asset.var[[i]] <- EstimateVar(fit.result[[i]], num.test.days)
    opt.weight[[i]]<-Optimize(asset.mean[[i]], asset.var[[i]], gamma)
  }
  return(opt.weight)
}

# Final return
initial.price <- list()
final.price <- list()
final.return <- c()

GetReturn <- function(weight){
  for (i in 1:num.experiment){
    initial.price[[i]] <- train.data[[i]][100,]
    final.price[[i]] <- test.data[[i]][5,]
    final.return[i] <- GetFinalReturn(initial.price[[i]], final.price[[i]], weight[[i]])
  }
  return(final.return)
}

final.weight<- list()
all.final.return <- list()
final.density = list()
for(i in 1:8){
  final.weight[[i]] <- GetWeight(gamma[i])
  all.final.return[[i]]<- GetReturn(final.weight[[i]])
  final.density[[i]] = density(all.final.return[[i]])
}

# Plot density graphs
plot(final.density[[8]], type = 'l', col = '1', main = 'Estimated parameters', xlab = 'Expected Reward')
for (i in 1:7) {
  lines(final.density[[i]], col =i)
  legend("bottomright",legend = c("0.0001", "0.0005", "0.001", "0.005", "0.01", "0.1", "1.0", "10"), 
         lty=c(1,1,1,1,1,1,1,1),col=c('1','2','3','4','5','6','7','8'), bty = 'n', cex = 0.7)
}


# ===== OPtimization for true parameters =====
mu <- rnorm(9, 0, 0.01)
sigma <- genPositiveDefMat(9)$Sigma

#Optimize 

GetTrueWeight <- function(gamma){
  asset.mean<- mu * num.test.days
  asset.var <- sigma * num.test.days
  opt.true.weight<-Optimize(asset.mean, asset.var, gamma)
  return(opt.true.weight)
}

# Final return
initial.price <- list()
final.price <- list()
final.true.return <- c()

GetTrueReturn <- function(weight){
  for (i in 1:num.experiment){
    initial.price[[i]] <- train.data[[i]][100,]
    final.price[[i]] <- test.data[[i]][5,]
    final.true.return[i] <- GetFinalReturn(initial.price[[i]], final.price[[i]], weight)
  }
  return(final.true.return)
}

final.true.weight <- list()
all.final.true.return <- list()
final.true.density <- list()
for(i in 1:8){
  final.true.weight[[i]] <- GetTrueWeight(gamma[i])
  all.final.true.return[[i]]<- GetTrueReturn(final.true.weight[[i]])
  final.true.density[[i]] = density(all.final.true.return[[i]])
}

#Plot density graphs
plot(final.true.density[[8]], type = 'l', col = '1', main = 'True parameters', xlab = 'Expected Reward')
for (i in 1:7) {
  lines(final.true.density[[i]], col =i)
  legend("bottomright",legend = c("0.0001", "0.0005", "0.001", "0.005", "0.01", "0.1", "1.0", "10"), 
         lty=c(1,1,1,1,1,1,1,1),col=c('1','2','3','4','5','6','7','8'), bty = 'n', cex = 0.7)
}

