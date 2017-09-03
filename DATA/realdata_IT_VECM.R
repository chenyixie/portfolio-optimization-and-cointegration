library(zoo)
library(R.matlab)
library(tsDyn)
library(vars)
setwd("~/Desktop/summer project/john k/DATA")

spx <- readMat("spx.mat")
spximport <- spx$pp

#aquire data
tickers=spximport[1,1,1]$tickers
names=spximport[2,1,1]$names
sector=spximport[3,1,1]$sector
dt=as.Date(spximport[4,1,1]$dt,origin = "0000-01-01")
Index=spximport[5,1,1]$Index

#prices for IT 
app.price = na.omit(px[,283])
intel.price = na.omit(px[,681])
msft.price = na.omit(px[,825])


# Construct IT dataframe 
IT.data = data.frame(app.price[1:6000],intel.price[1:6000],msft.price[1:6000])

# graph for IT 
plot(app.price, type = 'l', main = "prices for IT Group", ylab ="prices", xlab = "time")
lines(intel.price, col = 2)
lines(msft.price, col= 4)
legend("bottomright", c("APPLIED MATERIALS INC", "INTEL CORP", "MICROSOFT CORP"),lty = c(1,1,1),col=c(1,2,4), cex = 0.7, bty = "n")

# ranks checking 
varest = list()
lagLength = c()
cointest = list()
teststat = list()
cval = list()
rank = c()
rank.ratio = c()
# daily freq = 1
# window1: 1~100, 6~105,...
# weekly freq =5
# window2: 5*(1~100), 5*(6~105),...
# freq = 20
# monthly window3: 20*(1~100), 20*(6~105),...
GetRank<- function(all.data,num.test.days){
for(i in 1: ((dim(all.data)[1] - 100)/5)){ 
  data <- all.data[(num.test.days*(i-1)+1):(100+num.test.days*(i-1)),]
  varest[[i]] <- VAR(data, p=1, type="const",lag.max=12*(100/100)^0.25, ic="AIC")
  lagLength[i] <- max(12*(100/100)^0.25,varest[[i]]$p)
  cointest[[i]] <- ca.jo(all.data[(num.test.days*(i-1)+1):(100+num.test.days*(i-1)),], K= lagLength[i], type = "trace", ecdet = "const", spec = "transitory")
  teststat[[i]] = cointest[[i]]@teststat
  cval[[i]] = cointest[[i]]@cval
  #if teststat >= cval, reject the H0
  # teststat[[i]][3] vs cval[[1]][6] if > continue  a: 3,2,1; b: 6,5,4
  # teststat[[i]][2] vs cval[[1]][5] if <= rank = 1
  # teststat[[i]][1] vs cval[[1]][4]
  
  rank[i] <- 3
  
  for (r in 0:2) {
    a <- 3 - r
    b <- 6 - r
    if (teststat[[i]][a] > cval[[i]][b]) {
      next
    }
    else {
      rank[i] <- r
      break
    }
  }
}
  return(rank)
}
# Daily frequency 
Daily.IT.rank <- GetRank(IT.data,5)
# Week frequency 
weekly.IT.data <- GetWeeklyData(all.data)[[1]]
weekly.IT.rank <- GetRank(weekly.IT.data,5)
#
monthly.IT.data <- GetMonthlyData(all.data)[[1]]
monthly.IT.rank <- GetRank(monthly.IT.data,5)

# the number of ranks w.r.t different frequencies
rank.summary = list(Daily.IT.rank,weekly.IT.rank,monthly.IT.rank)
num.zero = c()
num.ratio = c()
for(i in 1:3){
  num.zero[i] = sum(rank.summary[[i]] == 0)
  num.ratio[i] = 1-num.zero[i]/length(rank.summary[[i]])
}


# Consumer Discretionary
lb.price = na.omit(px[,741])
pvh.price = na.omit(px[,965])
tjx.price = na.omit(px[,1100])
# Construct Consumer Discretionary dataframe 
cd.data = data.frame(lb.price[1:6000],pvh.price[1:6000],tjx.price[1:6000])

# Daily frequency 
Daily.cd.rank <- GetRank(cd.data,5)
# Week frequency 
weekly.cd.data <- GetWeeklyData(all.data)[[2]]
weekly.cd.rank <- GetRank(weekly.cd.data,5)
#
monthly.cd.data <- GetMonthlyData(all.data)[[2]]
monthly.cd.rank <- GetRank(monthly.cd.data,5)

# the number of ranks w.r.t different frequencies
rank.summary.cd = list(Daily.cd.rank,weekly.cd.rank,monthly.cd.rank)
num.zero.cd = c()
num.ratio.cd = c()
for(i in 1:3){
  num.zero.cd[i] = sum(rank.summary.cd[[i]] == 0)
  num.ratio.cd[i] = 1-num.zero.cd[i]/length(rank.summary.cd[[i]])
}

# Industrials
cat.price = na.omit(px[,384])
eaton.price = na.omit(px[,538])
north.price = na.omit(px[,872])
# Construct Consumer Discretionary dataframe 
id.data = data.frame(cat.price[1:6000],eaton.price[1:6000],north.price[1:6000])

# Daily frequency 
Daily.id.rank <- GetRank(id.data,5)
# Week frequency 
weekly.id.data <- GetWeeklyData(all.data)[[3]]
weekly.id.rank <- GetRank(weekly.id.data,5)
#
monthly.id.data <- GetMonthlyData(all.data)[[3]]
monthly.id.rank <- GetRank(monthly.id.data,5)

# the number of ranks w.r.t different frequencies
rank.summary.id = list(Daily.id.rank,weekly.id.rank,monthly.id.rank)
num.zero.id = c()
num.ratio.id = c()
for(i in 1:3){
  num.zero.id[i] = sum(rank.summary.id[[i]] == 0)
  num.ratio.id[i] = 1-num.zero.id[i]/length(rank.summary.id[[i]])
}



