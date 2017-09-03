library(zoo)
library(R.matlab)
library(tsDyn)
library(vars)
setwd("~/Desktop/summer project/john k/DATA")

rm(list=ls())

spx <- readMat("spx.mat")
spximport <- spx$pp

#aquire data
tickers=spximport[1,1,1]$tickers
names=spximport[2,1,1]$names
sector=spximport[3,1,1]$sector
dt=as.Date(spximport[4,1,1]$dt,origin = "0000-01-01")
Index=spximport[5,1,1]$Index
#weights
w=spximport[6,1,1]$w
px=spximport[7,1,1]$px
#closing prices
px.clean=spximport[8,1,1]$px.clean

# log prices for IT
app.price = na.omit(px[,283])
app.logprice = log(app.price)
intel.price = na.omit(px[,681])
intel.logprice = log(intel.price)
msft.price = na.omit(px[,825])
msft.logprice = log(msft.price)

# graph for IT 
plot(app.price, type = 'l', main = "prices for IT Group", ylab ="prices", xlab = "time")
lines(intel.price, col = 2)
lines(msft.price, col= 4)
legend("bottomright", c("APPLIED MATERIALS INC", "INTEL CORP", "MICROSOFT CORP"),lty = c(1,1,1),col=c(1,2,4), cex = 0.7, bty = "n")

#cointegration test for IT
s1= cbind(app.price,intel.price,msft.price)
varest <- VAR(s1, p=1, type="const",lag.max=12*(6102/100)^0.25, ic="AIC")
lagLength <- max(12*(6102/100)^0.25,varest$p)
cointest <- ca.jo(s1, K= lagLength, type = "trace", ecdet = "const", spec = "transitory")
cointest@teststat
cointest@cval


# log price for Consumer Discretionary
lb.price = na.omit(px[,741])
pvh.price = na.omit(px[,965])
tjx.price = na.omit(px[,1100])

# Cointegration test for Consumer Discretionary
s2= cbind(lb.price,pvh.price,tjx.price)
varest2 <- VAR(s2, p=1, type="const",lag.max=12*(6102/100)^0.25, ic="AIC")
lagLength2 <- max(12*(6102/100)^0.25,varest2$p)
cointest2 <- ca.jo(s2, K= lagLength2, type = "trace", ecdet = "trend", spec = "transitory")
cointest2@teststat
cointest2@cval


# graph for Consumer Discretionary
plot(lb.price, type = 'l', main = "prices for Consumer Discretionary Group", ylab ="prices", xlab = "time")
lines(pvh.price, col = 2)
lines(tjx.price, col= 4)
legend("bottomright", c("APACHE CORP", "HESS CORP", "HELMERICH & PAYNE"),lty = c(1,1,1),col=c(1,2,4), cex = 0.7, bty = "n") 


# log price for Industrials
cat.price = na.omit(px[,384])
cat.logprice = log(cat.price)
eaton.price = na.omit(px[,538])
eaton.logprice = log(eaton.price)
north.price = na.omit(px[,872])
north.logprice = log(north.price)

# cointegration test fot industrials 
s3= cbind(cat.price,eaton.price, north.price)
varest3 <- VAR(s3, p=1, type="const",lag.max=12*(6000/100)^0.25, ic="AIC")
lagLength3 <- max(12*(6000/100)^0.25,varest3$p)
cointest3 <- ca.jo(s3, K= lagLength3, type = "trace", ecdet = "trend", spec = "longrun")
cointest3@teststat
cointest3@cval

industrials.vecm = ca.jo(s3)
industrials.vecm@teststat
industrials.vecm@cval

# graph for Industrials
plot(cat.price, type = 'l', main = "prices for Industrials Group", ylab ="prices", xlab = "time")
lines(eaton.price, col = 2)
lines(north.price, col= 4)
legend("bottomright", c("CATERPILLAR INC", "EATON CORP PLC", "NORTHROP GRUMMAN CORP"),lty = c(1,1,1),col=c(1,2,4), cex = 0.7, bty = "n") 

# Construct data frame for three groups 
IT.price = data.frame(app.price, intel.price, msft.price)
Energy.price = data.frame(apa.price, hess.price, hp.price)
Industrials.price = data.frame(cat.price, eaton.price, north.price)

# Fit IT data into VECM (r = 2, lag = 0)
real.train.days = 5102
real.test.days = 1000
IT.trainData = IT.price[1:real.train.days, ]
IT.testData  = IT.price[(real.train.days+1):6102, ]

IT.vecm.fit = VECM(IT.trainData, r=2, lag=2, estim="ML", include="none")
IT.s = summary(IT.vecm.fit)
IT.R = cov(IT.s$residuals)
IT.alpha = IT.s$coefficients
IT.beta = IT.s$model.specific$beta
IT.Pi = IT.alpha %*% t(IT.beta)

#residual plot for IT
IT.residual = resid(IT.vecm.fit)
plot(IT.residual[,1])
plot(IT.residual[,2])
plot(IT.residual[,3])
abline(0,0, col = 2)

## Fit Energy data into VECM (r = 1, lag = 0)
Energy.trainData = Energy.price[1:real.train.days, ]
Energy.testData  = Energy.price[(real.train.days+1):6102, ]

Energy.vecm.fit = VECM(Energy.trainData, r=1, lag=0, estim="ML", include="none")
Energy.s = summary(Energy.vecm.fit)
Energy.R = cov(Energy.s$residuals)
Energy.alpha = Energy.s$coefficients
Energy.beta = Energy.s$model.specific$beta
Energy.Pi = Energy.alpha %*% t(Energy.beta)

#residual plot for IT
Energy.residual = resid(Energy.vecm.fit)
plot(Energy.residual[,1])
plot(Energy.residual[,2])
plot(Energy.residual[,3])
abline(0,0, col = 2)


## Fit Industrials data into VECM (r = 1, lag = 0)
Industrials.trainData = Industrials.price[1:real.train.days, ]
Industrials.testData  = Industrials.price[(real.train.days+1):6102, ]

Industrials.vecm.fit = VECM(Industrials.trainData, r=1, lag=0, estim="ML", include="none")
Industrials.s = summary(Industrials.vecm.fit)
Industrials.R = cov(Industrials.s$residuals)
Industrials.alpha = Industrials.s$coefficients
Industrials.beta = Industrials.s$model.specific$beta
Industrials.Pi = Industrials.alpha %*% t(Industrials.beta)

#residual plot for IT
Industrials.residual = resid(Industrials.vecm.fit)
plot(Industrials.residual[,1])
plot(Industrials.residual[,2])
plot(Industrials.residual[,3])
abline(0,0, col = 2)


########################### test 
u1<-rnorm(500)
u2<-arima.sim(list(ar=0.6),n=500) #????????????????????????????????????
u3<-arima.sim(list(ar=.4),n=500)
y1<-cumsum(u1) #???????????????????????? y1
y2<-0.4*y1+u2
y3<-0.8*y1+u3
plot(y1, type = 'l')
lines(y2)
lines(y3)
data<-data.frame(y1=y1,y2=y2,y3=y3)
model.vecm<-ca.jo(data)
varest4 <- VAR(data, p=1, type="const",lag.max=12*(500/100)^0.25, ic="AIC")
lagLength4 <- max(12*(500/100)^0.25,varest4$p)
cointest4 <- ca.jo(data, K= 2, type = "trace", ecdet = "trend", spec = "longrun")
cointest4@teststat
cointest4@cval

