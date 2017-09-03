ukx <- readMat("ukx.mat")
ukximport <- ukx$pp

#aquire data
Index=ukximport[1,1,1]$Index
tickers=ukximport[2,1,1]$tickers
names=ukximport[3,1,1]$names
sector=ukximport[4,1,1]$sector
#date
dt=as.Date(ukximport[5,1,1]$dt,origin = "0000-01-01")
#weights
w=ukximport[6,1,1]$w
px=ukximport[7,1,1]$px
#closing prices
px.clean=ukximport[8,1,1]$px.clean
#log prices
logp = log(px.clean)

logp[is.na(logp)]=0
logp[is.infinite(logp)]=0
logp=rbind(rep(0,ncol(logp)),logp)

s=na.omit(log(px.clean))

GetContinuousNonZeroData <- function(price.data) {
  # Get continuous none zero price data
  num.days = dim(price.data)
  for (t <- pr)
}

#Johansen test
for (i1 in 1:195){
  sector1 <- sector[[i1]][[1]][1, 1]
  if (sector1 == "NA") next
  
  for (i2 in 1:195){
    sector2 <- sector[[i2]][[1]][1, 1]
    if (sector2 == "NA") next
    
    for (i3 in 1:195){
      sector3 <- sector[[i3]][[1]][1, 1]
      if (sector3 == "NA") next
      
      if ((sector1 == sector2) & (sector2 == sector3)) {
      if ((i1 != i2) & (i1 != i3) & (i2 != i3)){
        ret1=ret[,i1]; ret2=ret[,i2]; ret3=ret[,i3]
        ret1<- ret1[ret1!=0]; ret2<- ret2[ret2!=0]; ret3<- ret3[ret3!=0]
        l=min(length(ret1),length(ret2),length(ret3))
        if(l>200){
          ret1=ret1[1:l]; ret2=ret2[1:l]; ret3 = ret3[1:l]
          data=data.frame(ret1,ret2,ret3)
          maxlag=12*(length(l)/100)^0.25
          varest=VAR(data,p=1,type="const",lag.max=maxlag,ic="AIC")
          lagLength=max(maxlag,varest$p)
          res=ca.jo(data,type="trace",ecdet="const",K=lagLength,spec="transitory")
          testStatistics=res@teststat
          criticalValues=res@cval
          if(testStatistics[length(testStatistics)] >= criticalValues[dim(criticalValues)[1],2]){
            print(c(i1,i2,i3))
            next
          }
          
        }
        
      }
    }
  }
}
}
