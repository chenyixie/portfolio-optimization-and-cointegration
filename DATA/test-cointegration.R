library(R.matlab)
library(zoo)
library(tseries)
library(vars)
library(urca)

ukx <- readMat("ukx.mat")
ukximport <- ukx$pp

#aquire data
Index=ukximport[1,1,1]$Index
tickers=ukximport[2,1,1]$tickers
names=ukximport[3,1,1]$names
sector=ukximport[4,1,1]$sector
#date
#dt=as.Date(ukximport[5,1,1]$dt,origin = "0000-01-01")
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

#calculate daily log-returns ( log-return = log(pt/p(t-1)) )
ret=log(px.clean[2:dim(px.clean)[1],])-log(px.clean[1:(dim(px.clean)[1]-1),])
ret[is.na(ret)]=0
ret[is.infinite(ret)]=0
ret=rbind(rep(0,ncol(ret)),ret)

s=na.omit(ret)

#days in each month#
days=diff(seq(as.Date(dt[1]), as.Date(dt[length(dt)]), by = "month"))
days=c(days,19)
w=rbind(w,rep(0,ncol(w)))
w.days=w[rep(1:nrow(w), times = days), ]

years=rep(0,12)
years[1]=sum(days[1:12])
for(i in 2:11){
  years[i]=sum(days[1:(i*12)])-sum(days[1:((i-1)*12)])
}
years[12]=sum(days[133:134])

# January weights replicated over entire year#
w.jan=matrix(nrow=12,ncol=dim(w)[2])
for(i in 1:12){
  w.jan[i,]=w[1+(i-1)*12,]
}
w.jan.days=w.jan[rep(1:nrow(w.jan), times = years), ]

#cumulative returns

index.ret2=vector("numeric") 
cumul.index.ret2=vector("numeric")
cumul.index.ret2[1,]=100
for(i in 1:nrow(ret)){
  index.ret2[i]=sum(ret[i,]*w.jan.days[i,])
  cumul.index.ret2[i+1]=cumul.index.ret2[i,]*(1+index.ret2[i,])
}

adf.test(cumul.index.ret2,k=0)$p.value


#spread functions
#find the index locations of N largest weights
indexN=function(N,mat){
  m=matrix(nrow=nrow(mat),ncol=N)
  for(i in 1:nrow(mat)){
    ndx <- order(mat[i,], decreasing = T)[1:N]
    m[i,]=ndx
  }
  return(m)
}

index.retN=function(N,weights,returns){
  index.retN=vector("numeric"); cumul.index.retN=vector("numeric")
  cumul.index.retN[1]=100
  index.matrix=indexN(N,weights)
  for(i in 1:nrow(weights)){
    index.retN[i]=sum(returns[i,][index.matrix[i,]]*weights[i,][index.matrix[i,]])
    cumul.index.retN[i+1]=cumul.index.retN[i]*(1+index.retN[i])
  }
  return(list(index.retN,cumul.index.retN))
}

# 30 stock in chosen index
index2.4=index.retN(30,w.jan.days,ret)
ret2.30=index2.30[[1]]
cumul2.30=index2.30[[2]]

#diagnostic check

sd(ret2.30-index.ret2)

par(mfrow=c(2,1))
ts.plot(index.ret2,main="FTSE 100 1.1.1990-19.3.2014",ylab="Daily log-returns of index")
lines(ret2.30,col="red",lty=3)
legend("bottomright",c("All stocks","Largest 30 stocks"),lty=c(1,3),col=c("black","red"),cex=0.5)


ts.plot(cumul.index.ret2,main="FTSE 100 1.1.1990-19.3.2014",ylab="Daily cumulative log-returns of index")
lines(cumul2.30,col="red",lty=3)
legend("bottomright",c("All stocks","Largest 30 stocks"),lty=c(1,3),col=c("black","red"),cex=0.5)
dev.off()

#check for cointegration

create.spread.opt<-function(s1,s2,cor_thresh,p_thresh){
  l=min(length(s1),length(s2))
  s1=s1[1:l]; s2=s2[1:l]
  s1=log(s1); s2=log(s2)
  
  if(adf.test(s1,k=0)$p.value<p_thresh || adf.test(s2,k=0)$p.value<p_thresh){
    return("S1 or S2 not I(1)")}
  if(abs(cor(diff(s1),diff(s2)))<cor_thresh) return("Corr of log-returns below threshold")
  if(abs(cor(diff(s1),diff(s3)))<cor_thresh) return("Corr of log-returns below threshold")
  
  lm1=lm(s1~s2)
  u1=lm1$coefficients[1]; u2=lm2$coefficients[1]
  g1=lm1$coefficients[2]; g2=lm2$coefficients[2]
  z1=lm1$residuals; z2=lm2$residuals
  
  lag1=floor(12*(length(z1)/100)^0.25)
  lag2=floor(12*(length(z2)/100)^0.25)
  
  while(abs(adf.test(z1,k=lag1)@test$statistic)<1.6 && lag1>0) lag1=lag1-1
  while(abs(adf.test(z2,k=lag2)@test$statistic)<1.6 && lag2>0) lag2=lag2-1
  
  list1=adf.test(z1,k=lag1)
  list2=adf.test(z2,k=lag2)
  
  statistic1=list1@test$statistic; statistic2=list2@test$statistic
  p.val1=list1@test$p.val; p.val2=list2@test$p.val
  
  if(statistic1<statistic2 && p.val1<=0.05 && kpss.test(z1)$p.value>=0.05){
    #cat("s1=g*s2+u\n")
    return(list(c(u1,g1,statistic1),"z1"=z1))
  }else if(statistic2<statistic1 && p.val2<=0.05 && kpss.test(z2)$p.value>=0.05){
    #cat("s2=g*s1+u\n")
    return(list(c(u2,g2,statistic2),"z2"=z2))
  }else{
    return("No Cointegration between S1 and S2")
  }
}

#check for all combinations
from2=4; to2=20; by2=4; days2=180; cor_thresh2=0.9; pval_thresh2=0.01
for(j in seq(from2,to2,by2)){
  cat("\nNumber of stocks=",j,"\n\n")
  #culmulative return 
  index.j2=index.retN(j,w.jan.days,ret)[[2]]
  for(i in 1:floor(length(cumul.index.ret2)/days2)){
    spread2=create.spread.opt(cumul.index.ret2[(days2*(i-1)):(days2*i)],index.j2[(days2*(i-1)):(days2*i)],cor_thresh2,pval_thresh2)[[1]]
    print(c(spread2,i))
  }
  rm(spread2,index.j2,i,j)
}

#using Johansen Test
maxlag2=12*(length(dt)/100)^0.25
for(j in seq(from2,to2,by2)){
  cat("\nNumber of stocks=",j,"\n\n")
  index.j2=index.retN(j,w.jan.days,ret)[[2]]
  for(i in 1:floor(length(cumul.index.ret2)/days2)){
    cumul.ret.all2=log(cumul.index.ret2[(days2*(i-1)):(days2*i)])
    cumul.ret.j2=log(index.j2[(days2*(i-1)):(days2*i)])
    data2=data.frame(cumul.ret.all2,cumul.ret.j2)
    
    varest2=VAR(data2,p=1,type="const",lag.max=maxlag2,ic="AIC")
    lagLength2=max(maxlag2,varest2$p)
    res2=ca.jo(data2,type="trace",ecdet="const",K=lagLength2,spec="transitory")
    testStatistics2=res2@teststat
    criticalValues2=res2@cval
    if(testStatistics2[length(testStatistics2)] >= criticalValues2[length(testStatistics2)+4]){
      print(c(res2@V[1:(ncol(data2)+1),which.max(res2@lambda)],i))
    }
  }
  rm(index.j2,cumul.ret.all2,cumul.ret.j2,data2,varest2,lagLength2,res2,testStatistics2,criticalValues2,i,j)
} 


#Johansen test
coint_set <- NULL
for (i1 in 1:195) {
  for (i2 in 1:195) {
    for (i3 in 1:195) {
      if (sector[i1] == sector[i2]) {
        if ((i1 != i2) &(i1 != i3) & (i2 != i3) & (!any((c(i1, i2, i3) %in% coint_set)))) {
          ret1 = ret[, i1]
          ret2 = ret[, i2]
          ret3 = ret[, i3]
          ret1 <-ret1[-ret1 != 0]
          ret2 <- ret2[-ret2 != 0]
          ret3 <- ret3[-ret3 != 0]
          l = min(length(ret1), length(ret2), length(ret3))
          if (l > 200) {
            ret1 = ret1[1:l]
            ret2 = ret2[1:l]
            ret3 = ret3[1:l]
            data = data.frame(ret1, ret2, ret3)
            varest = VAR(data, p = 1, type = "const", lag.max = maxlag, ic = "AIC" )
            maxlag = 12 * (l / 100) ^ 0.25
            lagLength = max(maxlag, varest$p)
            res = ca.jo(data, type = "trace", ecdet = "const", K = lagLength, spec = "transitory")
            testStatistics = res@teststat
            criticalValues = res@cval
            if (testStatistics[length(testStatistics)] >= criticalValues[dim(criticalValues)[1], 2]) {
              print(c(i1, i2, i3))
              coint_set <- unique(c(coint_set, i1, i2, i3))
            }
          }
        }
      }
    }
  }
}


length(ret[,1])
adf=rep(0,ncol(ret))
f<-function(x){
  ur.df(x, type = "none", selectlags = "AIC")
}

ret1=ret[,1]; ret2=ret[,2]; ret3=ret[,4]
ret1<- ret1[-ret1!=0]; ret2<- ret2[-ret2!=0]; ret3<- ret3[-ret3!=0]
l=min(length(ret1),length(ret2),length(ret3))
ret1=ret1[1:l]; ret2=ret2[1:l]; ret3 = ret3[1:l]
data=data.frame(ret1,ret2,ret3)
varest=VAR(data,p=1,type="const",lag.max=maxlag,ic="AIC")
lagLength=max(maxlag,varest$p)
res=ca.jo(data,type="trace",ecdet="const",K=2,spec="transitory")
testStatistics=res@teststat
criticalValues=res@cval
testStatistics
criticalValues

plot(ret1, type='l', ylim= c(-0.2,0.1),xlim = c(0,500))
lines(ret2, col= 'red' )
lines(ret3, col = 'blue')
plot(logp[,1],type='l',xlim = c(0,500))
lines(logp[,2], col='red')
lines(logp[,4], col='blue')

plot(ret1, type = 'l')
ret1<- ret[,1] 
ret1<- ret1[-ret1!=0]
ret1
length(ret1)



