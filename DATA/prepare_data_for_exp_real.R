library(zoo)
library(R.matlab)

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
intel.price = na.omit(px[,681])
msft.price = na.omit(px[,825])

# log price for Consumer Discretionary
lb.price = na.omit(px[,741])
pvh.price = na.omit(px[,965])
tjx.price = na.omit(px[,1100])

# log price for Industrials
cat.price = na.omit(px[,384])
eaton.price = na.omit(px[,538])
north.price = na.omit(px[,872])

# cleaning data for 6000 days 

app.price = app.price[1:6000]
intel.price = intel.price[1:6000]
msft.price = msft.price[1:6000]

lb.price = lb.price[1:6000]
pvh.price = pvh.price[1:6000]
tjx.price = tjx.price[1:6000]

cat.price = cat.price[1:6000]
eaton.price = eaton.price[1:6000]
north.price = north.price[1:6000]

save.image("~/Desktop/summer project/john k/DATA/data_for_exp_real.RData")



