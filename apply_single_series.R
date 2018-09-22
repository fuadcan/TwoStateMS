setwd("TwoStateMS")
source("MSArfimaFit.R")
source("dlvPath.R")

# Read data
dat  <- read.csv("data/madisonFrom-1930.csv",sep=";")
dat  <- log(dat)
yearrange <- seq(1930,1930+nrow(dat)-1)

# Log differences of each pair
pdat <- apply(combn(ncol(dat),2),2,function(c) dat[,c[1]]-dat[,c[2]] )
pnames <- apply(combn(ncol(dat),2),2,function(c) paste0(colnames(dat)[c[1]],"-",colnames(dat)[c[2]] ))
colnames(pdat) <- pnames

# The first series of the data
ser     <- pdat[,1]
sername <- pnames[1]


  
## only d switch states
opt  <- MSArfima.fit(ser,"D")
path <- dlvPath_d(opt$par,ser)
dplot <- c(opt$par[1:2] %*% path)

plot(yearrange, dplot, type='l',xlab = "Year",ylab = "d")
title(sername)



## only d and mu switch states
opt  <- MSArfima.fit(ser,"DM")
path <- dlvPath_dm(opt$par,ser)
dplot  <- c(opt$par[1:2] %*% path)
muplot <- c(tail(opt$par,2) %*% path)

plot(yearrange, dplot, type='l',xlab = "Year",ylab = "d")
title(sername)

plot(yearrange, muplot, type='l',xlab = "Year",ylab = "mu")
title(sername)


## only d and sigma switch states
opt  <- MSArfima.fit(ser,"DS")
path <- dlvPath_ds(opt$par,ser)
dplot     <- c(opt$par[1:2] %*% path)
sigmaplot <- c(tail(opt$par,5:6) %*% path)

plot(yearrange, dplot, type='l',xlab = "Year",ylab = "d")
title(sername)

plot(yearrange, sigmaplot, type='l',xlab = "Year",ylab = "sigma")
title(sername)



## d, sigma and mu switch states
opt  <- MSArfima.fit(ser,"DMS")
path <- dlvPath_dsm(opt$par,ser)
dplot <- c(opt$par[1:2] %*% path)
muplot    <- c(tail(opt$par,7:8) %*% path)
sigmaplot <- c(tail(opt$par,5:6) %*% path)

plot(yearrange, dplot, type='l',xlab = "Year",ylab = "d")
title(sername)

plot(yearrange, muplot, type='l',xlab = "Year",ylab = "mu")
title(sername)

plot(yearrange, sigmaplot, type='l',xlab = "Year",ylab = "sigma")
title(sername)
