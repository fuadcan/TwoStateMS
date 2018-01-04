# source("lnviDM2.R")
source("mclapplyhack.R")
library("urca")

convDLV_d <- function(yearOrRegion){
  
  if(is.numeric(yearOrRegion)){year      <- yearOrRegion
  filename  <- paste0("data/madisonFrom-",year,".csv")
  z         <- read.table(filename,header = T,sep = ";")
  z         <- data.matrix(z)} else if(yearOrRegion=="Maddison"){
    filename  <- "data/madisonFromMaxNA2-1950.csv"
    z         <- read.table(filename,header = T,sep = " ")
    z         <- data.matrix(z)
  } else if(yearOrRegion == "Europe") {z <- data.matrix(read.table("data/mds_Europe-1950.csv",header = T,sep = ";"))} else {
    fname1 <- paste0("data/mds_G7-1950.csv"); fname2 <- paste0("data/mds_Europe-1950.csv"); fname3 <- paste0("data/mds_S&P-1950.csv")
  z_g7         <- read.table(fname1,header = T,sep = ";"); z_g7 <- data.matrix(z_g7)
  z_eu         <- read.table(fname2,header = T,sep = ";"); z_eu <- data.matrix(z_eu)
  z_sp         <- read.table(fname3,header = T,sep = ";"); z_sp <- data.matrix(z_sp)
  
  
  if(yearOrRegion=="Europe+G7"|yearOrRegion=="G7+Europe")
  {z<- matrix(c(z_g7,z_eu),dim(z_g7)[1],); colnames(z)<- c(colnames(z_g7),colnames(z_eu))} else
    if(yearOrRegion=="Europe+S&P"|yearOrRegion=="S&P+Europe")
    {z<- matrix(c(z_eu,z_sp),dim(z_eu)[1],); colnames(z)<- c(colnames(z_eu),colnames(z_sp))} else
      if(yearOrRegion=="G7+S&P"|yearOrRegion=="S&P+G7")
      {z<- matrix(c(z_g7,z_sp),dim(z_g7)[1],); colnames(z)<- c(colnames(z_g7),colnames(z_sp))} else {stop("Unknown Country List")}
  if(sum(duplicated(t(z)))!=0){z<- z[,-which(duplicated(t(z)))]}
  }
  
  
  
  z <- log(z)
  n <- ncol(z)
  cmbn<- combn(n,2)
  
  pPanel <- sapply(1:ncol(cmbn), function(i) z[,cmbn[1,i]] - z[,cmbn[2,i]])
  
  lowerV <- c(-2,-2,.8,.8,0.001,-50)
  upperV <- c(2,2,.999,.999,50,50)
  inits  <- c(.75,  1.2,  0.90,  0.99,  0.01, 0.1)
  
  nc    <- ncol(pPanel)
  step  <- floor(nc/5)
  from  <- 1 + (0:4)*step
  to    <- (1:5)*step
  to[5] <- nc

  
  res1 <- simplify2array(mclapply.hack(from[1]:to[1], function(n){
    cat(paste0(n,"\n"))
    temp <- tryCatch(optim(inits, function(p) -lnviD2(p,pPanel[,n]),lower = lowerV,
                           upper = upperV,method="L-BFGS-B"), error=function(e){NA})
    temp <- tryCatch(c(temp$par,temp$value),error=function(e){rep(NA,7)}); return(temp)}))
  
  cat(paste0("step 1 with ",to[1]-from[1]+1, " series is done\n"))
  save(res1, file=paste0("output/d_",yearOrRegion,"_D_res1.rda"))
  
  res2 <- simplify2array(mclapply.hack(from[2]:to[2], function(n){
    temp <- tryCatch(optim(inits, function(p) -lnviD2(p,pPanel[,n]),lower = lowerV,
                           upper = upperV,method="L-BFGS-B"), error=function(e){NA})
    temp <- tryCatch(c(temp$par,temp$value),error=function(e){rep(NA,7)}); return(temp)}))
  
  save(res2, file=paste0("output/d_",yearOrRegion,"_D_res2.rda"))
  cat(paste0("step 2 with ",to[1]-from[1]+1, " series is done\n"))
  
  res3 <- simplify2array(mclapply.hack(from[3]:to[3], function(n){
    temp <- tryCatch(optim(inits, function(p) -lnviD2(p,pPanel[,n]),lower = lowerV,
                           upper = upperV,method="L-BFGS-B"), error=function(e){NA})
    temp <- tryCatch(c(temp$par,temp$value),error=function(e){rep(NA,7)}); return(temp)}))
  
  save(res3, file=paste0("output/d_",yearOrRegion,"_D_res3.rda"))
  cat(paste0("step 3 with ",to[1]-from[1]+1, " series is done\n"))
  
  res4 <- simplify2array(mclapply.hack(from[4]:to[4], function(n){
    temp <- tryCatch(optim(inits, function(p) -lnviD2(p,pPanel[,n]),lower = lowerV,
                           upper = upperV,method="L-BFGS-B"), error=function(e){NA})
    temp <- tryCatch(c(temp$par,temp$value),error=function(e){rep(NA,7)}); return(temp)}))
  
  save(res4, file=paste0("output/d_",yearOrRegion,"_D_res4.rda"))
  cat(paste0("step 4 with ",to[1]-from[1]+1, " series is done\n"))
  
  res5 <- simplify2array(mclapply.hack(from[5]:to[5], function(n){
    temp <- tryCatch(optim(inits, function(p) -lnviD2(p,pPanel[,n]),lower = lowerV,
                           upper = upperV,method="L-BFGS-B"), error=function(e){NA})
    temp <- tryCatch(c(temp$par,temp$value),error=function(e){rep(NA,7)}); return(temp)}))
  
  save(res5, file=paste0("output/d_",yearOrRegion,"_D_res5.rda"))
  cat(paste0("step 5 with ",to[1]-from[1]+1, " series is done\n"))
  
  res   <- cbind(res1,res2,res3,res4,res5)
  
  resdf <- apply(pPanel,2, function(p){r <- ur.df(p,"drift"); hp <- r@cval[1,2] > r@teststat[1]; return(c(r@teststat[1],hp))})
  res   <- t(rbind(res,resdf))
  
  save(res, file=paste0("output/d_",yearOrRegion,"_D_resALL.rda"))
  
  return(res)
  
}
