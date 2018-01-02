source("lnviDM2.R")
source("lnviD2.R")
source("mclapplyhack.R")


convDLV <- function(yearOrRegion,type){
  # type = 1 for allowing only d switching states,
  # type = 2 for allowing both d and mu switching states
  # yearOrRegion <- 1930; type <- 1
  
  # Reading data
  if(is.numeric(yearOrRegion)){year      <- yearOrRegion
  filename  <- paste0("data/madisonFrom-",year,".csv")
  z         <- read.table(filename,header = T,sep = ";")
  z         <- data.matrix(z)} else if(yearOrRegion=="Maddison"){
    filename  <- "data/madisonFromMaxNA2-1950.csv"
    z         <- read.table(filename,header = F,sep = " ",skip = 1)
    # colnames(z) <- as.character(z[1,]); z <- z[-1,]
    # z         <- data.matrix(z)
  } else if(yearOrRegion=="Europe"){z <- data.matrix(read.table("data/mds_Europe-1950.csv",header = T,sep = ";"))} else {fname1 <- paste0("data/mds_G7-1950.csv"); fname2 <- paste0("data/mds_Europe-1950.csv"); fname3 <- paste0("data/mds_S&P-1950.csv")
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
  
  # Constructing Pairs Panel
  pPanel <- sapply(1:ncol(cmbn), function(i) z[,cmbn[1,i]] - z[,cmbn[2,i]])
  
  # Parameters
  list_dm <- list(lowerV = c(-2,-2,.8,.8,0.001,-5,-5),
                  upperV = c(2,2,.999,.999,5,5,5),
                  inits1 = c(1.2,  1.7,  0.9,  0.9,  0.01, 0.1,0.1),
                  inits2 = c(.8,  1.2,  0.9,  0.9,  0.01, 0.1,0.1),
                  inits3 = c(.5,  1,  0.9,  0.9,  0.01, 0.1,0.1))
  list_d  <- list(lowerV = c(-2,-2,.8,.8,0.001,-5),
                  upperV = c(2,2,.999,.999,5,5),
                  inits1 = c(1.2,  1.7,  0.9,  0.9,  0.01, 0.1),
                  inits2 = c(.8,  1.2,  0.9,  0.9,  0.01, 0.1),
                  inits3 = c(.5,  1,  0.9,  0.9,  0.01, 0.1))
  
  lst     <- if(type==1) list_d else list_dm
  DLValgo <- if(type==1) lnviD2 else lnviDM2
  outname <- if(type==1) "D" else "DM"
  
  # Constraints
  const_mat <- matrix(0,length(lst$lowerV),length(lst$lowerV))
  diag(const_mat) <- 1
  const_mat <- rbind(const_mat,-const_mat)
  const_mat <- cbind(const_mat,c(lst$lowerV,-lst$upperV))
  
  # Function for optimizing DLValgo
  optimizator <- function(i){
    inits1  <- lst$inits1
    temp1 <- constrOptim(inits1, function(p) -DLValgo(p,pPanel[,i]), NULL, ui = const_mat[,-ncol(const_mat)], const_mat[,ncol(const_mat)])
    if(temp1$convergence!=0){out <- temp1} else {
      inits2 <- lst$inits2
      temp2  <- constrOptim(inits2, function(p) -DLValgo(p,pPanel[,i]), NULL, ui = const_mat[,-ncol(const_mat)], const_mat[,ncol(const_mat)])
      if(temp2$convergence != 0) {out <- temp2} else {
        inits3 <- lst$inits3
        temp3  <- constrOptim(inits3, function(p) -DLValgo(p,pPanel[,i]), NULL, ui = const_mat[,-ncol(const_mat)], const_mat[,ncol(const_mat)])
        if(temp3$convergence != 0) {out <- temp3} else {
          templist <- list(temp1,temp2,temp3)
          lkls     <- sapply(templist, function(t) t$value)
          out      <- templist[[which(min(lkls)==lkls)[1]]]
        }
      }
    }
    return(out)
  }
  
  # The outputs are saved in 5 parts and a combined one to prevent any loss of computations
  
  nc    <- ncol(pPanel)
  step  <- floor(nc/5)
  from  <- 1 + (0:4)*step
  to    <- (1:5)*step
  to[5] <- nc
  
  
  res1 <- mclapply.hack(from[1]:to[1], function(n){
    cat(paste0(n,"\n"))
    temp <- tryCatch(optimizator(n), error=function(e){rep(NA,length(lst$inits1))})
    return(temp)})
  
  cat(paste0("step 1 with ",to[1]-from[1]+1, " series is done\n"))
  save(res1, file=paste0("output/d_",yearOrRegion,"_",outname,"_res1.rda"))
  
  res2 <- mclapply.hack(from[2]:to[2], function(n){
    cat(paste0(n,"\n"))
    temp <- tryCatch(optimizator(n), error=function(e){rep(NA,length(lst$inits1))})
    return(temp)})
  
  save(res2, file=paste0("output/d_",yearOrRegion,"_",outname,"_res2.rda"))
  cat(paste0("step 2 with ",to[1]-from[1]+1, " series is done\n"))
  
  res3 <- mclapply.hack(from[3]:to[3], function(n){
    cat(paste0(n,"\n"))
    temp <- tryCatch(optimizator(n), error=function(e){rep(NA,length(lst$inits1))})
    return(temp)})
  
  save(res3, file=paste0("output/d_",yearOrRegion,"_",outname,"_res3.rda"))
  cat(paste0("step 3 with ",to[1]-from[1]+1, " series is done\n"))
  
  res4 <- mclapply.hack(from[4]:to[4], function(n){
    cat(paste0(n,"\n"))
    temp <- tryCatch(optimizator(n), error=function(e){rep(NA,length(lst$inits1))})
    return(temp)})
  
  save(res4, file=paste0("output/d_",yearOrRegion,"_",outname,"_res4.rda"))
  cat(paste0("step 4 with ",to[1]-from[1]+1, " series is done\n"))
  
  res5 <- mclapply.hack(from[5]:to[5], function(n){
    cat(paste0(n,"\n"))
    temp <- tryCatch(optimizator(n), error=function(e) {rep(NA,length(lst$inits1))})
    return(temp)})
  
  save(res5, file=paste0("output/d_",yearOrRegion,"_",outname,"_res5.rda"))
  cat(paste0("step 5 with ",to[1]-from[1]+1, " series is done\n"))
  
  # 
  res <- c(res1,res2,res3,res4,res5)
  res <- res[1:ncol(pPanel)]
  
  save(res, file=paste0("output/d_",yearOrRegion,"_",outname,"_resALL.rda"))
  
  return(res)
  
}
