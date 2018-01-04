source("dlvPath.R")

plotRejs <- function(yearOrRegion){
  # Readang data
  if(is.numeric(yearOrRegion)){year      <-yearOrRegion
  filename  <- paste0("data/madisonFrom-",year,".csv")
  z         <- read.table(filename,header = T,sep = ";")
  z         <- data.matrix(z)} else
  {fname1 <- paste0("data/mds_G7-1950.csv"); fname2 <- paste0("data/mds_Europe-1950.csv"); fname3 <- paste0("data/mds_S&P-1950.csv")
  z_g7         <- read.table(fname1,header = T,sep = ";"); z_g7 <- data.matrix(z_g7)
  z_eu         <- read.table(fname2,header = T,sep = ";"); z_eu <- data.matrix(z_eu)
  z_sp         <- read.table(fname3,header = T,sep = ";"); z_sp <- data.matrix(z_sp)
  
  
  if(yearOrRegion=="Europe+G7"|yearOrRegion=="G7+Europe")
  {z<- matrix(c(z_g7,z_eu),dim(z_g7)[1],); colnames(z)<- c(colnames(z_g7),colnames(z_eu))} else
    if(yearOrRegion=="Europe+S&P"|yearOrRegion=="S&P+Europe")
    {z<- matrix(c(z_eu,z_sp),dim(z_eu)[1],); colnames(z)<- c(colnames(z_eu),colnames(z_sp))} else
      if(yearOrRegion=="G7+S&P"|yearOrRegion=="S&P+G7")
      {z<- matrix(c(z_g7,z_sp),dim(z_g7)[1],); colnames(z)<- c(colnames(z_g7),colnames(z_sp))} else 
        if(yearOrRegion=="Maddison"){
          filename  <- "data/madisonFromMaxNA2-1950.csv"
          z         <- read.table(filename,header = T,sep = " ")
          z         <- data.matrix(z)
        } else
        {stop("Unknown Country List")}
  if(sum(duplicated(t(z)))!=0){z<- z[,-which(duplicated(t(z)))]}
  }
  
  # Processing data
  z      <- log(z)
  n      <- ncol(z)
  cmbn   <- combn(n,2)
  cnames <- colnames(z)
  
  pPanel <- sapply(1:(n*(n-1)/2), function(i) z[,cmbn[1,i]] - z[,cmbn[2,i]])
  
  res      <- get(load(paste0("results/d_",yearOrRegion,"_resALL.rda")))
  allPaths <- lapply(1:ncol(pPanel), function(n) dlvPath(res[n,1:6], pPanel[,n]))
  oneRouth <- t(sapply(allPaths, function(p) apply(p,1,sum)==0))
  inds     <- cbind(oneRouth,oneRouth,matrix(rep(FALSE,n*(n-1)),,2))
  
  ds <- sapply(1:ncol(pPanel), function(i){
    # 1 if estimated d is less than 0, 1 otherwise
    d <- allPaths[[i]][2,]
    d[d==0] <- res[i,1]; d[d==1] <- res[i,2]
    return(d)
  })

  # rates of pairs with d<1
  rejs <- apply(ds,1, function(d) mean(d<1))
  
  # Plotting
  year <- if(is.numeric(yearOrRegion)){yearOrRegion:2010} else {1950:2010}
  heading <- paste0(yearOrRegion," Rejection Rates")
  pdf(paste0("plots/",yearOrRegion,"rejRates.pdf"))
  
  plot(year, rejs, main=heading,pch=16)

  dev.off()
  
  return(plot(year, rejs,main=heading))
}
    

  
