
# Converting outputs to table
totable <- function(res){
  t(sapply(res, function(r) c(r$par,r$value)))
}

# Reshaping tables into similar forms
correctRes <- function(res){
  # Reformats the table of estimates for it to be in "dH,dL,PHH,PLL, ..." order
  if(ncol(res) == 7){
    res[res[,1]>res[,2],] <- res[res[,1]>res[,2],][,c(2,1,4,3,5:ncol(res))]
  } else {
    res[res[,1]>res[,2],] <- res[res[,1]>res[,2],][,c(2,1,4,3,5,7,6,8)]}
  
  
  return(res)
}


report <- function(res){
  # reports the statistics of the output table
  nr        <- nrow(res)
  noChange  <- res[res[,1]==res[,2],]
  noChangeS <- noChange[noChange[,1]<1,]
  noChangeN <- noChange[noChange[,1]>=1,]
  change    <- res[res[,1]!=res[,2],]
  changeSS  <- change[change[,1] < 1 & change[,2] < 1,]
  changeNN  <- change[change[,1] >= 1 & change[,2] >= 1,]
  changeNS  <- change[(change[,1] >= 1 & change[,2] <= 1) | (change[,1] <= 1 & change[,2] >= 1),]
  
  repp <- t(sapply(list(noChangeS,noChangeN,changeSS,changeNN,changeNS), function(d)
    c(nrow(d)/nr, apply(d[,c(1:4)],2,function(x) mean(x,na.rm=T)))))
  repp <- repp[c(3,5,4),]
  rownames(repp) <- c("C - C","C - D","D - D")
  colnames(repp) <- c("Perc.","Avg. d_1","Avg. d_2","Avg. P_11","Avg. P_22", "Rej. ADF")
  
  return(repp)
}

gen_pdat <- function(yearOrRegion){
  data_file <- "data/"
  if(suppressWarnings(!is.na(as.numeric(yearOrRegion)))){year      <- yearOrRegion
  filename  <- paste0(data_file, "madisonFrom-",year,".csv")
  z         <- read.table(filename,header = T,sep = ";")
  z         <- data.matrix(z)} else if(yearOrRegion=="Maddison"){
    filename  <- paste0(data_file,"madisonFromMaxNA2-1950.csv")
    z         <- read.table(filename,header = T,sep = " ")
    z         <- data.matrix(z)
  } else if(yearOrRegion=="Europe"){
    fname <- paste0(data_file,"mds_Europe-1950.csv")
    z     <- read.table(fname,header = T,sep = ";"); z <- data.matrix(z)
  } else {
    fname1 <- paste0(data_file,"mds_G7-1950.csv")
    fname2 <- paste0(data_file,"mds_Europe-1950.csv")
    fname3 <- paste0(data_file,"mds_S&P-1950.csv")
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
  
  
  pdat   <- apply(combn(ncol(z),2),2, function(x) log(z[,x[1]]) - log(z[,x[2]]))
  pnames <- apply(combn(colnames(z),2), 2, function(n) paste0(sort(n),collapse = " - ")) 
  
  
  colnames(pdat) <- pnames
  return(pdat)
  
}

report <- function(res){
  # reports the statistics of the output table
  ind        <- (res[,1] < res[,2]) & (res[,1]!=0 | res[,2]!=0)
  if(ncol(res)==7) {res[ind,] <- res[ind, c(2,1,4,3,5,7,6)]} else {
    res[ind,] <- res[ind, c(2,1,4,3,5,6)]
  }
  
  nr         <- nrow(res)
  noChange   <- res[(res[,1]==0) | (res[,2]==0),]
  noChangeS  <- noChange[noChange[,1]<1,]
  noChangeN  <- noChange[noChange[,1]>=1,]
  change     <- res[(res[,1]!=0) & (res[,2]!=0),]
  changeSS   <- change[change[,1] < 1 & change[,2] < 1,]
  changeNN   <- change[change[,1] >= 1 & change[,2] >= 1,]
  changeNS   <- change[(change[,1] >= 1 & change[,2] <= 1) | (change[,1] <= 1 & change[,2] >= 1),]
  
  
  repp <- t(sapply(list(noChangeS,noChangeN,changeSS,changeNN,changeNS), function(d){
    d[d==0] <- NA
    c(nrow(d)/nr, apply(d[,c(1:ncol(res))],2,function(x) mean(x,na.rm=T))
    )
  }))
  repp <- repp[c(1,2,3,5,4),]
  rownames(repp) <- c("C","D","C - C","C - D","D - D")
  colnames(repp) <- c("Perc.","Avg. d_1","Avg. d_2","Avg. P_11","Avg. P_22", "Avg. Sigma", "Avg. mu_1", "Avg. mu_2")[1:(ncol(res)+1)]
  repp <- round(repp,3)
  return(repp)
}
