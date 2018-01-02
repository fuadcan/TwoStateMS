data_file <- "output"

source("dlvPath.R")
gen_pdat <- function(yearOrRegion){
  
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
  pnames <- apply(combn(colnames(z),2), 2, function(n) paste0(n,collapse = " - ")) 
  
  
  colnames(pdat) <- pnames
  return(pdat)
  
}

# Read outputs
results <- dir("output/",pattern = "resALL.rda")
resss   <- lapply(results, function(r) get(load(paste0("output/",r))))
dnames  <- gsub("d_|_resALL.rda","",results)
pdats   <- lapply(dnames, function(n) gen_pdat(n))

pathss   <- lapply(1:length(resss), function(i) lapply(1:length(resss[[i]]), function(x) dlvPath(resss[[i]][[x]]$par,pdats[[i]][,x])))
ischange <- lapply(pathss, function(paths) t(sapply(paths, function(p) {temp <- apply(p,1,sum) > 0; temp <- c(temp,temp,T,temp); return(temp)})))
parss    <- lapply(resss,  function(ress)  t(sapply(ress, function(res) head(res$par,7))))
parss    <- lapply(1:length(resss),  function(i)  parss[[i]] * ischange[[i]])

# res <- parss[[1]]
report <- function(res){
  # reports the statistics of the output table
  ind        <- (res[,1] < res[,2]) & (res[,1]!=0 | res[,2]!=0)
  res[ind,] <- res[ind, c(2,1,4,3,5,7,6)]
  
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
      c(nrow(d)/nr, apply(d[,c(1:7)],2,function(x) mean(x,na.rm=T))
        )
    }))
  repp <- repp[c(1,2,3,5,4),]
  rownames(repp) <- c("C","D","C - C","C - D","D - D")
  colnames(repp) <- c("Perc.","Avg. d_1","Avg. d_2","Avg. P_11","Avg. P_22", "Avg. Sigma", "Avg. mu_1", "Avg. mu_2")
  repp <- round(repp,3)
  return(repp)
}

reports <- lapply(parss, report)
reports <- lapply(1:length(dnames), function(i) {temp <- cbind(dnames[i], reports[[i]])
                                                colnames(temp)[1] <- "Data"; return(temp)})
reports <- do.call(rbind,reports)
write.csv(reports, "~/Documents/twostate/results/reports.csv")


results <- dir("~/Documents/twostate_d_results/",pattern = "resALL.rda")
resss   <- lapply(results, function(r) get(load(paste0("~/Documents/twostate_d_results/",r))))
dnames  <- gsub("d_|_resALL.rda","",results)
pdats   <- lapply(dnames, function(n) gen_pdat(n))

pathss   <- lapply(1:length(resss), function(i) lapply(1:nrow(resss[[i]]), function(x) dlvPath(resss[[i]][x,],pdats[[i]][,x])))
ischange <- lapply(pathss, function(paths) t(sapply(paths, function(p) {temp <- apply(p,1,sum) > 0; temp <- c(temp,temp,T,T); return(temp)})))
parss    <- lapply(resss,  function(ress)  t(apply(ress, 1, function(res) head(res,6))))
parss    <- lapply(1:length(resss),  function(i)  parss[[i]] * ischange[[i]])

# res <- parss[[1]]
report <- function(res){
  # reports the statistics of the output table
  ind        <- (res[,1] < res[,2]) & (res[,1]!=0 | res[,2]!=0)
  res[ind,] <- res[ind, c(2,1,4,3,5,7,6)]
  
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
    c(nrow(d)/nr, apply(d[,c(1:7)],2,function(x) mean(x,na.rm=T))
    )
  }))
  repp <- repp[c(1,2,3,5,4),]
  rownames(repp) <- c("C","D","C - C","C - D","D - D")
  colnames(repp) <- c("Perc.","Avg. d_1","Avg. d_2","Avg. P_11","Avg. P_22", "Avg. Sigma", "Avg. mu_1", "Avg. mu_2")
  repp <- round(repp,3)
  return(repp)
}


reports <- lapply(parss, report)
reports <- lapply(1:length(dnames), function(i) {temp <- cbind(dnames[i], reports[[i]])
colnames(temp)[1] <- "Data"; return(temp)})
reports <- do.call(rbind,reports)
write.csv(reports, "~/Documents/twostate/results/reports.csv")






lapply(1:length(resss), function(i) {temp <- ischange[[i]] * parss[[i]]; t})
t(parss[[1]])
dss     <- lapply(1:length(resss), function(i) sapply(1:length(resss[[i]]), function(x) matrix(resss[[i]][[x]]$par[1:2],1,) %*% dlvPath(resss[[i]][[x]]$par,pdats[[i]][,x])))

d_plots     <- lapply(dss, function(ds) apply(ds<1,1,mean)) 

for(i in 1:length(d_plots)){
  titl <- paste0(dnames[i]," ","Rejection Rates")
  png(paste0("Documents/twostate/plots/",dnames[i],"_rejRates_DM.png"))
  plot((2010 - length(d_plots[[i]])+ 1):2010, d_plots[[i]],xlab="year",ylab="rejs",main=titl)
  dev.off()
}
# n1 <- "Finland"; n2 <- "Greece"; datname <- 1930

plot_specific <- function(n1,n2,datname){
  pairname <- paste0(n1," - ",n2)
  datind   <- which(datname==dnames)
  pnames   <- colnames(gen_pdat(datname))
  serind   <- which(pairname == pnames)
  ser      <- dss[[datind]][,serind]; maintext <- paste0(pnames[serind], " ", dnames[datind])
  pdf(paste0("Documents/twostate/plots/",dnames[datind],"_",pnames[serind],".pdf"))
  plot((2010 - length(ser)+ 1):2010, ser,type="l",xlab="year",ylab="d",main=maintext)
  dev.off()
  plot((2010 - length(ser)+ 1):2010, ser,type="l",xlab="year",ylab="d",main=maintext)
}




for(i in 1:length(dnames)){
  pdf(paste0("~/Documents/twostate/plots/d_lt_1_",dnames[i],".pdf"))
  plot((2010 - length(d_plots[[i]]) + 1):2010, d_plots[[i]],ylab="prop(d < 1)",xlab="year",main=dnames[i])
  dev.off()
}

datind <- 4; pnames <- colnames(gen_pdat(dnames[datind]))
serind <- 1

shiftcount <- sapply(1:length(pnames), function(serind) sum(dss[[datind]][,serind] != c(tail(dss[[datind]][,serind],-1),0))-1)

for (serind in (1:length(pnames))[shiftcount>=2]) {
  ser <- dss[[datind]][,serind]; maintext <- paste0(pnames[serind], ", ", dnames[datind])
  pdf(paste0("Documents/twostate/plots/",dnames[datind],"_",pnames[serind],"_scount", shiftcount[serind],".pdf"))
  plot((2010 - length(ser)+ 1):2010, ser,type="l",xlab="year",ylab="d",main=maintext)
  dev.off()
}

convendcount <- sapply(1:length(pnames), function(serind) tail(dss[[datind]][,serind],1) < 1)

for (serind in (1:length(pnames))[convendcount]) {
  ser <- dss[[datind]][,serind]; maintext <- paste0(pnames[serind], ", ", dnames[datind])
  pdf(paste0("Documents/twostate/plots/",dnames[datind],"_",pnames[serind],"_convend",".pdf"))
  plot((2010 - length(ser)+ 1):2010, ser,type="l",xlab="year",ylab="d",main=maintext)
  dev.off()
}

pnames <- colnames(gen_pdat(dnames[datind]))

for (serind in (1:length(pnames))[shiftcount>=3]) {
  ser <- dss[[datind]][,serind]; maintext <- paste0(pnames[serind], ", ", dnames[datind])
  png(paste0("Documents/twostate/plots_png/",dnames[datind],"_",pnames[serind],"_scount", shiftcount[serind],".png"))
  plot((2010 - length(ser)+ 1):2010, ser,type="l",xlab="year",ylab="d",main=maintext)
  dev.off()
}

convendcount <- sapply(1:length(pnames), function(serind) tail(dss[[datind]][,serind],1) < 1)

for (serind in (1:length(pnames))[convendcount]) {
  ser <- dss[[datind]][,serind]; maintext <- paste0(pnames[serind], ", ", dnames[datind])
  png(paste0("Documents/twostate/plots_png/",dnames[datind],"_",pnames[serind],"_convend",".png"))
  plot((2010 - length(ser)+ 1):2010, ser,type="l",xlab="year",ylab="d",main=maintext)
  dev.off()
}



# res <- parss[[6]]
reportD <- function(res){
  # reports the statistics of the output table
  ind        <- (res[,1] < res[,2]) & (res[,1]!=0 | res[,2]!=0)
  res[ind,] <- res[ind, c(2,1,4,3,5,6)]
  
  nr         <- nrow(res)
  noChange   <- res[(res[,1]==0) | (res[,2]==0),]
  noChangeS  <- noChange[noChange[,1]<1,]
  noChangeN  <- noChange[noChange[,1]>=1,]
  change     <- res[(res[,1]!=0) & (res[,2]!=0),]
  changeSS   <- matrix(change[change[,1] < 1 & change[,2] < 1,],,6)
  changeNN   <- change[change[,1] >= 1 & change[,2] >= 1,]
  changeNS   <- change[(change[,1] >= 1 & change[,2] <= 1) | (change[,1] <= 1 & change[,2] >= 1),]
  
  
  repp <- t(sapply(list(noChangeS,noChangeN,changeSS,changeNN,changeNS), function(d){
    d[d==0] <- NA
    d <- matrix(d[,1:6],,6)
    c(nrow(d)/nr, apply(d,2, function(x) mean(x,na.rm=T))
    )
  }))
  repp <- repp[c(1,2,3,5,4),]
  rownames(repp) <- c("C","D","C - C","C - D","D - D")
  colnames(repp) <- c("Perc.","Avg. d_1","Avg. d_2","Avg. P_11","Avg. P_22", "Avg. Sigma", "Avg. mu_1")
  repp <- round(repp,3)
  return(repp)
}

reports <- lapply(parss, reportD)
reports <- lapply(1:length(dnames), function(i) {temp <- cbind(dnames[i], reports[[i]])
colnames(temp)[1] <- "Data"; return(temp)})
reports <- do.call(rbind,reports)
write.csv(reports, "~/Documents/twostate_d_results/reports.csv")
