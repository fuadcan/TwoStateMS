data_file <- "~/Documents/twostate/data/"
source("Documents/threestate/return_path.R")
gen_pdat <- function(yearOrRegion){
  
  if(suppressWarnings(!is.na(as.numeric(yearOrRegion)))){year      <- yearOrRegion
  filename  <- paste0(data_file, "madisonFrom-",year,".csv")
  z         <- read.table(filename,header = T,sep = ";")
  z         <- data.matrix(z)} else if(yearOrRegion=="Maddison"){
    filename  <- paste0(data_file,"madisonFromMaxNA2-1950.csv")
    z         <- read.table(filename,header = T,sep = " ")
    z         <- data.matrix(z)
  } else {fname1 <- paste0(data_file,"mds_G7-1950.csv")
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

results <- dir("~/Documents/threestate/modified_results/",pattern = "modress")
resss   <- lapply(results, function(r) get(load(paste0("~/Documents/threestate/modified_results/",r))))
dnames  <- gsub("_modress.rda","",results)
pdats   <- lapply(dnames, function(n) gen_pdat(n))

dss     <- lapply(1:length(resss), function(i) sapply(1:length(resss[[i]]), function(x) matrix(resss[[i]][[x]]$par[1:3],1,) %*% return_path(resss[[i]][[x]]$par,pdats[[i]][,x])))

d_plots     <- lapply(dss, function(ds) apply(ds<1,1,mean)) 
s_plots     <- lapply(dss, function(ds) apply(ds<.5,1,mean)) 
mr_plots    <- lapply(dss, function(ds) apply(ds>=.5 & ds<1,1,mean))
nc_plots    <- lapply(dss, function(ds) apply(ds>=1,1,mean))
conf_plots  <- lapply(dss, function(ds) apply(ds<1,1,mean))


for(i in 1:length(dnames)){
  pdf(paste0("~/Documents/threestate/plots/stat_",dnames[i]),".pdf")
  plot((2010 - length(d_plots[[i]]) + 1):2010, s_plots[[i]],ylab="prop(d < 0.5)",xlab="year",main=dnames[i])
  dev.off()
}

for(i in 1:length(dnames)){
  pdf(paste0("~/Documents/threestate/plots/meanrev_",dnames[i],".pdf"))
  plot((2010 - length(d_plots[[i]]) + 1):2010, mr_plots[[i]],ylab="prop(0.5 <= d <1)",xlab="year",main=dnames[i])
  dev.off()
}



for(i in 1:length(dnames)){
  pdf(paste0("~/Documents/threestate/plots/nc_",dnames[i],".pdf"))
  plot((2010 - length(d_plots[[i]]) + 1):2010, nc_plots[[i]],ylab="prop(d >= 1)",xlab="year",main=dnames[i])
  dev.off()
}

for(i in 1:length(dnames)){
  pdf(paste0("~/Documents/threestate/plots/conv_",dnames[i],".pdf"))
  plot((2010 - length(d_plots[[i]]) + 1):2010, conv_plots[[i]],ylab="prop(d < 1)",xlab="year",main=dnames[i])
  dev.off()
}

source("~/Documents/threestate/return_pairs.R")
# for(i in 1:length())
datind <- 1; pnames <- return_pairs(dnames[datind])
serind <- 1

# serind <- which(pnames == "Austria - Germany")
# ser <- dss[[datind]][,serind]; maintext <- paste0(pnames[serind], ", ", dnames[datind])
# pdf(paste0("~/Documents/ThreeStateConvergence/Sep_14_2017/plots/",dnames[datind],"_",pnames[serind],".pdf"))
# plot((2010 - length(ser)+ 1):2010, ser,type="l",xlab="year",ylab="d",main=maintext)
# dev.off()


shiftcount <- sapply(1:length(pnames), function(serind) sum(dss[[datind]][,serind] != c(tail(dss[[datind]][,serind],-1),0))-1)

for (serind in (1:length(pnames))[shiftcount>=5]) {
  ser <- dss[[datind]][,serind]; maintext <- paste0(pnames[serind], ", ", dnames[datind])
  pdf(paste0("plots/",dnames[datind],"_",pnames[serind],"_scount", shiftcount[serind],".pdf"))
  plot((2010 - length(ser)+ 1):2010, ser,type="l",xlab="year",ylab="d",main=maintext)
  dev.off()
}

convendcount <- sapply(1:length(pnames), function(serind) tail(dss[[datind]][,serind],1) < 1)

for (serind in (1:length(pnames))[convendcount]) {
  ser <- dss[[datind]][,serind]; maintext <- paste0(pnames[serind], ", ", dnames[datind])
  pdf(paste0("plots/",dnames[datind],"_",pnames[serind],"_convend",".pdf"))
  plot((2010 - length(ser)+ 1):2010, ser,type="l",xlab="year",ylab="d",main=maintext)
  dev.off()
}

datind <- 4; pnames <- return_pairs(dnames[datind])


for (serind in (1:length(pnames))[shiftcount>=5]) {
  ser <- dss[[datind]][,serind]; maintext <- paste0(pnames[serind], ", ", dnames[datind])
  png(paste0("plots_png/",dnames[datind],"_",pnames[serind],"_scount", shiftcount[serind],".png"))
  plot((2010 - length(ser)+ 1):2010, ser,type="l",xlab="year",ylab="d",main=maintext)
  dev.off()
}

convendcount <- sapply(1:length(pnames), function(serind) tail(dss[[datind]][,serind],1) < 1)

for (serind in (1:length(pnames))[convendcount]) {
  ser <- dss[[datind]][,serind]; maintext <- paste0(pnames[serind], ", ", dnames[datind])
  png(paste0("plots_png/",dnames[datind],"_",pnames[serind],"_convend",".png"))
  plot((2010 - length(ser)+ 1):2010, ser,type="l",xlab="year",ylab="d",main=maintext)
  dev.off()
}

pathss   <- lapply(1:length(resss), function(i) lapply(1:length(resss[[i]]), function(x) return_path(resss[[i]][[x]]$par,pdats[[i]][,x])))
ischange <- lapply(pathss, function(paths) t(sapply(paths, function(p) {temp <- apply(p,1,sum) > 0; temp <- c(temp,temp,temp,T,temp); return(temp)})))
parss    <- lapply(resss,  function(ress)  t(sapply(ress, function(res) res$par)))
parss    <- lapply(1:length(resss),  function(i)  parss[[i]] * ischange[[i]])



# res <- parss[[1]]

report <- function(res){
  # reports the statistics of the output table
  res[res==0] <- NA
  ranks <- t(apply(res[,1:3],1,function(r) order(r,decreasing = T,na.last = T)))
  
  ranks <- cbind(ranks,ranks+3,ranks+6,10,ranks+10)
  res <- t(sapply(1:nrow(ranks), function(r) res[r,ranks[r,]]))
  # head(res,20)
  
  nr          <- nrow(res)
  ind         <- apply(is.na(res[,1:3]),1,sum) == 2
  noChange    <- res[ind,]
  noChangeS   <- noChange[noChange[,1]<.5,]
  noChangeM   <- noChange[noChange[,1]>=.5 & noChange[,1]<1,]
  noChangeN   <- noChange[noChange[,1]>=1,]
  change      <- res[!ind,]
  change_dmin <- apply(change[,1:3],1,function(c) min(c,na.rm = T))
  change_dmax <- apply(change[,1:3],1,function(c) max(c,na.rm = T))
  changeSS    <- change[change_dmin < .5 & change_dmax < .5,]
  changeNN    <- change[change_dmin >= 1 & change_dmax >= 1,]
  changeMM    <- change[change_dmin >= .5 & change_dmin < 1 & change_dmax >= .5 & change_dmax < 1,]
  changeNS    <- change[change_dmin < .5 & change_dmax >= 1,]
  changeMS    <- change[change_dmin < .5 & change_dmax < 1 & change_dmax >= .5,]
  changeNM    <- change[change_dmin >= .5 & change_dmin < 1 & change_dmax >= 1,]
  
  
  caselist <- list(noChangeS,noChangeM,noChangeN,changeSS, changeNS ,changeMS,changeMM,changeNM,changeNN)
  if(any(sapply(caselist,class)!="matrix")){
    caselist[sapply(caselist,class)!="matrix"] <- lapply(caselist[sapply(caselist,class)!="matrix"], function(c) matrix(c,,13))
  }
  repp <- t(sapply(caselist, function(d){
    d[d==0] <- NA
    d <- matrix(d[,1:13],,13)
    c(nrow(d)/nr, apply(d,2,function(x) mean(x,na.rm=T))
    )
  }))
  
  rownames(repp) <- c("S","M","N","S - S","N - S","M - S","M - M", "N - M", "N - N")
  colnames(repp) <- c("Perc.","Avg. d_1","Avg. d_2","Avg. d_3","Avg. P_11","Avg. P_22","Avg. P_33","Avg. P_12","Avg. P_21","Avg. P_31", "Avg. Sigma", "Avg. mu_1", "Avg. mu_2", "Avg. mu_3")
  
  return(repp)
}

reports <- lapply(parss, report)
reports <- lapply(1:length(dnames), function(i) {temp <- cbind(dnames[i], reports[[i]])
colnames(temp)[1] <- "Data"; return(temp)})
reports <- do.call(rbind,reports)

write.csv(reports, "~/Documents/threestate/results/reports.csv")

