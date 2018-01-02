data_file <- "~/Documents/twostate/data/"

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
  pnames <- apply(combn(colnames(z),2), 2, function(n) paste0(sort(n),collapse = " - ")) 
  
  
  colnames(pdat) <- pnames
  return(pdat)
  
}

results_DM <- dir("~/Documents/twostate/results/",pattern = "resALL.rda")
resss_DM   <- lapply(results_DM, function(r) get(load(paste0("~/Documents/twostate/results/",r))))
dnames     <- gsub("d_|_resALL.rda","",results_DM)
pdats      <- lapply(dnames, function(n) gen_pdat(n))

source("~/Documents/twostate/dlvPath.R")
pathss_DM   <- lapply(1:length(resss_DM), function(i) lapply(1:length(resss_DM[[i]]), function(x) dlvPath(resss_DM[[i]][[x]]$par,pdats[[i]][,x])))
change_DM   <- lapply(pathss_DM, function(paths) sapply(paths, function(p) sum(apply(p,1,sum)>0)>1))

results_D <- dir("~/Documents/twostate_d_results/",pattern = "resALL.rda")
resss_D   <- lapply(results_D, function(r) get(load(paste0("~/Documents/twostate_d_results/",r))))
source("~/Documents/threestate/ReplciationFiles/dlvPath.R")
pathss_D  <- lapply(1:length(resss_D), function(i) lapply(1:nrow(resss_D[[i]]), function(x) dlvPath(resss_D[[i]][x,1:6],pdats[[i]][,x])))
change_D  <- lapply(pathss_D, function(paths) sapply(paths, function(p) sum(apply(p,1,sum)>0)>1))

results_1 <- dir("~/Documents/onestate/results/",pattern = "resALL.rda")
resss_1   <- lapply(results_1, function(r) get(load(paste0("~/Documents/onestate/results/",r))))

lklh_D    <- lapply(resss_D,  function(r) -r[,7])[-c(2,3)]
lklh_DM   <- lapply(resss_DM, function(ress) -sapply(ress, function(r) r$value))[-c(2,3)] 
lklh_1    <- lapply(resss_1,  function(ress) -sapply(ress, function(r) r$value)) 


pnames <- lapply(dnames[-c(2,3)], function(n) colnames(gen_pdat(n)))
lklhs     <- lapply(1:length(dnames[-c(2,3)]), function(i) {temp <- cbind(lklh_DM[[i]],lklh_D[[i]],lklh_1[[i]])
                                      colnames(temp) <- c("DM","D","one")
                                      rownames(temp) <- pnames[[i]]; temp})

temp <- resss_D[[3]]
rownames(temp) <- colnames(gen_pdat(1950))
needed_DM <- resss_DM[[3]][names(resss_DM[[3]]) %in% needed]
needed_D  <- temp[rownames(temp) %in% needed,] 
names(res_needed)
rownames(needed_D)
needed_DM   <- needed_DM[order(names(needed_DM))]
needed_D    <- needed_D[order(rownames(needed_D)),]
res_needed  <- res_needed[order(names(res_needed))]
names(res_needed)
rownames(needed_D)
names(needed_DM)

lklh_needed <- cbind(sapply(needed_DM,function(r) r$value), needed_D[,7] ,sapply(res_needed,function(r) r$value))
colnames(lklh_needed) <- c("DM","D","one")
write.csv(lklh_needed, file = paste0("~/Documents/twostate/results/lklhs_1950_needed.csv"))
write.csv(lklhs[[5]], file = paste0("~/Documents/twostate/results/lklhs_",dnames[-c(2,3)][5]))

source("~/Documents/onestate/lnviD2.R")

calculate_specific <- function(p1,p2,datname){
  pdat  <- gen_pdat(datname)
  pname <- paste0(sort(c(p1,p2)), collapse = " - ")
  ser   <- pdat[,colnames(pdat)==pname]
  
}

lklh_diff_DM <- sapply(1:length(lklh_DM), function(i) (lklh_DM[[i]][change_DM[[i]]] - lklh_1[[i]][change_DM[[i]]]))
lklh_diff_D  <- sapply(1:length(lklh_DM), function(i) (lklh_D[[i]][change_D[[i]]] - lklh_1[[i]][change_D[[i]]]))

pnames <- lapply(dnames, function(n) colnames(gen_pdat(n)))
lklh_diff_DM <- lapply(1:length(dnames), function(i) {temp <- lklh_diff_DM[[i]]; names(temp) <- pnamess[[i]][change_DM[[i]]]; temp})
lklh_diff_D  <- lapply(1:length(dnames), function(i) {temp <- lklh_diff_D[[i]]; names(temp) <- pnamess[[i]][change_D[[i]]]; temp})

save(lklh_diff_DM, file = "results/likelihood_diff_DM.rda")
save(lklh_diff_D , file = "results/likelihood_diff_D.rda")


# lklh_rep <- cbind(lklh_diff_DM,lklh_diff_D)

# lklhs     <- lapply(1:length(lklh_DM), function(i) cbind(lklh_DM[[i]],lklh_D[[i]],lklh_1[[i]]))
# lklh_diff <- lapply(lklhs, function(l) {temp <- abs(cbind(l[,1]-l[,2],l[,1]-l[,3],l[,2]-l[,3])); colnames(temp) <- c("DM-D","DM-1","D-1"); return(temp)})
# lklh_rep  <- t(sapply(lklh_diff, function(d) apply(d>3.84,2,mean)))

rownames(lklh_rep) <- dnames
write.csv(lklh_rep, "~/Documents/twostate/results/lklh_rep.csv")

mean(lklh_DM[[4]] > lklh_D[[4]])
mean(lklh_DM[[2]] > lklh_1[[2]])
mean(lklh_DM[[4]] <= lklh_D[[4]])

