setwd("~/TwoStateMS/")
source("convDLV.R")
# source("plotAll.R")
# source("plotRejs.R")

# Estimating d's for only d switching analysis
res1930   <- convDLV(1930,1)
Sys.time()
resg7sp   <- convDLV("G7+S&P",1)
Sys.time()
Sys.sleep(300)
reseu     <- convDLV("Europe",1)
Sys.time()
resg7eu   <- convDLV("Europe+G7",1)
Sys.time()
Sys.sleep(300)
Sys.time()
resspeu   <- convDLV("Europe+S&P",1)
Sys.time()
ress_d <- list(res1930,resg7sp,resg7eu,resspeu)

# Estimating d's for both d and mu switching analysis
res1930   <- convDLV(1930,2)
Sys.time()
resg7sp   <- convDLV("G7+S&P",2)
Sys.time()
Sys.sleep(300)
reseu     <- convDLV("Europe",2)
Sys.time()
resg7eu   <- convDLV("Europe+G7",2)
Sys.time()
Sys.sleep(300)
Sys.time()
resspeu   <- convDLV("Europe+S&P",2)
Sys.time()
ress_dm   <- list(res1930,resg7sp,resg7eu,resspeu)


correctRes <- function(res){
  # Reformats the table of estimates for it to be in "dH,dL,PHH,PLL, ..." order
  res[res[,1]>res[,2],] <- res[res[,1]>res[,2],][,c(2,1,4,3,5:ncol(res))]
  
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
    c(nrow(d)/nr, apply(d[,c(1:4)],2,function(x) mean(x,na.rm=T)), sum(d[,9])/nr)))
  repp <- repp[c(3,5,4),]
  rownames(repp) <- c("N - N","N - Y","Y - Y")
  colnames(repp) <- c("Perc.","Avg. d_1","Avg. d_2","Avg. P_11","Avg. P_22", "Rej. ADF")
  
  return(repp)
}

# Reading outputs
res1930 <- get(load("results/d_1930_resALL.rda"))
resg7eu <- get(load("results/d_Europe+G7_resALL.rda"))
resg7sp <- get(load("results/d_G7+S&P_resALL.rda"))
res1940 <- get(load("results/d_1940_resALL.rda"))
resspeu <- get(load("results/d_Europe+S&P_resALL.rda"))

# Procesing outputs
res1930 <- correctRes(res1930)
save(res1930,file="results/d_1930_resALL.rda")
resg7eu <- correctRes(resg7eu)
save(resg7eu,file="results/d_Europe+G7_resALL.rda")
resg7sp <- correctRes(resg7sp)
save(resg7sp,file="results/d_G7+S&P_resALL.rda")
res1940 <- correctRes(res1940)
save(res1940,file="results/d_1940_resALL.rda")
resspeu <- correctRes(resspeu)
save(resspeu,file="results/d_Europe+S&P_resALL.rda")


# Reporting outputs
rep <- suppressWarnings(data.frame(rbind(report(res1930),report(res1940),report(resg7eu),report(resg7sp),report(resspeu))))
rep <- cbind(unlist(lapply(c("1930","1940","g7+eu","g7+sp","eu+sp"),function(d) rep(d, 3))),rep)
rep <- cbind(unlist(lapply(c("S - S","N - Y","Y - Y"),function(d) rep(d, 5))),rep)
colnames(rep)[1:2] <- c("class","data")
write.table(rep, "table.csv", col.names = T, row.names = T)

# Plotting all path graphs
lapply(c(1930, 1940), plotAll)
lapply(c("G7+S&P","Europe+G7","Europe+S&P"),plotAll)
lapply(c(1930, 1940), plotRejs)
lapply(c("G7+S&P","Europe+G7","Europe+S&P"), plotRejs)
