setwd("~/TwoStateMS/")
source("convDLV.R")
source("utils.R")
source("dlvPath.R")
# source("plotAll.R")
# source("plotRejs.R")

# Estimating d's for only d switching analysis
res1930   <- convDLV(1930,1)
res1940   <- convDLV(1940,1)
res1950   <- convDLV(1950,1)
reseu     <- convDLV("Europe",1)
resg7eu   <- convDLV("Europe+G7",1)
resg7sp   <- convDLV("G7+S&P",1)
resspeu   <- convDLV("Europe+S&P",1)

# Estimating d's for both d and mu switching analysis
res1930   <- convDLV(1930,2)
res1940   <- convDLV(1940,2)
res1950   <- convDLV(1950,2)
reseu     <- convDLV("Europe",2)
resg7eu   <- convDLV("Europe+G7",2)
resg7sp   <- convDLV("G7+S&P",2)
resspeu   <- convDLV("Europe+S&P",2)




# Reading outputs
ress_d  <- lapply(dir("output","_D_.*resALL.rda"), function(d) get(load(paste0("output/",d))))
ress_dm <- lapply(dir("output","_DM_.*resALL.rda"), function(d) get(load(paste0("output/",d))))

# Converting to tables
ress_d  <- lapply(ress_d, totable)
ress_dm <- lapply(ress_dm, totable)

# Reformating results
ress_d  <- lapply(ress_d, correctRes)
ress_dm <- lapply(ress_dm, correctRes)

# data names and pair panels
dnames     <- gsub("d_|_D_resALL.rda","",dir("output","_D_.*resALL.rda"))
pdats      <- lapply(dnames, gen_pdat)

# State switching series for dm and d
pathss_DM   <- lapply(1:length(ress_dm), function(i) lapply(1:nrow(ress_dm[[i]]), function(x) dlvPath_dm(ress_dm[[i]][x,-8],pdats[[i]][,x])))
change_DM   <- lapply(pathss_DM, function(paths) sapply(paths, function(p) sum(apply(p,1,sum)>0)>1))

pathss_D    <- lapply(1:length(ress_d), function(i) lapply(1:nrow(ress_d[[i]]), function(x) dlvPath_d(ress_d[[i]][x,-7],pdats[[i]][,x])))
change_D    <- lapply(pathss_D, function(paths) sapply(paths, function(p) sum(apply(p,1,sum)>0)>1))

# Correcting results
ischange <- lapply(pathss_DM, function(paths) t(sapply(paths, function(p) {temp <- apply(p,1,sum) > 0; temp <- c(temp,temp,T,temp); return(temp)})))
parss_dm <- lapply(1:length(ress_dm),  function(i)  ress_dm[[i]][,-8] * ischange[[i]])

ischange <- lapply(pathss_D, function(paths) t(sapply(paths, function(p) {temp <- apply(p,1,sum) > 0; temp <- c(temp,temp,T,T); return(temp)})))
parss_d  <- lapply(1:length(ress_d),  function(i)  ress_d[[i]][,-7] * ischange[[i]])


# Reporting outputs
rep <- suppressWarnings(data.frame(rbind(report(res1930),report(res1940),report(res1940),report(resg7eu),report(resg7sp),report(resspeu))))
rep <- cbind(unlist(lapply(c("1930","1940","g7+eu","g7+sp","eu+sp"),function(d) rep(d, 3))),rep)
rep <- cbind(unlist(lapply(c("S - S","N - Y","Y - Y"),function(d) rep(d, 5))),rep)
colnames(rep)[1:2] <- c("class","data")
write.table(rep, "table.csv", col.names = T, row.names = T)

# Plotting all path graphs
lapply(c(1930, 1940), plotAll)
lapply(c("G7+S&P","Europe+G7","Europe+S&P"),plotAll)
lapply(c(1930, 1940), plotRejs)
lapply(c("G7+S&P","Europe+G7","Europe+S&P"), plotRejs)
