source("lnviD2.R")
source("lnviDM2.R")
source("lnviDSM2.R")
source("lnviDS2.R")

# Initial and Boundary Parameters
pars <- list(
  DMS = list(lowerV = c(-2,-2,.8,.8,0.001,0.001,-5,-5),
             upperV = c(2,2,.999,.999,5,5,5,5),
             inits1 = c(1.2,  1.7,  0.9,  0.9,  0.01,0.01, 0.1,0.1),
             inits2 = c(.8,  1.2,  0.9,  0.9,  0.01,0.01, 0.1,0.1),
             inits3 = c(.5,  1,  0.9,  0.9,  0.01,0.01, 0.1,0.1)),
  DM=list(lowerV = c(-2,-2,.8,.8,0.001,-5,-5),
          upperV = c(2,2,.999,.999,5,5,5),
          inits1 = c(1.2,  1.7,  0.9,  0.9,  0.01, 0.1,0.1),
          inits2 = c(.8,  1.2,  0.9,  0.9,  0.01, 0.1,0.1),
          inits3 = c(.5,  1,  0.9,  0.9,  0.01, 0.1,0.1)),
  DS  = list(lowerV = c(-2,-2,.8,.8,0.001,0.001,-5),
             upperV = c(2,2,.999,.999,5,5,5),
             inits1 = c(1.2,  1.7,  0.9,  0.9,  0.01,0.01, 0.1),
             inits2 = c(.8,  1.2,  0.9,  0.9,  0.01,0.01, 0.1),
             inits3 = c(.5,  1,  0.9,  0.9,  0.01,0.01, 0.1)),
  D  = list(lowerV = c(-2,-2,.8,.8,0.001,-5),
            upperV = c(2,2,.999,.999,5,5),
            inits1 = c(1.2,  1.7,  0.9,  0.9,  0.01, 0.1),
            inits2 = c(.8,  1.2,  0.9,  0.9,  0.01, 0.1),
            inits3 = c(.5,  1,  0.9,  0.9,  0.01, 0.1))
  
)

# Algorithms
DLValgo <- list(DMS=lnviDMS2,DM=lnviDM2,DS=lnviDS2,D=lnviD2)


MSArfima.fit <- function(ser, type){
### Inputs
## ser is the series to fit
## type is the type of the algorithm, whether:
# (d,mu,sigma) state switching (type="DMS") or
# (d,mu)       state switching (type="DM")  or
# (d,sigma)    state switching (type="DS")  or
# (d)          state switching (type="D").


# Parameters
DLValgo  <- DLValgo[type]
pars     <- pars[type]

# Constraints
const_mat <- matrix(0,length(pars$lowerV),length(pars$lowerV))
diag(const_mat) <- 1
const_mat <- rbind(const_mat,-const_mat)
const_mat <- cbind(const_mat,c(pars$lowerV,-pars$upperV))

# Function for optimizing DLValgo
optimizator <- function(i){
  inits1  <- pars$inits1
  temp1 <- constrOptim(inits1, function(p) -DLValgo(p,pPanel[,i]), NULL, ui = const_mat[,-ncol(const_mat)], const_mat[,ncol(const_mat)])
  if(temp1$convergence!=0){out <- temp1} else {
    inits2 <- pars$inits2
    temp2  <- constrOptim(inits2, function(p) -DLValgo(p,pPanel[,i]), NULL, ui = const_mat[,-ncol(const_mat)], const_mat[,ncol(const_mat)])
    if(temp2$convergence != 0) {out <- temp2} else {
      inits3 <- pars$inits3
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
}
