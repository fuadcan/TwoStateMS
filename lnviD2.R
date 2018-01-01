lnviD2 <- function(b,w){

  # b <- c(.4,.8,.99,.86,.4,.9,.1,-.3)
  
  size <- length(w)
  
  d1hat = b[1]      ; d2hat = b[2]
  p11hat = b[3]     ; p22hat = b[4];
  sigma1hat =  b[5] ; sigma2hat =  b[5];
  mu1hat = b[6]     ; mu2hat = b[6];
  
  rho1hat=0; theta1hat=0;
  
  tranmax = matrix(c(p11hat,1-p22hat,1-p11hat,p22hat),2) # CONTROL
  
  mu     <- c(mu1hat,mu2hat)
  sigma  <- c(sigma1hat,sigma2hat)
  coefd1 = matrix(cumprod((d1hat+(0:(size-2)))/(1:(size-1))),,1)
  coefd2 = matrix(cumprod((d2hat+(0:(size-2)))/(1:(size-1))),,1)
  coeff = cbind(coefd1,coefd2)
  
  # muhat is to store the selected path of state: from date 2 */
  eta1  = matrix(0,4,size); # row 1:2 for path1, 3:4 for path2 */
  eta   = matrix(0,4,size); # same as above */
  xihat = matrix(0,4,size); # same as above */
  ehat  = matrix(0,2,size); # row 1, for path1, row 2 for path2 */
  
  # to keep the long memory part, not iid */
  zsigmahat = matrix(0,2,size) # row 1, for path1, row 2 for path2 */
  zhat = matrix(0,2,size) # row 1, for path1, row 2 for path2 */
  lnlk <- matrix(0,2,size)
  
  path1 = matrix(0,2,size)
  path2 = matrix(0,2,size)
  path1[,1] <- c(1,0)
  path2[,1] <- c(0,1)
  #   path1[1,1:2] = matrix(c(1,0),1)
  #   path2[1,1:2] = matrix(c(0,1),1)
  zsigmahat[,1] = c(w[1] - mu[1], w[1] - mu[2])
  # zsigmahat is de-drifted series for each given path; or it is z*sigmahat
  
  zhat[,1] <- zsigmahat[,1] / sigma 
  
  ehat[,1] = zhat[,1]
  
  lnlk[,1] <- log(exp(-ehat[,1]^2/2 + 1e-50 ) / (sqrt(2*pi))*
                    c( (1-p22hat)/(2-p11hat-p22hat) +1e-50  , (1-p11hat)/(2-p11hat-p22hat)) +1e-50)
  
  j = 2
  
  
  ####################!!!!!!!!!!!!!!!!!!!!!
  while(j<=size){
    # tempsMUZ and temps are needed for likelihood possibility
    tempsMUZ <- c(t(matrix(c(rev(zsigmahat[1,1:(j-1)]),rev(zsigmahat[2,1:(j-1)])),,2))%*%coeff[1:(j-1),1:2])
    
    temps    <- rep(mu,2)[c(1,3,2,4)]+rho1hat*rep(sigma,2)[c(1,3,2,4)]*rep(zhat[1:2,j-1],2)+
      theta1hat*rep(sigma,2)[c(1,3,2,4)]*rep(ehat[1:2,j-1],2)
    
    # likelihood possibility
    eta1[,j]  <- exp(-(w[j] - temps - tempsMUZ)^2 / (2*rep(sigma,2)[c(1,3,2,4)]^2+1e-50))/
      (sqrt(2*pi*rep(sigma,2)[c(1,3,2,4)]^2))
    
    eta[,j] = eta1[,j] * c(tranmax)
    
    # Case is 1 for 1->1; 2 for 2->1; 3 for 1->2; 4 for 2->2
    case <- c(log(eta[1,j]+1e-50) + lnlk[1,j-1] > log(eta[2,j]+1e-50) + lnlk[2,j-1],
              log(eta[3,j]+1e-50) + lnlk[1,j-1] > log(eta[4,j]+1e-50) + lnlk[1,j-1])
    
    case <- c(sign(!case[1])+1,sign(!case[2])+3) # to convert cases to numbers
    
    # holding previous two path
    oldlnlk <- lnlk
    oldPath <- list(path1,path2)
    oldzsigmahat <- zsigmahat
    oldzhat <- zhat
    oldehat <- ehat
    
    # recalculating values
    # path ends up with state 1:
    lnlk[1,j]   <- log(eta[case[1],j]+1e-50)+oldlnlk[case[1],j-1]
    path1[,1:j] <- cbind(oldPath[[case[1]]][,1:(j-1)],c(1,0))
    zsigmahat[1,1:j] <- c(oldzsigmahat[case[1],1:(j-1)],  w[j] - mu[1] - tempsMUZ[case[1]])
    zhat[1,1:j] <- c(oldzhat[case[1],1:(j-1)], zsigmahat[1,j]/(sigma1hat))
    ehat[1,1:j] <- c(oldehat[case[1],1:(j-1)], (w[j] - tempsMUZ[case[1]] - temps[case[1]])/sigma1hat) # ONEMLI!!!
    
    # path ends up with state 2:
    lnlk[2,j]   = log(eta[case[2],j]+1e-50)+oldlnlk[case[2]-2,j-1];
    path2[,1:j] <- cbind(oldPath[[case[2]-2]][,1:(j-1)],c(0,1))
    zsigmahat[2,1:j] = c(oldzsigmahat[case[2]-2,1:(j-1)], w[j] - mu[2] - tempsMUZ[case[2]])
    zhat[2,1:j] = c(oldzhat[case[2]-2,1:(j-1)],zsigmahat[2,j]/(sigma2hat));
    ehat[2,1:j] = c(oldehat[case[2]-2,1:(j-1)],(w[j] - tempsMUZ[case[2]] -temps[case[2]])/sigma2hat)
    
    j=j+1
  }
  #   lnlipath1 = t(lnlk1);
  #   lnlipath2 = t(lnlk2);
  if(lnlk[1,size] > lnlk[2,size]){
    lnlkRes = lnlk[1,size]} else {lnlkRes = lnlk[2,size]}
  
  return(lnlkRes)
}
