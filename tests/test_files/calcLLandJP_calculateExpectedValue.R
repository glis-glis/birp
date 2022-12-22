calcLLandJPforTest <- function(gamma, data, rho){
  #calculate p_ik
  change <- rho %*% gamma;
  change <- exp(change - mean(change));
  sigma <- as.matrix(data$efforts * rep(change, each=data$num.data.sets));
  sigma <- sigma / mean(sigma);
  p <- sigma / rowSums(sigma); #normalize
  
  #calculate likelihood
  logLik <- 0;
  for(i in 1:data$num.data.sets){
    logLik <- logLik + sum(data$counts[i,]*log(p[i,]))
  }
  
  #calc prior density
  n <- rowSums(data$counts);
  num.epochs <- length(gamma);
  I <- matrix(0, ncol=num.epochs, nrow=num.epochs);
  
  for(i in 1:data$num.data.sets){
    denom <- sum(sigma[i,])^2;
    for(k in 1:data$num.time.points){
      if(data$efforts[i,k] > 0){
        r <- rep(0, num.epochs);
        for(e in 1:num.epochs){
          r[e] <- sum((rho[k,e] - rho[,e])*sigma[i,]);
        }
        r <- r * sigma[i,k] / denom
        I <- I + r %*% t(r) * n[i] / p[i,k];
      }
    }
  }
  
  suppressWarnings(logPrior <- log(det(I)));
  if(is.na(logPrior) | is.infinite(logPrior)){
    logPrior <- log(det(I+diag(10^-100, dim(I)[1])));
  }
  
  #return posterior density
  return(c(logLik, logPrior));
}

set.seed(1000)
# test 1
gamma = c(-2,1.5)
x <- simulate_birp(efforts = 100, gamma = gamma, times = 1:4, times.of.change = 2)
data <- x
# calc rho
rho <- .calculateRho.birp_epochs(times = 1:4, times.of.change = 2)
# calc LL and JP in log
calcLLandJPforTest(gamma=gamma, data = data, rho = rho)

# test 2
gamma = c(-1,1)
x <- simulate_birp(efforts = 100, gamma = gamma, times = 1:10, times.of.change = 4)
data <- x
# calc rho
rho <- .calculateRho.birp_epochs(times = 1:10, times.of.change = 4)
# calc LL and JP in log
calcLLandJPforTest(gamma=gamma, data = data, rho = rho)
