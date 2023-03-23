# calc Rho test 1
times <- c(2015,2016,2017,2019,2020)
rho <- .calculateRho.birp_epochs(times = times, times.of.change = 2018)

# calc Rho test 2
times <- c(2000,2001,2002,2004,2005)
rho <- .calculateRho.birp_epochs(times = times, times.of.change = 2002)

# calc Rho test 3 - TOC same as two of the timepoints
times <- 1:9
rho <- .calculateRho.birp_epochs(times = times, times.of.change = c(3,7))

# calc Rho test 4 - TOC empty
times <- 1:5
rho <- .calculateRho.birp_epochs(times = times, times.of.change = NULL)


####################################################
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


# test NB LL & JP -> use file method_diggingInDirt_LLandJP1.csv
gamma <- c(-1, 1)
times.of.change = 3
mu <- c(0.2, 0.1, 0.5, 0.2)
a <- 1
effort <- rep(100,4)
rho <- .calculateRho.birp_epochs(1:4, times.of.change)
phi <- exp(rho %*% gamma)

# use equation from grant for DM distribution to get counts
library(dppmix)
alpha <- a * mu * phi * effort
set.seed(420)
probs <- rdirichlet(n=1, alpha)
# now use these probability to draw counts from a multinomial distr
rmultinom(n=1, size=100, prob=probs)
counts <- c(58, 2, 25, 15)

# calculate the density, but do it by hand, because we have to subtract the terms that drop out in the hastings ratio
LL <- log(gamma(counts + probs)) - log(gamma(probs)) - log(gamma())




