library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(MASS)
set.seed(99999)

### define data generation functions-with predictors
m <- 50
X=t(rbind(rep(1,m),rnorm(m,0,1),rbinom(m,1,0.2)))
beta=t(rbind(c(4,3,5),c(1,1,0)))

data.reg <-  function(tautheta=3,tausigma=1,
                      r1=0.3,r2=0.3, sigmas=1,...){
  # r1 - individual level
  # r2 - population level
  n <- dim(X)[1]
  theta = vector()
  sigma2 = vector()
  y = vector()
  s2 = vector()
  Sigma = list()
  temp = matrix(rep(NA),ncol = 2, nrow = n)
  indiv = matrix(rep(NA),ncol = 2, nrow = n)
  # generate individual true mean and variance from population
  mu <- X%*%beta
  varmat <- matrix(c(tautheta^2,r2*tausigma*tautheta,r2*tausigma*tautheta,tausigma^2), nrow=2)
  for (k in 1:n) temp[k,] <- MASS::mvrnorm(1,mu[k,],varmat)
  theta <- temp[,1]
  sigma2 <- (exp(temp[,2]))^2
  # generate observed value from true individual
  for (i in 1:n){
    Sigma[[i]] <- c(sigma2[i],r1*sqrt(sigma2[i])*sigmas,r1*sqrt(sigma2[i])*sigmas,sigmas^2)
    indiv[i,] <- mvrnorm(1, temp[i,],matrix(Sigma[[i]], nrow = 2))
  }
  y <- indiv[,1]
  s <- exp(indiv[,2])
  dataset <- data.frame(cbind(y, s, X))
  return(list(data = dataset,
              para = c(tautheta=tautheta, tausigma=tausigma, r1=r1, r2=r2, sigmas=sigmas)))
}

data <- data.reg(r1=0,r2=0)$data
res <- unimod(data$y, data$s, with_cov = TRUE, X = X,
              seed = 9999, iter = 1000, chain = 1, verbose = FALSE,
              warmup= 500, thin=1)$pop_par
