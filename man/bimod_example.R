library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(MASS)
set.seed(99999)

### define data generation functions-no predictors
data.nonreg <-  function(n=50,mutheta=5,musigma=2,
                         tautheta=3, tausigma=1, r1=0.3, r2=0.3, sigmas=1,...){
  # mutheta,musigma,tautheta, tausigma -population level
  # r1 - individual level
  # r2 - population level
  theta = vector()
  sigma2 = vector()
  y = vector()
  s2 = vector()
  Sigma = list()
  indiv = matrix(rep(NA),ncol = 2, nrow = n)
  # generate individual true mean and variance from population
  mu <- c(mutheta,musigma)
  varmat <- matrix(c(tautheta^2,r2*tausigma*tautheta,r2*tausigma*tautheta,tausigma^2), nrow=2)
  temp <- MASS::mvrnorm(n, mu, varmat)
  theta <- temp[,1]
  sigma2 <- (exp(temp[,2]))^2
  # generate observed value from true individual
  for (i in 1:n){
    Sigma[[i]] <- c(sigma2[i],r1*sqrt(sigma2[i])*sigmas,r1*sqrt(sigma2[i])*sigmas,sigmas^2)
    indiv[i,] <- mvrnorm(1, temp[i,],matrix(Sigma[[i]], nrow = 2))
  }
  y <- indiv[,1]
  s <- exp(indiv[,2])
  dataset <- data.frame(cbind(y, s))
  return(list(data = dataset,
              para = c(n=n,mutheta=mutheta,musigma=musigma,
                       tautheta=tautheta, tausigma=tausigma, r1=r1, r2=r2, sigmas=sigmas)))
}


data <- data.nonreg(r1=0,r2=0)$data
res <- bimod(data$y, data$s, with_cov = FALSE, bi_option = TRUE,
             seed = 9999, iter = 3000, chain = 3, verbose = FALSE,
             warmup= 1000, thin=3, adapt_delta = 0.99)$pop_par
