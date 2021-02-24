#' Univariate Bayestian hierarchical model for combining summaries from multiple sources.
#' @description \code{unimod} is mainly modeling the measures for combining summaries from multiple sources. It returns the estimates and posterior distributions of parameters of interest fitted with data.
#'
#' @param y A n-by-1 vector, consisting of the measures of n observations.
#' @param s A n-by-1 vector, consisting of the standard error of the measures for n observations.
#' @param with_cov if TRUE, include covariates in the model.
#' @param X A n-by-p matrix of p covariates including intercept, which is all 1 for the first column.
#' @param seed The seed for random number generation. The default is generated from 1 to the maximum integer supported by R on the machine. Even if multiple chains are used, only one seed is needed, with other chains having seeds derived from that of the first chain to avoid dependent samples.The default is 9999.
#' @param iter a number to specify the iteration times of Stan. The default is 3000.
#' @param chains A positive integer specifying the number of Markov chains. The default is 4.
#' @param warmup a number of warmup iterations to be excluded when computing the summaries. The default is 1000.
#' @param thin A positive integer specifying the period for saving samples. The default value is 3.
#' @param adapt_delta A parameters to control the sampler’s behavior. The default is 0.8.
#' @param max_treedepth A parameters to control the sampler’s behavior in For algorithm NUTS.
#' @param verbose TRUE or FALSE: flag indicating whether to print intermediate output from Stan on the console, which might be helpful for model debugging.
#' @return a list of information of the fitted univariate Bayestian hierarchical model including
#' \describe{
#'   \item{\code{pop_par}}{the estimates and credible intervals for the population-level parameters}
#'   \item{\code{ind_par}}{the estimates and credible intervals for the individual-level parameters}
#'   \item{\code{poster_dist}}{posterior distributions of all parameters fitted in the model}
#' }
#' @example man/unimod_example.R
#' @import rstan
#' @import MASS
#' @export


unimod <- function(y = data$y, s = data$s, with_cov = FALSE, X = NULL,
                   seed = 9999, iter = 3000, chain = 4, verbose = FALSE, warmup= 1000, thin=3,
                   adapt_delta = 0.99, max_treedepth =3, ...){
# 0. check and warnings
  if(is.vector(y)==FALSE) {stop("y has to be a vector.")}
  if(is.vector(s)==FALSE) {stop("s has to be a vector.")}
  if(with_cov == TRUE & is.null(X) == 1) {stop("covariate X has to be included")}
  if(with_cov == FALSE & is.null(X) == 0) {stop("covariate X has to be null")}
# 1. first define the stan models
# 1.1. univariate without covariates
  uni.stan =
  "
  data {
    int<lower=0> N;
    real y[N];
    real<lower=0> sigma[N];
  }
  parameters {
      real mu;
      real<lower=0> tau;
      vector[N] eta;
  }
  transformed parameters {
      vector[N] theta;
      theta = mu + eta*tau;
  }
  model {
      tau ~ cauchy(0,2.5);
      eta ~ normal_lpdf(0,1);
      y ~ normal_lpdf(theta, sigma);
  }

  "
# 1.2. univariate with covariates
  uni.wcov.stan =
  "
  data {
      int<lower=0> N;
      int<lower=1> K;
      real y[N];
      real<lower=0> sigma[N];
      matrix[N,K] X;
  }
  parameters {
      vector[K] beta;
      real<lower=0> tau;
      vector[N] eta;
  }
  transformed parameters {
      vector[N] theta;
      theta = X*beta + eta*tau;
  }
  model {
      tau ~ cauchy(0,2.5);
      eta  ~ normal(0,1);
      to_vector(beta) ~ normal(0, 1e6);
      y ~ normal(theta, sigma);
  }

  "
  # 2. then prepare the dataset, and run the stan model
  if(with_cov == F){
    dat = list(N=length(y),y=y,sigma=s)
    fit <- stan(model_code = uni.stan,
                data = dat,seed = seed, iter = iter, chains = chain, warmup= warmup,
                thin = thin, verbose = TRUE,control = list(adapt_delta = adapt_delta, max_treedepth =  max_treedepth))
    fit.result <- data.frame(as.matrix(fit))
    # result
    n <- length(y)
    uni_mu_dist <- fit.result$mu
    uni_theta_dist <- matrix(rep(NA),nrow=n,ncol=dim(fit.result)[1])
    for (i in 1:n) uni_theta_dist[i,] <- fit.result[,1+1+n+i]
    pop_par <- c(mean = mean(uni_mu_dist),
                 lower.bound = quantile(uni_mu_dist, 0.05),
                 upper.bound = quantile(uni_mu_dist, 0.95))
    ind_par <- do.call(rbind, lapply(1:n, function(x)
      c(mean = mean(uni_theta_dist[x,]),
        lower.bound = quantile(uni_theta_dist[x,], 0.05),
        upper.bound = quantile(uni_theta_dist[x,], 0.95))))
  }
  if(with_cov == T){
    dat=list(N=length(y),y=y,sigma=s,X=X, K=ncol(X))
    fit <- stan(model_code = uni.wcov.stan,
                data = dat,seed = seed, iter = iter, chains = chain, warmup= warmup,
                thin = thin, verbose = verbose,control = list(adapt_delta = adapt_delta))
    fit.result <- data.frame(as.matrix(fit))
    n <- length(y)
    p <- ncol(X)
    uni_beta_dist <- matrix(rep(NA),nrow=p,ncol=dim(fit.result)[1])
    uni_theta_dist <- matrix(rep(NA),nrow=n,ncol=dim(fit.result)[1])
    for (i in 1:p) uni_beta_dist[i,] <- fit.result[,i]
    for (i in 1:n) uni_theta_dist[i,] <- fit.result[,p+1+n+i]
    pop_par <- do.call(rbind, lapply(1:p, function(x)
      c(mean = mean(uni_beta_dist[x,]),
        lower.bound = quantile(uni_beta_dist[x,], 0.05),
        upper.bound = quantile(uni_beta_dist[x,], 0.95))))
    ind_par <- do.call(rbind, lapply(1:n, function(x)
      c(mean = mean(uni_theta_dist[x,]),
        lower.bound = quantile(uni_theta_dist[x,], 0.05),
        upper.bound = quantile(uni_theta_dist[x,], 0.95))))
  }
  return(list(pop_par = pop_par,
              ind_par = ind_par,
              poster_dist = fit.result))
}

#' Bivariate Bayestian hierarchical model for combining summaries from multiple sources.
#' @description \code{bimod} is modeling the measures and its uncertainty jointly for combining summaries from multiple sources. It returns the estimates and posterior distributions of parameters of interest fitted with data.
#' @inheritParams unimod
#' @param bi_option if TRUE, use empirical value for standard deviation of uncertainty of measure
#' @return a list of information of the fitted bivariate Bayestian hierarchical model including
#' \describe{
#'   \item{\code{pop_par}}{the estimates and credible intervals for the population-level parameters}
#'   \item{\code{ind_par}}{the estimates and credible intervals for the individual-level parameters}
#'   \item{\code{poster_dist}}{posterior distributions of all parameters fitted in the model}
#' }
#' @example man/bimod_example.R
#' @export


bimod <- function(y = data$y, s = data$s, with_cov = FALSE, X = NULL, bi_option = TRUE,
                  seed = 9999, iter = 3000, chain = 4, verbose = FALSE, warmup= 1000, thin=3,
                  adapt_delta = 0.99, max_treedepth =  3, ...){
  # 0. check and warnings
  if(is.vector(y)==FALSE) {stop("y has to be a vector.")}
  if(is.vector(s)==FALSE) {stop("s has to be a vector.")}
  if(with_cov == TRUE & is.null(X) == 1) {stop("covariate X has to be included")}
  if(with_cov == FALSE & is.null(X) == 0) {stop("covariate X has to be null")}
  # 1. first define the stan models
  # 1.1 bivariate without covariates
  bi.stan =
    "
  data {
        int<lower=1> N;
        vector[2] x[N];
  }
  parameters {
        vector[2] mu;
        vector[2] theta[N];
        real<lower=-1,upper=1> r1;
        vector<lower=0>[2] st_devs;
        cholesky_factor_corr[2] L_corr;
        real<lower=0> sigma_s[N];
  }
  transformed parameters {
        real<lower=0> sigma[N];
        matrix[2,2] Sigma[N];
        for (i in 1:N){
        sigma[i] = exp(theta[i,2]);
        Sigma[i][1,1] = square(sigma[i]);
        Sigma[i][2,2] = square(sigma_s[i]);
        Sigma[i][1,2] = r1*sigma[i]*sigma_s[i];
        Sigma[i][2,1] = Sigma[i][1,2];
        }
  }
  model {
        st_devs ~ cauchy(0,2.5);
        L_corr ~ lkj_corr_cholesky(4);
        theta ~ multi_normal_cholesky(mu,diag_pre_multiply(st_devs, L_corr));
        for (i in 1:N){
          sigma_s[i]~ cauchy(0,2.5);
          x[i] ~ multi_normal(theta[i], Sigma[i]);
        }
  }

  "
  # 1.2 bivariate with covariates
  bi.wcov.stan =
    "
    data {
    int<lower=1> N;
    int<lower=1> K;
    vector[2] x[N];
    vector[K] X[N];
    }
  parameters {
    matrix[2,K] beta;
    vector[2] theta[N];
    real<lower=-1,upper=1> r1;
    vector<lower=0>[2] st_devs;
    cholesky_factor_corr[2] L_corr;
    real<lower=0> sigma_s[N];
    }
  transformed parameters {
    vector[2] mu[N];
    real<lower=0> sigma[N];
    matrix[2,2] Sigma[N];
    for (i in 1:N){
      mu[i] = beta * X[i];
      sigma[i] = exp(theta[i,2]);
      Sigma[i][1,1] = square(sigma[i]);
      Sigma[i][2,2] = square(sigma_s[i]);
      Sigma[i][1,2] = r1*sigma[i]*sigma_s[i];
      Sigma[i][2,1] = Sigma[i][1,2];
      }
    }
  model {
    to_vector(beta) ~ normal(0, 1e3);
    st_devs ~ cauchy(0,2.5);
    L_corr ~ lkj_corr_cholesky(4);
    theta ~ multi_normal_cholesky(mu,diag_pre_multiply(st_devs, L_corr));
    for (i in 1:N) {
      sigma_s[i]~ cauchy(0,2.5);
      x[i] ~ multi_normal(theta[i], Sigma[i]);
    }
  }

  "
  # 1.3 bivariate without covariates, fixed s_i
  bi.f.stan =
    "
  data {
        int<lower=1> N;
        vector[2] x[N];
  }
  parameters {
        vector[2] mu;
        vector[2] theta[N];
        real<lower=-1,upper=1> r1;
        vector<lower=0>[2] st_devs;
        cholesky_factor_corr[2] L_corr;
  }
  transformed parameters {
        real<lower=0> sigma[N];
        matrix[2,2] Sigma[N];
        for (i in 1:N){
        sigma[i] = exp(theta[i,2]);
        Sigma[i][1,1] = square(sigma[i]);
        Sigma[i][1,2] = r1*sigma[i];
        Sigma[i][2,1] = Sigma[i][1,2];
        Sigma[i][2,2] = 1;
        }
  }
  model {
        st_devs ~ cauchy(0,2.5);
        L_corr ~ lkj_corr_cholesky(4);
        theta ~ multi_normal_cholesky(mu,diag_pre_multiply(st_devs, L_corr));
        for (i in 1:N) x[i] ~ multi_normal(theta[i], Sigma[i]);
  }

  "

  # 1.4 bivariate with covariates, fixed s_i
  bi.wcov.f.stan =
    "
  data {
    int<lower=1> N;
    int<lower=1> K;
    vector[2] x[N];
    vector[K] X[N];
    }
  parameters {
    matrix[2,K] beta;
    vector[2] theta[N];
    real<lower=-1,upper=1> r1;
    vector<lower=0>[2] st_devs;
    cholesky_factor_corr[2] L_corr;
    }
  transformed parameters {
    vector[2] mu[N];
    real<lower=0> sigma[N];
    matrix[2,2] Sigma[N];
    for (i in 1:N){
      mu[i] = beta * X[i];
      sigma[i] = exp(theta[i,2]);
      Sigma[i][1,1] = square(sigma[i]);
      Sigma[i][1,2] = r1*sigma[i];
      Sigma[i][2,1] = Sigma[i][1,2];
      Sigma[i][2,2] = 1;
      }
    }
  model {
    to_vector(beta) ~ normal(0, 1e3);
    st_devs ~ cauchy(0,2.5);
    L_corr ~ lkj_corr_cholesky(4);
    theta ~ multi_normal_cholesky(mu,diag_pre_multiply(st_devs, L_corr));
    for (i in 1:N)
      x[i] ~ multi_normal(theta[i], Sigma[i]);
  }

  "
  # 2. then prepare the dataset, and run the stan model
  if(with_cov == F & bi_option == T){
    dat=list(N=length(y),x=cbind(y,log(s)))
    fit <- stan(model_code = bi.stan,
                data = dat,seed = seed, iter = iter, chains = chain, warmup= warmup,
                thin = thin, verbose = TRUE,
                control = list(adapt_delta = adapt_delta, max_treedepth =  max_treedepth))
    fit.result <- data.frame(as.matrix(fit))
    # result
    n <- length(y)
    bi_mu_dist <- fit.result$mu.1.
    bi_theta_dist <- matrix(rep(NA),nrow=n,ncol=dim(fit.result)[1])
    for (i in 1:n) bi_theta_dist[i,] <- fit.result[,2+i]
    pop_par <- c(mean = mean(bi_mu_dist),
                 lower.bound = quantile(bi_mu_dist, 0.05),
                 upper.bound = quantile(bi_mu_dist, 0.95))
    ind_par <- do.call(rbind, lapply(1:n, function(x)
      c(mean = mean(bi_theta_dist[x,]),
        lower.bound = quantile(bi_theta_dist[x,], 0.05),
        upper.bound = quantile(bi_theta_dist[x,], 0.95))))
  }
  if(with_cov == T & bi_option == T){
    dat=list(N=length(y),x=cbind(y,log(s)),X=X, K=ncol(X))
    fit <- stan(model_code = bi.wcov.stan,
                data = dat,seed = seed, iter = iter, chains = chain, warmup= warmup,
                thin = thin, verbose = TRUE,
                control = list(adapt_delta = adapt_delta, max_treedepth =  max_treedepth))
    fit.result <- data.frame(as.matrix(fit))
    # result
    n <- length(y)
    p <- ncol(X)
    bi_beta_dist <- matrix(rep(NA),nrow=p,ncol=dim(fit.result)[1])
    bi_theta_dist <- matrix(rep(NA),nrow=n,ncol=dim(fit.result)[1])
    for (i in 1:p) bi_beta_dist[i,] <- fit.result[,i]
    for (i in 1:n) bi_theta_dist[i,] <- fit.result[,2*p+i]
    pop_par <- do.call(rbind, lapply(1:p, function(x)
      c(mean = mean(bi_beta_dist[x,]),
        lower.bound = quantile(bi_beta_dist[x,], 0.05),
        upper.bound = quantile(bi_beta_dist[x,], 0.95))))
    ind_par <- do.call(rbind, lapply(1:n, function(x)
      c(mean = mean(bi_theta_dist[x,]),
        lower.bound = quantile(bi_theta_dist[x,], 0.05),
        upper.bound = quantile(bi_theta_dist[x,], 0.95))))
  }
  if(with_cov == F & bi_option == F){
    dat=list(N=length(y),x=cbind(y,log(s)))
    fit <- stan(model_code = bi.f.stan,
                data = dat,seed = seed, iter = iter, chains = chain, warmup= warmup,
                thin = thin, verbose = TRUE,
                control = list(adapt_delta = adapt_delta, max_treedepth =  max_treedepth))
    fit.result <- data.frame(as.matrix(fit))
    # result
    n <- length(y)
    bi_mu_dist <- bi.result$mu.1.
    bi_theta_dist <- matrix(rep(NA),nrow=n,ncol=dim(fit.result)[1])
    for (i in 1:n) bi_theta_dist[i,] <- fit.result[,2+i]
    pop_par <- c(mean = mean(bi_mu_dist),
                 lower.bound = quantile(bi_mu_dist, 0.05),
                 upper.bound = quantile(bi_mu_dist, 0.95))
    ind_par <- do.call(rbind, lapply(1:n, function(x)
      c(mean = mean(bi_theta_dist[x,]),
        lower.bound = quantile(bi_theta_dist[x,], 0.05),
        upper.bound = quantile(bi_theta_dist[x,], 0.95))))
  }
  if(with_cov == T & bi_option == F){
    dat=list(N=length(y),x=cbind(y,log(s)),X=X, K=ncol(X))
    fit <- stan(model_code = bi.wcov.f.stan,
                data = dat,seed = seed, iter = iter, chains = chain, warmup= warmup,
                thin = thin, verbose = TRUE,
                control = list(adapt_delta = adapt_delta, max_treedepth =  max_treedepth))
    fit.result <- data.frame(as.matrix(fit))
    # result
    n <- length(y)
    p <- ncol(X)
    bi_beta_dist <- matrix(rep(NA),nrow=p,ncol=dim(fit.result)[1])
    bi_theta_dist <- matrix(rep(NA),nrow=n,ncol=dim(fit.result)[1])
    for (i in 1:p) bi_beta_dist[i,] <- fit.result[,i]
    for (i in 1:n) bi_theta_dist[i,] <- fit.result[,2*p+i]
    pop_par <- do.call(rbind, lapply(1:p, function(x)
      c(mean = mean(bi_beta_dist[x,]),
        lower.bound = quantile(bi_beta_dist[x,], 0.05),
        upper.bound = quantile(bi_beta_dist[x,], 0.95))))
    ind_par <- do.call(rbind, lapply(1:n, function(x)
      c(mean = mean(bi_theta_dist[x,]),
        lower.bound = quantile(bi_theta_dist[x,], 0.05),
        upper.bound = quantile(bi_theta_dist[x,], 0.95))))
  }
  return(list(pop_par = pop_par,
              ind_par = ind_par,
              poster_dist = fit.result))
}


