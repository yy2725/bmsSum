#' Bayesian hierarchical models for the combining the summay measures
#' \code{sumbms} returns the posterior distributions of parameters fitted with data \code{data}
#' @param data A data frame.
#' @param dim A number 1 or 2.
#' @param with_cov logical, TRUE or FALSE
#' @param bi_option logical, TRUE or FALSE
#' @return The posterior distributions of parameters fitted with data \code{data} with univariate or bivariate model.
#' @examples
#' add(data, 1)
#' add(data, 2)

bmsSum::unimod <- function(data, with_cov = FALSE, bi_option = TRUE,...){
# first define the stan models
# (1) univariate without covariates
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
# (2) univariate with covariates
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

  # 1
  if(dim==1 & with_cov == F){
    dat=list(N=dim(data)[1],y=data$y,sigma=sqrt(data$s2))
    fit <- stan(model_code = uni.stan,
                data = dat,seed = 9999, iter = 5000, chains = 3, warmup= 2000,
                thin=10,verbose = TRUE,control = list(adapt_delta = 0.99))
    fit.result <- data.frame(as.matrix(fit))
  }
  # 2
  if(dim==1 & with_cov == T){
    X = data[,5:dim(data)[2]]
    dat=list(N=dim(data)[1],y=data$y,sigma=sqrt(data$s2),X=X, K=ncol(X))
    fit <- stan(model_code = uni.wcov.stan,
                data = dat, seed = 9999, iter = 5000, chains = 3, warmup= 2000,
                thin=10, verbose = TRUE, control = list(adapt_delta = 0.99))
    fit.result <- data.frame(as.matrix(fit))
    }

  return(fit.result)
}

bmsSum::bimod <- function(data, with_cov = FALSE, bi_option = TRUE,...){
  # first define the stan models
  # (3) bivariate without covariates
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
  # (4) bivariate with covariates
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
    to_vector(beta) ~ normal(0, 1e6);
    st_devs ~ cauchy(0,2.5);
    L_corr ~ lkj_corr_cholesky(4);
    theta ~ multi_normal_cholesky(mu,diag_pre_multiply(st_devs, L_corr));
    for (i in 1:N) {
      sigma_s[i]~ cauchy(0,2.5);
      x[i] ~ multi_normal(theta[i], Sigma[i]);
    }
  }

  "
  # (5) bivariate without covariates, fixed s_i
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

  # (6) bivariate with covariates, fixed s_i
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
    to_vector(beta) ~ normal(0, 1e6);
    st_devs ~ cauchy(0,2.5);
    L_corr ~ lkj_corr_cholesky(4);
    theta ~ multi_normal_cholesky(mu,diag_pre_multiply(st_devs, L_corr));
    for (i in 1:N)
      x[i] ~ multi_normal(theta[i], Sigma[i]);
  }

  "
  # 3
  if(dim==2 & with_cov == F & bi_option == T){
    dat=list(N=dim(data)[1],x=cbind(data$y,log(sqrt(data$s2))))
    fit <- stan(model_code = bi.stan,
                data = dat,seed = 9999, iter = 5000, chains = 3, warmup= 2000,
                thin=10, verbose = TRUE, control = list(adapt_delta = 0.9999, max_treedepth=15))
    fit.result <- data.frame(as.matrix(fit))
  }
  # 4
  if(dim==2 & with_cov == T & bi_option == T){
    X = data[,5:dim(data)[2]]
    dat=list(N=dim(data)[1],x=cbind(data$y,log(sqrt(data$s2))),X=X, K=ncol(X))
    fit <- stan(model_code = bi.wcov.stan,
                data = dat,seed = 9999, iter = 5000, chains = 3, warmup= 2000,
                thin=10, verbose = TRUE,control = list(adapt_delta = 0.9999, max_treedepth=15))
    fit.result <- data.frame(as.matrix(fit))
  }
  # 5
  if(dim==2 & with_cov == F & bi_option == F){
    X = data[,5:dim(data)[2]]
    dat=list(N=dim(data)[1],x=cbind(data$y,log(sqrt(data$s2))),X=X, K=ncol(X))
    fit <- stan(model_code = bi.f.stan,
                data = dat,seed = 9999, iter = 5000, chains = 3, warmup= 2000,
                thin=10, verbose = TRUE,control = list(adapt_delta = 0.9999, max_treedepth=15))
    fit.result <- data.frame(as.matrix(fit))
  }
  # 6
  if(dim==2 & with_cov == T & bi_option == F){
    X = data[,5:dim(data)[2]]
    dat=list(N=dim(data)[1],x=cbind(data$y,log(sqrt(data$s2))),X=X, K=ncol(X))
    fit <- stan(model_code = bi.wcov.f.stan,
                data = dat,seed = 9999, iter = 5000, chains = 3, warmup= 2000,
                thin=10, verbose = TRUE,control = list(adapt_delta = 0.9999, max_treedepth=15))
    fit.result <- data.frame(as.matrix(fit))
  }
  return(fit.result)
}

# test
# data <- data.nonreg(r1=0,r2=0)$data
# res <- sumbms(data, dim=1, with_cov = FALSE, bi_option = TRUE)
