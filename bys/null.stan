/*  Variable naming:
 obs       = observed
 cen       = (right) censored
 N         = number of samples
 tau       = scale parameter
*/
data {
  int<lower=0> Nobs;
  int<lower=0> Ncen;
  int<lower=0> J; //cohort
  vector[Nobs] yobs;
  vector[Ncen] ycen;
  vector[Nobs] Jobs;
  vector[Ncen] Jcen;
}

transformed data {
  real<lower=0> tau_mu;
  real<lower=0> tau_al;

  tau_mu = 10.0;
  tau_al = 10.0;
}

parameters {
  real alpha_raw;
  vector[J] theta;
  real<lower=0> tau_theta;
  real mu_theta;
}

transformed parameters {
  real alpha;
  alpha = exp(tau_al * alpha_raw);
}

model {
  for (i in 1:Nobs){
    yobs[i] ~ weibull(alpha, exp(-(theta[Jobs[i]])/alpha));
  }
  for (i in 1:Ncen){
    target += weibull_lccdf(ycen[i] | alpha, exp(-(theta[Jcen[i]])/alpha));
  }

  alpha_raw ~ normal(0.0, 1.0);
  theta ~ normal(loc_mu, tau_mu);
  mean_theta ~ normal(0.0, 1.0);
  tau_theta ~ cauchy(0 , 2);
}