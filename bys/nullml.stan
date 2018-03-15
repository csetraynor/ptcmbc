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
  int Jobs[Nobs];
  int Jcen[Ncen];
}

transformed data {
  real<lower=0> tau_al;
  tau_al = 10.0;
}

parameters {
  real alpha_raw;
  real<lower=0> tau_theta;
  vector[J] theta;
  real mu_theta;
}

transformed parameters {
  real<lower=0> alpha;
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
  theta ~ normal(mu_theta, tau_theta);
  mu_theta ~ normal(0 ,1);
  tau_theta ~ exponential(.1);
}
