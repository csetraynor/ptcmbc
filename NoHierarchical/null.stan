/*  Variable naming:
 obs       = observed
 cen       = (right) censored
 N         = number of samples
 tau       = scale parameter
*/
functions{
  real surv_log(vector yobs, vector v, real alpha, real sigma){
    real lpdf;
    real pdf;
    real cdf;
    vector[num_elements(yobs)] prob;
    real lprob;
    
    for (i in 1:num_elements(yobs)){
          lpdf = weibull_lpdf(yobs[i] | alpha, sigma);
          pdf = exp(lpdf);
          cdf = weibull_cdf(yobs[i], alpha, sigma );
          prob[i] = (pdf^v[i])*((1-cdf)^(1-v[i]));
    }
    lprob = sum(log(prob));
    return lprob;
  }
}



data {
  int<lower=0> N;
  vector[N] yobs;
  vector[N] v;
}

transformed data {
  real<lower=0> tau_mu;
  real<lower=0> tau_al;

  tau_mu = 10.0;
  tau_al = 10.0;
}

parameters {
  real alpha_raw;
  real mu;
}

transformed parameters {
  real alpha;
  real<lower=0> sigma;
  alpha = exp(tau_al * alpha_raw);
  sigma = exp(-(mu)/alpha);
}

model {

  alpha_raw ~ normal(0.0, 1.0);
  mu ~ normal(0.0, tau_mu);
  
  yobs ~ surv(v,   alpha, sigma);
}

generated quantities{
    real log_lik[N];
    real lpdf;
    real pdf;
    real cdf;

    for (i in 1:N) {
      
          lpdf = weibull_lpdf(yobs[i] | alpha, sigma);
          pdf = exp(lpdf);
          cdf = weibull_cdf(yobs[i], alpha, sigma );
          log_lik[i] =log((pdf^v[i])*((1-cdf)^(1-v[i]))) ;
        
    }
}
