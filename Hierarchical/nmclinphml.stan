/*  Variable naming:
  obs       = observed
  cen       = (right) censored
  N         = number of samples
  M         = number of covariates
  clin        = established risk (or protective) factors
  gene      = candidate genearkers (candidate risk factors)
  tau       = scale parameter
  
  sequence_cum_sum : 
    position : position of the vector
  x : phi parameter (metastasis rate)
  te : theta parameter of proportional hazard
  
  mixcurelog : This function calculates the likelihood for the nonmixcure model:
    y obs time
  ai alpha weibull
  em mu weibull 
  pe phi rate of metastasis
  v censoring indicator
  theta proportional hazard prognostic index
  num index of individuals, (individuals must be ordered in descending order of observed time)
  Xclin must include the intercept
*/
functions {
  vector sqrt_vec(vector x) {
    vector[dims(x)[1]] res;

    for (m in 1:dims(x)[1]){
      res[m] = sqrt(x[m]);
    }

    return res;
  }
//Gaussian Prior
  vector g_prior_lp(real r_global, vector r_local) {
    r_global ~ normal(0.0, 10.0);
    r_local ~ inv_chi_square(1.0);

    return r_global * sqrt_vec(r_local);
  }
  //Likelihood function for the cure model
  real ptcm_log(vector yobs, vector v, vector lf, real alpha, vector sigma, int[] type){
      real lpdf;
      real pdf;
      real cdf;
      real term1;
      real term2;
      vector[num_elements(yobs)] prob;
      real lprob;
      
     for(i in 1:num_elements(yobs)){
          lpdf = weibull_lpdf(yobs[i] | alpha, sigma[type[i]]);
          pdf = exp(lpdf);
          term1 = - (log(lf[i]))*pdf;
          cdf = weibull_cdf(yobs[i], alpha, sigma[type[i]]);
          term2 = exp(log(lf[i]) * cdf);
          prob[i] = term1^v[i] * term2;
        }
     lprob = sum(log(prob));
     return lprob;
  }
}

data {
  int<lower=0> N;  //training N
  int<lower=0> J; //cluster of cancer
  int<lower=0> M; //number of clinical covariates
  vector[N] yobs;   // observed time (Training)
  vector[N] v;    //censor indicator (Training)
  matrix[N, M] Z_c; //matrix of covariate values
  int type[N];    //cancer type
  real varm[N]; //var of interest
  
}
  
transformed data {
  real<lower=0> tau_al;
  vector[N] icept;
  vector[N] varoi;
  real beta_unit;
  
  beta_unit = 5;

  tau_al = 10.0;

  for (i in 1:N){
     icept[i] = 1; 
  }
  varoi = to_vector(varm);
}
  
parameters {
  vector[J] beta0;
  real loc0; //b0 expectation
  real<lower=0> s0; //b0 scale
  
  vector[J] betaI;
  real locI; //b0 expectation
  real<lower=0> sI; //b0 scale
  
  vector[M] beta_c;
  

  real alpha_raw;
  
  real mu_mu;
  real mu_J;
  vector<lower=0>[J] s_mu;
  
  real<lower=0> tau_s_g_raw;
  vector<lower=0>[M] tau_g_raw;
  vector[M] gamma_g_raw;
  
  real gI_mu;
  real gI_J;
  vector<lower=0>[J] s_gI;
    
}
  
transformed parameters {
  real<lower=0> alpha;
  vector<lower=0>[N] sigma;
  vector<lower=0, upper=1>[N] lf;
  vector[J] mu;
  vector[J] gamma_I;
  vector[M] gamma_g;
  
  gamma_g = g_prior_lp(tau_s_g_raw, tau_g_raw) .* gamma_g_raw;
  for(i in 1:J){
    gamma_I[i] = gI_J + gI_mu * s_gI[i];
  }
    
  alpha = exp(tau_al * alpha_raw);
  for(i in 1:J){
    mu[i] = mu_J + mu_mu * s_mu[i];
  }
  for (i in 1:N){
      sigma[i] = exp(-(mu[type[i]] + Z_c[i] * gamma_g + varoi[i]*gamma_I[type[i]]) / alpha);
    }
    
  for (i in 1:N){
    lf[i]=1 ./ (1 +exp(-(icept[i]*beta0[type[i]]+ Z_c[i,]*beta_c + varoi[i]*betaI[type[i]] ))); //link function 
  }
}
  
model {
  beta0 ~ cauchy(loc0 , s0);
  loc0 ~ normal(0 ,1);
  s0 ~ cauchy(0, 2);
  
  betaI ~ cauchy(locI , sI);
  locI ~ normal(0 ,1);
  sI ~ exponential(1);

  beta_c ~ cauchy(0, 2.5);
  alpha_raw ~ normal(0.0, 1.0);
  
  mu_J ~ normal(0, 1);
  mu_mu ~ normal(0, 1);
  s_mu ~ exponential(1);
  
  gamma_g_raw ~ normal(0,1);
  
  gI_J ~ normal(0, 1);
  gI_mu ~ normal(0, 1);
  s_gI ~ exponential(1);
  
  yobs ~ ptcm(v, lf, alpha, sigma, type);

}
generated quantities{
    real lik[N];
    real lpdf;
    real pdf;
    real cdf;
    
    //loglik train
    for (i in 1:N) {
      //estimate log_lik
        lpdf = weibull_lpdf(yobs[i] | alpha, sigma[i] );
        pdf = exp(lpdf);
        cdf = weibull_cdf(yobs[i], alpha, sigma[i] );
        lik[i] = ((-(log(lf[i]))*pdf)^v[i]) * (exp(log(lf[i]) * cdf));
    }
    
}



          


  

  
  
