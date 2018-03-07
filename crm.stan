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
  vector prior_lp(real r_global, vector r_local) {
    r_global ~ normal(0.0, 10.0);
    r_local ~ inv_chi_square(1.0);
        
    return r_global * sqrt_vec(r_local);
  }
  //Likelihood function for the cure model
  real ptcm_log(vector yobs, vector v, vector b, vector g, matrix X, real alpha, real mu, real a){
      real lpdf;
      real pdf;
      real cdf;
      real term1;
      real term2;
      vector[num_elements(yobs)] prob;
      vector[num_elements(yobs)] f;   //link function
      real lprob;
     f = exp(X * b) ./ (1 + a * exp(X * b) );
     
     for(i in 1:num_elements(yobs)){
        if( ( exp(-(mu + X[i,2:] * g) / alpha) <= 0) || (alpha <= 0) ){
                    prob[i] = 1e-10;  //catch errors for limit of computation
        }else{
          lpdf = weibull_lpdf(yobs[i] | alpha, exp(-(mu + X[i,2:] * g) / alpha) );
          pdf = exp(lpdf);
          term1 = (f[i] * pdf)^v[i];
          cdf = weibull_cdf(yobs[i] , alpha, exp(-(mu + X[i,2:] * g) / alpha) );
          term2 = (1 - ( a * f[i] * cdf))^(1 / (a - v[i]));
          prob[i] = term1 * term2;
        }
     }
     lprob = sum(log(prob));
     return lprob;
    }
}  

data {
  int<lower=0> N;
  int<lower=0> M_clinical;
  vector[N] yobs;   // observed time    
  vector[N] v_i;    //censor indicator
  matrix[N, M_clinical] X_clin;
}
  
transformed data {
  real<lower=0> tau_al;
  real<lower=0> tau_mu;

  tau_al = 10.0;
  tau_mu = 10.0;
  real<lower=0, upper=1> a;
    
  a = 0.7;  
}
  
parameters {
  real<lower=0> tau_s_cb_raw;
  vector<lower=0>[M_clinical] tau_cb_raw;
  vector[M_clinical] beta_clin_raw;
  
  real<lower=0> tau_s_cg_raw;
  vector<lower=0>[M_clinical-1] tau_cg_raw;
  vector[M_clinical - 1] gamma_clin_raw;
  
  real alpha_raw;
  real mu;
  
    
}
  
transformed parameters {
  vector[M_clinical] beta_clin;
  vector[M_clinical-1] gamma_clin;
  real<lower=0> alpha;

  beta_clin = prior_lp(tau_s_cb_raw, tau_cb_raw) .* beta_clin_raw;
  gamma_clin = prior_lp(tau_s_cg_raw, tau_cg_raw) .* gamma_clin_raw;
  
  alpha = exp(tau_al * alpha_raw);
    
}
  
model {
  yobs ~ ptcm(v_i, beta_clin, gamma_clin, X_clin, alpha, mu, a);
    
  beta_clin_raw ~ normal(0.0, 1.0);
  gamma_clin_raw ~ normal(0.0, 1.0);
  
  alpha_raw ~ normal(0.0, 1.0);
    
  mu ~ normal(0.0, tau_mu);
  
  a ~ beta(2,2);
    
}
  

  
  
