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
  real ptcm_log(vector yobs, vector v, vector p, real alpha, real mu){
      real lpdf;
      real pdf;
      real cdf;
      real term1;
      real term2;
      vector[num_elements(yobs)] prob;
      real lprob;
        
        
     for(i in 1:num_elements(yobs)){
        if(alpha/exp(-(mu/alpha)) != 0 && yobs[i]/exp(-(mu/alpha)) != 0 && exp(-(mu/alpha)) != 0){ //catch errors for limit of computation
          if(v[i] == 0){ 
              term1 = 1;
            } else{
              lpdf = weibull_lpdf(yobs[i] | alpha, exp(-(mu/alpha)));
              pdf = exp(lpdf);
              term1 = (-log(p[i])*pdf);
            }
          cdf = weibull_cdf(yobs[i] , alpha, exp(-(mu/alpha)));
          term2 = exp(log(p[i])*cdf);
          prob[i] = term1 * term2;
        }else{
          prob[i] = 1e-10;
        }
      }
      lprob = sum(log(prob));
      return lprob;
    }
      
  vector sqrt_vec(vector x) {
      vector[dims(x)[1]] res;
        
    for (m in 1:dims(x)[1]){
        res[m] = sqrt(x[m]);
    }
        
      return res;
  }
      
  vector prior_lp(real r_global, vector r_local) {
    r_global ~ normal(0.0, 10.0);
    r_local ~ inv_chi_square(1.0);
        
    return r_global * sqrt_vec(r_local);
  }
      
}  

data {
  int<lower=0> N;
  int<lower=0> M_clinical;
  vector[N] yobs;   // observed time    
  vector[N] v_i;    //censor indicator
  matrix[N, M_clinical] X_clin;
  int<lower=1> num_id[N];   //identification of number
}
  
transformed data {
  real<lower=0> tau_al;
  real<lower=0> tau_mu;
    
  tau_al = 10.0;
  tau_mu = 10.0;
    
}
  
parameters {
  real<lower=0> tau_s_clin_raw;
  vector<lower=0>[M_clinical] tau_clin_raw;
    
  real alpha_raw;
  vector[M_clinical] beta_clin_raw;
    
  real mu_raw;
    
}
  
transformed parameters {
  vector[M_clinical] beta_clin;
  vector[M_clinical] p;
  real<lower=0> alpha;
  real<lower=0> sigma;
  real mu;
    
  mu = tau_mu * mu;
  beta_clin = prior_lp(tau_s_clin_raw, tau_clin_raw) .* beta_clin_raw;
  alpha = exp(tau_al * alpha_raw);
    
  sigma = exp(-(mu/alpha));
  p = exp(-exp(X_clin * beta_clin));
    
}
  
model {
  yobs ~ ptcm(v_i, p, alpha, sigma);
    
  beta_clin_raw ~ normal(0.0, 1.0);
  alpha_raw ~ normal(0.0, 1.0);
    
  mu_raw ~ normal(0.0, 1.0);
    
}
  

  
  
