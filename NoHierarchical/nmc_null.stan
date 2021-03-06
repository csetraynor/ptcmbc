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
  real ptcm_log(vector yobs, vector v, vector f, real alpha, real mu){
      real lpdf;
      real pdf;
      real cdf;
      real term1;
      real term2;
      vector[num_elements(yobs)] prob;
      real lprob;
     
     for(i in 1:num_elements(yobs)){
        if( exp(-(mu) / alpha) <= 0 || exp(-(mu) / alpha) > 10e10 || alpha <= 0 || alpha > 10e10){
                    prob[i] = 1e-10;  //catch errors
        }else{
          lpdf = weibull_lpdf(yobs[i] | alpha, exp(-(mu) / alpha) );
          pdf = exp(lpdf);
          term1 = - (log(f[i]))*pdf;
          cdf = weibull_cdf(yobs[i], alpha, exp(-(mu) / alpha) );
          term2 = exp(log(f[i]) * cdf);
          prob[i] = term1^v[i] * term2;
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
  matrix[N, M_clinical] Z_clin;
}
  
transformed data {
  real<lower=0> tau_al;
  real<lower=0> tau_mu;
  
  tau_al = 10.0;
  tau_mu = 10.0;

}
  
parameters {

  real alpha_raw;
  real mu;
  
  vector<lower=0, upper= 1>[N] p;
    
}
  
transformed parameters {
  real<lower=0> alpha;
  
  
  alpha = exp(tau_al * alpha_raw);
  

}
  
model {
  p ~ beta(1.0, 1.0);
  alpha_raw ~ normal(0.0, 1.0);
    
  mu ~ normal(0.0, tau_mu);
  
  yobs ~ ptcm(v_i, p, alpha, mu);
  
}


  

  
  
