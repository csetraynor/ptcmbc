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
        if(exp(-(mu/alpha)) > 0){ //catch errors for limit of computation
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
  
  vector hs_prior_lp(real r1_global, real r2_global, vector r1_local, vector r2_local, real nu) {
    r1_global ~ normal(0.0, 1.0);
    r2_global ~ inv_gamma(0.5, 0.5);

    r1_local ~ normal(0.0, 1.0);
    r2_local ~ inv_gamma(0.5 * nu, 0.5 * nu);

    return (r1_global * sqrt(r2_global)) * r1_local .* sqrt_vec(r2_local);
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
  int M_genomic;
  vector[N] yobs;   // observed time    
  vector[N] v_i;    //censor indicator
  matrix[N, M_clinical] X_clin;
  matrix[N, M_genomic] X_gene;

}
  
transformed data {
  real<lower=0> tau_al;
  real<lower=1> nu;
    
  tau_al = 10.0;
  nu = 3.0;
    
}
  
parameters {
  vector[M_clinical] beta_clin_raw;
  real<lower=0> tau_s_clin_raw;
  vector<lower=0>[M_clinical] tau_clin_raw;
  
  vector[M_genomic] beta_gene_raw;
  real<lower=0> tau_s1_gene_raw;
  real<lower=0> tau_s2_gene_raw;
  vector<lower=0>[M_genomic] tau1_gene_raw;
  vector<lower=0>[M_genomic] tau2_gene_raw;
    
  real alpha_raw;
    
  real<lower=0> tau_mu;
  real mu;
    
}
  
transformed parameters {
  vector[M_clinical] beta_clin;
  vector[M_genomic] beta_gene;
  vector[N] p;
  real<lower=0> alpha;
    
  beta_gene = hs_prior_lp(tau_s1_gene_raw, tau_s2_gene_raw, tau1_gene_raw, tau2_gene_raw, nu) .* beta_gene_raw;
  beta_clin = prior_lp(tau_s_clin_raw, tau_clin_raw) .* beta_clin_raw;
  
  alpha = exp(tau_al * alpha_raw);
  
  p = exp(-exp(X_clin * beta_clin + X_gene * beta_gene));
    
}
  
model {
  yobs ~ ptcm(v_i, p, alpha, mu);
    
  beta_clin_raw ~ normal(0.0, 1.0);
  beta_gene_raw ~ normal(0.0, 1.0);
  alpha_raw ~ normal(0.0, 1.0);
    
  mu ~ normal(0.0, tau_mu);
  tau_mu ~ cauchy(0.0, 2.0);
    
}

generated quantities{
  vector[N] log_lik;
  real term1;
  real term2;
  real cdf;
  real pdf;
  real lpdf;
  real prob;

  for(n in 1:N){
      if(exp(-(mu/alpha)) > 0){ //catch errors for limit of computation
        if(v_i[n] == 0){ 
            term1 = 1;
          } else{
            lpdf = weibull_lpdf(yobs[n] | alpha, exp(-(mu/alpha)));
            pdf = exp(lpdf);
            term1 = (-log(p[n])*pdf);
          }
        cdf = weibull_cdf(yobs[n] , alpha, exp(-(mu/alpha)));
        term2 = exp(log(p[n])*cdf);
        prob = term1 * term2;
      }else{
        prob = 1e-10;
        }
      log_lik[n] = log(prob);
    }
}
  

  
  
