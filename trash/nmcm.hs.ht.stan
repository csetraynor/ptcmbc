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
  real ptcm_log(vector yobs, vector v, vector f, real alpha, real mu){
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
              lpdf = weibull_lpdf(yobs[i] | alpha, exp(-(mu)/alpha));
              pdf = exp(lpdf);
              term1 = (-log(f[i])*pdf);
            }
          cdf = weibull_cdf(yobs[i] , alpha, exp(-(mu)/alpha));
          term2 = exp(log(f[i])*cdf);
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
  
  vector hs_prior_lp(real r1_global, real r2_global, vector r1_local, vector r2_local, real scale_global, real nu_global, real nu_local) {
    vector[num_elements(r1_local)] lambda;
    real tau;

    //half-t prior for tau
    r1_global ~ normal(0.0, scale_global);
    r2_global ~ inv_gamma(0.5 * nu_global, 0.5 * nu_global);

    //half-t prior for lambdas
    r1_local ~ normal(0.0, 1.0);
    r2_local ~ inv_gamma(0.5 * nu_local, 0.5 * nu_local);
    
    lambda = r1_local .* sqrt_vec(r2_local);
    tau = r1_global * sqrt(r2_global);

    return (lambda * tau);
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
  real<lower=1> nu_local;      //degrees of freedom for half-t tau
  real<lower=1> nu_global;     //degrees of freedom for half-t for lambdas
  real<lower=0> scale_global;  //scale for the half-t prior for tau
  vector[N] icept;
  
  for(i in 1:N){
    icept[i] = 1; 
  }
    
  tau_al = 10.0;
  
  nu_local = 1.0;    //horshor prior
  scale_global = 3e-3;
  nu_global = 1.0;  //half-cauchy prior  
    
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
  
  real beta0;
  real<lower=0> scale_icept;
    
  real alpha_raw;
    
  real<lower=0> tau_mu;
  real mu;
  
  real<lower=0> scale_global_raw;
  real<lower=1> nu_local_raw;
  real<lower=1> nu_global_raw;
    
}
  
transformed parameters {
  vector[M_clinical] beta_clin;
  vector[M_genomic] beta_gene;
  vector[N] f;                   //latent values
  real<lower=0> alpha;  

  beta_gene = hs_prior_lp(tau_s1_gene_raw, tau_s2_gene_raw, tau1_gene_raw, tau2_gene_raw, scale_global, nu_global, nu_local) .* beta_gene_raw;
  beta_clin = prior_lp(tau_s_clin_raw, tau_clin_raw) .* beta_clin_raw;
  
  alpha = exp(tau_al * alpha_raw);
  
  f = exp(-exp(icept * beta0 + X_clin * beta_clin + X_gene * beta_gene));
}
  
model {

  // gaussian prior for intercept
  scale_icept ~ exponential(0.1);
  beta0 ~ normal(0.0, scale_icept);
    
  beta_clin_raw ~ normal(0.0, 1.0);
  beta_gene_raw ~ normal(0.0, 1.0);
  
  alpha_raw ~ normal(0.0, 1.0);
    
  tau_mu ~ exponential(0.1);  
  mu ~ normal(0.0, tau_mu);
  
  //observation
  yobs ~ ptcm(v_i, f, alpha, mu);
    
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
            term1 = (-log(f[n])*pdf);
          }
        cdf = weibull_cdf(yobs[n] , alpha, exp(-(mu/alpha)));
        term2 = exp(log(f[n])*cdf);
        prob = term1 * term2;
      }else{
        prob = 1e-10;
        }
      log_lik[n] = log(prob);
    }
}
  

  
  
