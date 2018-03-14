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
      if(x[m] >=0){
        res[m] = sqrt(x[m]);
      }else{
        print("Not a positive real number");
      }
    }
      return res;
  }
  //Horseshoe Prior
  vector prior_hs_lp(real r1_global, real r2_global, vector r1_local, vector r2_local, real nu_local, real scale_global, real nu_global){

    //half-t prior for lambda
    r1_local ~ normal(0.0, 1.0);
    r2_local ~ inv_gamma(0.5 * nu_local, 0.5 * nu_local);
    //half-t prior for tau
    r1_global ~ normal(0.0, scale_global);
    r2_global ~ inv_gamma(0.5 * nu_global, 0.5 * nu_global);

    return  r1_local .* sqrt(r2_local) * r1_global * sqrt(r2_global);
  }
  //Likelihood function for the cure model
  real ptcm_log(vector yobs_t, vector v_t, vector lf_t, real alpha, real sigma){
      real lpdf;
      real pdf;
      real cdf;
      real term1;
      real term2;
      vector[num_elements(yobs_t)] prob;
      real lprob;
      
     for(i in 1:num_elements(yobs_t)){
          lpdf = weibull_lpdf(yobs_t[i] | alpha, sigma);
          pdf = exp(lpdf);
          term1 = - (log(lf_t[i]))*pdf;
          cdf = weibull_cdf(yobs_t[i], alpha, sigma );
          term2 = exp(log(lf_t[i]) * cdf);
          prob[i] = term1^v_t[i] * term2;
        }
     lprob = sum(log(prob));
     return lprob;
  }
}

data {
  int<lower=0> N_t;  //training N
  int<lower=0> M_clinical;
  int<lower=0> M_g;
  vector[N_t] yobs_t;   // observed time (Training)
  vector[N_t] v_t;    //censor indicator (Training)
  matrix[N_t, M_clinical] Z_clin_t;
  vector[N_t] icept_t;
  matrix[N_t, M_g] Z_gene_t;
  
  int<lower=0> N_h;  //(Holdout)
  vector[N_h] yobs_h;   //(Holdout)
  matrix[N_h, M_clinical] Z_clin_h; //(Holdout)
  vector[N_h] v_h;    //censor indicator (Training)
  vector[N_h] icept_h;
  matrix[N_h, M_g] Z_gene_h;
  

  real<lower=1> nu_local; //degrees of freedom for lambda
  real<lower=1> nu_global; //degrees of freedom for tau 
  real<lower=0> scale_global;
}
  
transformed data {
  real<lower=0> tau_al;
  real<lower=0> tau_mu;

  tau_al = 10.0;
  tau_mu = 10.0;
}
  
parameters {
  vector[M_clinical] beta_clin;
  real beta0;

  real alpha_raw;
  real mu;
  
  vector[M_g] beta_g_raw;
  real<lower=0> tau1_global;
  real<lower=0> tau2_global;
  vector<lower=0>[M_g] tau1_local;
  vector<lower=0>[M_g] tau2_local;
    
}
  
transformed parameters {
  real<lower=0> alpha;
  real<lower=1e-10> sigma;
  vector<lower=1e-10>[N_t] lf_t;
  vector[M_g] beta_g;
  
  beta_g = prior_hs_lp(tau1_global, tau2_global, tau1_local, tau2_local, nu_local, scale_global, nu_global) .* beta_g_raw ;
  
  alpha = exp(tau_al * alpha_raw);
  sigma = exp(-(mu) / alpha);
  lf_t = 1 ./ (1 + exp( -(icept_t * beta0 + Z_clin_t * beta_clin + Z_gene_t * beta_g) )); //link function
  
}
  
model {
  
  beta_g_raw ~ normal(0.0, 1.0);
  
  beta0 ~ cauchy(0.0, 10.0);
  
  beta_clin ~ cauchy(0.0, 2.5);

  alpha_raw ~ normal(0.0, 1.0);
    
  mu ~ normal(0.0, tau_mu);
  
  yobs_t ~ ptcm(v_t, lf_t, alpha, sigma);

}
generated quantities{
    real lik_t[N_t];
    real lik_h[N_h];
    vector[N_h] lf_h;
    real lpdf;
    real pdf;
    real cdf;
    
    //loglik train
    for (i in 1:N_t) {
      //estimate log_lik
        lpdf = weibull_lpdf(yobs_t[i] | alpha, sigma );
        pdf = exp(lpdf);
        cdf = weibull_cdf(yobs_t[i], alpha, sigma );
        lik_t[i] = ((-(log(lf_t[i]))*pdf)^v_t[i]) * (exp(log(lf_t[i]) * cdf));
    }
    
  //calculate lf_h
  lf_h = 1 ./ (1 + exp( -(icept_h * beta0 +  Z_clin_h * beta_clin + Z_gene_h * beta_g) )); 
    //loglik holdout
    for (i in 1:N_h) {
      //estimate log_lik
        lpdf = weibull_lpdf(yobs_h[i] | alpha, sigma);
        pdf = exp(lpdf);
        cdf = weibull_cdf(yobs_h[i], alpha, sigma );
        lik_h[i] = ((-(log(lf_h[i]))*pdf)^v_h[i]) * (exp(log(lf_h[i]) * cdf));
    }
}


          


  

  
  
