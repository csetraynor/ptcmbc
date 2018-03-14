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
  vector[N_t] yobs_t;   // observed time (Training)
  vector[N_t] v_t;    //censor indicator (Training)
  matrix[N_t, M_clinical] Z_clin_t;
  vector[N_t] icept_t;
  vector[N] type;  //cancer type
  
  int<lower=0> N_h;  //(Holdout)
  vector[N_h] yobs_h;   //(Holdout)
  matrix[N_h, M_clinical] Z_clin_h; //(Holdout)
  vector[N_h] v_h;    //censor indicator (Training)
  vector[N_h] icept_h;
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
    
}
  
transformed parameters {
  real<lower=0> alpha;
  real<lower=1e-10> sigma;
  vector<lower=1e-10>[N_t] lf_t;
  
  alpha = exp(tau_al * alpha_raw);
  sigma = exp(-(mu) / alpha);
  lf_t = 1 ./ (1 + exp( -(icept_t * beta0 + Z_clin_t * beta_clin) )); //link function
  
}
  
model {
 
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
    lf_h = 1 ./ (1 + exp( -(icept_h * beta0 +  Z_clin_h * beta_clin) )); 
    //loglik holdout
    for (i in 1:N_h) {
      //estimate log_lik
        lpdf = weibull_lpdf(yobs_h[i] | alpha, sigma);
        pdf = exp(lpdf);
        cdf = weibull_cdf(yobs_h[i], alpha, sigma );
        lik_h[i] = ((-(log(lf_h[i]))*pdf)^v_h[i]) * (exp(log(lf_h[i]) * cdf));
    }
}


          


  

  
  
