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
  real ptcm_log(vector yobs, vector v, vector lf, real alpha, real sigma){
      real lpdf;
      real pdf;
      real cdf;
      real term1;
      real term2;
      vector[num_elements(yobs)] prob;
      real lprob;
      
     for(i in 1:num_elements(yobs)){
          lpdf = weibull_lpdf(yobs[i] | alpha, sigma);
          pdf = exp(lpdf);
          term1 = - (log(lf[i]))*pdf;
          cdf = weibull_cdf(yobs[i], alpha, sigma);
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
  
}
  
transformed data {
  real<lower=0> tau_al;
  vector[N] icept;
  real<lower=0> tau_mu;
  
  tau_mu = 10.0;
  tau_al = 10.0;

  for (i in 1:N){
     icept[i] = 1; 
  }

}
  
parameters {
  vector[J] beta0;
  real loc0; //b0 expectation
  real<lower=0> s0; //b0 scale
  
  vector[M] beta_c;

  real alpha_raw;
  real mu;
  
}
  
transformed parameters {
  real<lower=0> alpha;
  real<lower=0> sigma;
  vector<lower=0, upper=1>[N] lf;
    
  alpha = exp(tau_al * alpha_raw);
  sigma = exp(-(mu) / alpha);
  
  for (i in 1:N){
    lf[i]=1 ./ (1 +exp(-(icept[i]*beta0[type[i]]+ Z_c[i,]*beta_c))); //link function 
  }
}
  
model {
  beta0 ~ cauchy(loc0 , s0);
  loc0 ~ normal(0 ,1);
  s0 ~ exponential(.2);
  beta_c ~ student_t(5, 0, 2.5);

  alpha_raw ~ normal(0.0, 1.0);
  
  mu ~ normal(0.0, tau_mu);
  yobs ~ ptcm(v, lf, alpha, sigma);

}
generated quantities{
    real lik[N];
    real lpdf;
    real pdf;
    real cdf;
    
    //loglik train
    for (i in 1:N) {
      //estimate log_lik
        lpdf = weibull_lpdf(yobs[i] | alpha, sigma );
        pdf = exp(lpdf);
        cdf = weibull_cdf(yobs[i], alpha, sigma );
        lik[i] = ((-(log(lf[i]))*pdf)^v[i]) * (exp(log(lf[i]) * cdf));
    }
    
}


          


  

  
  
