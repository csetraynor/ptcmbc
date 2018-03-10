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
  //Horseshoe Prior
  vector prior_hs_lp(real r1_global, real r2_global, vector r1_local, vector r2_local, real nu_local, real scale_global, real nu_global){
    vector[num_elements(r1_local)] lambda;
    real tau;
    
    //half-t prior for lambda
    r1_local ~ normal(0.0, 1.0);
    r2_local ~ inv_gamma(.5 * nu_local, .5 * nu_local);
    //half-t prior for tau
    r1_global ~ normal(0.0, scale_global);
    r2_global ~ inv_gamma(.5 * nu_global, .5 * nu_global);
    
    lambda =  r1_local .* sqrt_vec(r2_local);
    tau = r1_global * sqrt(r2_global);
    
    return lambda *  tau;
  }
  //Gaussian Prior
  vector prior_lp(real r_global, vector r_local) {
    r_global ~ normal(0.0, 10.0);
    r_local ~ inv_chi_square(1.0);
        
    return r_global * sqrt_vec(r_local);
  }
  //Likelihood function for the cure model
  real ptcm_log(vector yobs, vector v, vector beta_c, vector beta_g, matrix Z_clin, matrix Z_gene, real alpha, real mu){
      real lpdf;
      real pdf;
      real cdf;
      real term1;
      real term2;
      real prob[num_elements(yobs)];
      vector[num_elements(yobs)] f;
      real lprob;
      
      f = exp( Z_clin * beta_c + Z_gene * beta_g) ./ (1 + exp(Z_clin * beta_c + Z_gene * beta_g)); //link function
     
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
  int<lower=0> M_c; //clinical covariates
  int<lower=0> M_g; //genomic covariates
  vector[N] yobs;   // observed time    
  vector[N] v;    //censor indicator
  matrix[N, M_c] Z_clin;
  matrix[N, M_g] Z_gene;
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
  real<lower=0> tau_s_cb_raw;
  vector<lower=0>[M_c] tau_cb_raw;
  vector[M_c] beta_c_raw;

  real alpha_raw;
  real mu;
  
  vector[M_g] beta_g_raw;
  real<lower=0> tau1_global;
  real<lower=0> tau2_global;
  vector<lower=0>[M_g] tau1_local;
  vector<lower=0>[M_g] tau2_local;
}
  
transformed parameters {
  vector[M_c] beta_c;
  vector[M_g] beta_g;
  real<lower=0> alpha;
  
  beta_g = prior_hs_lp(tau1_global, tau2_global, tau1_local, tau2_local, nu_local, scale_global, nu_global) .* beta_g_raw ;
  
  beta_c = prior_lp(tau_s_cb_raw, tau_cb_raw) .* beta_c_raw;
  
  alpha = exp(tau_al * alpha_raw);
}
  
model {
  beta_g_raw ~ normal(0.0, 1.0);
  
  beta_c_raw ~ normal(0.0, 1.0);

  alpha_raw ~ normal(0.0, 1.0);
    
  mu ~ normal(0.0, tau_mu);
  
  yobs ~ ptcm(v, beta_c, beta_g, Z_clin,  Z_gene, alpha, mu);

}

generated quantities{

    real log_lik[N];
    vector[N] lp;
    real lpdf;
    real pdf;
    real cdf;
    real term1;
    real term2;

    //calculate lp
    lp = exp( Z_clin * beta_c + Z_gene * beta_g) ./ (1 + exp(Z_clin * beta_c + Z_gene * beta_g)); 
    
    for (i in 1:N) {
      //estimate log_lik
        lpdf = weibull_lpdf(yobs[i] | alpha, exp(-(mu) / alpha) );
        pdf = exp(lpdf);
        term1 = - (log(lp[i]))*pdf;
        cdf = weibull_cdf(yobs[i], alpha, exp(-(mu) / alpha) );
        term2 = exp(log(lp[i]) * cdf);
        log_lik[i] = log(term1^v[i] * term2);
        
    }
}


          


  

  
  