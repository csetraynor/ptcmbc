#### Run Stan Clinical Model
#setwd('c:/RFactory/ptcm_project/ptcmbc')
#load(data)
############################################################################
gen_stan_data <- function(data, formula = as.formula(~ 1)) {
  
  if (!inherits(formula, 'formula'))
    formula <- as.formula(formula)
  
  data <- data %>% mutate(
    dfs_progression = (dfs_status == 'Recurred/Progressed')
  )
  
  Z <- data %>% 
    model.matrix(formula, data = . )
  
  M <- ncol(Z)
  
  stan_data <- list(
    N = nrow(data),
    yobs = data$dfs_obs,
    v_i = as.numeric(data$dfs_progression),
    M_clinical = M,
    Z_clin = array(Z, dim = c(nrow(data), M))
  )
}

into_data <- gen_stan_data(simd, '~ z1 + z2')
into_data %>% glimpse


 rstan::stan_rdump(ls(into_data), file = "checking.data.R",
                   envir = list2env(into_data))
#--- Set initial Values ---#
gen_inits <- function(M){
  function()
    list(
      alpha_raw = 0.01*rnorm(1),
      mu = rnorm(1, 0, 10),
      
      tau_s_cb_raw = 0.1*abs(rnorm(1)),
      tau_cb_raw = array(abs(rnorm(M)), dim = c(M)),
      beta_clin_raw = array(rnorm(M), dim = c(M))
    )
}
 inits <- gen_inits(M = 3, N = 1000)
 rstan::stan_rdump(ls(inits), file = "checking.inits.R",
                   envir = list2env(inits))
#-----Run Stan-------#
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
nChain <- 1
stanfile <- 'nmc.stan'
sim_fit <- stan(stanfile,
                data = gen_stan_data(simd, '~ z1 + z2'),
                init = gen_inits(M = 3),
                iter = 2000,
                thin = 1,  #applying a thin of 10 for avoiding the autocorrelation in baseline hazard
                cores = min(nChain, parallel::detectCores()),
                seed = 7327,
                chains = nChain)
                #  control = list(adapt_delta = 0.95),
               # pars = c("alpha", "mu",  "lp__", "beta_gene", "beta_clin"))

print(sim_fit)


