#-----Run Stan-------#
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
nChain <- 2
stanfile <- 'ptcmbc/nmc_loglik.stan'
sim_fit <- stan(stanfile,
                data = gen_stan_data(md, '~ stage + nodes'),
                init = gen_inits(M = 3),
                iter = 2000,
                thin = 1,  
                cores = min(nChain, parallel::detectCores()),
                seed = 7327,
                chains = nChain)



#  control = list(adapt_delta = 0.95),
# pars = c("alpha", "mu",  "lp__", "beta_gene", "beta_clin"))

print(sim_fit)
