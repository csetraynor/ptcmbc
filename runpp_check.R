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


#-----Run Stan PH-------#
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
nChain <- 4
stanfile <- 'ptcmbc/nmcph.stan'
ph_fit <- stan(stanfile,
                data = gen_stan_data(md, '~ stage + nodes'),
                init = gen_inits(M = 6),
                iter = 1000,
                thin = 1,  
                cores = min(nChain, parallel::detectCores()),
                seed = 7327,
                chains = nChain)
 if (interactive())
   shinystan::launch_shinystan(ph_fit)


#---Posterior predictive checks---#
pp_beta <- as.data.frame.array(rstan::extract(ph_fit,pars = 'beta_clin', permuted = TRUE)$beta_clin) 
pp_gamma <- as.data.frame.array(rstan::extract(ph_fit,pars = 'gamma_clin', permuted = TRUE)$gamma_clin)
pp_alpha <- as.data.frame.array(rstan::extract(ph_fit,pars = 'alpha', permuted = TRUE)$alpha) 
pp_mu <- as.data.frame.array(rstan::extract(ph_fit,pars = 'mu', permuted = TRUE)$mu) 
# create list
pp_beta <-  split(pp_beta, seq(nrow(pp_beta)))
pp_gamma <-  split(pp_gamma, seq(nrow(pp_gamma)))
pp_mu <-  split(pp_mu, seq(nrow(pp_mu)))
pp_alpha <-  split(pp_alpha, seq(nrow(pp_alpha)))
Z_clin <- as.matrix(into_data$Z_clin, ncol = into_data$M_clinical)

pp_newdata <- 
  purrr::pmap(list(pp_beta, pp_gamma, pp_alpha, pp_mu),
              function(pp_beta, pp_gamma, pp_alpha, pp_mu) {fun_sim_data(
                                                         b = pp_beta,
                                                         g = pp_gamma, 
                                                         alpha = pp_alpha,
                                                         mu = pp_mu,
                                                         N = nrow(md),
                                                         Z = Z_clin)} )
ggplot(pp_newdata %>%
         dplyr::bind_rows() %>%
         dplyr::mutate(type = 'posterior predicted values') %>%
         bind_rows(md %>% dplyr::mutate(type = 'actual data')),
       aes(x = dfs_months, group = dfs_status,
             colour = dfs_status, fill = dfs_status)) +
  geom_density(alpha = 0.5) + xlim(c(0, 1000)) +
  facet_wrap(~type, ncol = 1)

pp_predict_surv <- function(beta, gamma, alpha, mu, N, Z,
                            level = 0.9, 
                            plot = F, data = NULL,
                            fun_sim = fun_sim_data,
                            formula = as.formula(~1)) {

  pp_newdata <- purrr::pmap(list(beta, gamma, alpha, mu),
                function(beta, gamma, alpha, mu) {fun_sim(
                  b = beta,
                  g = gamma, 
                  alpha = alpha,
                  mu = mu,
                  N = N,
                  Z = Z)} )
  
  pp_survdata <-
    pp_newdata %>%
    purrr::map(~ dplyr::mutate(., dfs_deceased = dfs_status == "Recurred/Progressed")) %>%
    purrr::map(~ survival::survfit(Surv(dfs_months, dfs_deceased) ~ 1, data = .)) %>%
    purrr::map(fortify)
  ## compute quantiles given level 
  lower_p <- 0 + ((1 - level)/2)
  upper_p <- 1 - ((1 - level)/2)
  pp_survdata_agg <- 
    pp_survdata %>%
    purrr::map(~ dplyr::mutate(.,
                               time_group = floor(time))) %>%
    dplyr::bind_rows() %>%
    dplyr::group_by(time_group) %>%
    dplyr::summarize(surv_mean = mean(surv)
                     , surv_p50 = median(surv)
                     , surv_lower = quantile(surv,
                                             probs = lower_p)
                     , surv_upper = quantile(surv,
                                             probs = upper_p) ) %>%
    dplyr::ungroup()
  if (plot == FALSE) {
    return(pp_survdata_agg)
  } 
  ggplot_data <- pp_survdata_agg %>%
    dplyr::mutate(type = 'posterior predicted values') %>%
    dplyr::rename(surv = surv_p50,
                  lower = surv_lower,
                  upper = surv_upper, time = time_group)
  if (!is.null(data))
    ggplot_data <- 
    ggplot_data %>% 
    bind_rows(
      fortify(
        survival::survfit(
          Surv(dfs_months, dfs_deceased) ~ 1, 
          data = data %>% 
            dplyr::mutate(
              dfs_deceased = dfs_status == "Recurred/Progressed") )) %>%
        dplyr::mutate(lower = surv,
                      upper = surv, type = 'actual data') )
  pl <- ggplot(ggplot_data,
               aes(x = time, group = type, linetype = type)) + 
    geom_line(aes(y = surv, colour = type)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2)
  
  pl 
}

pl <- pp_predict_surv(beta = pp_beta,
                      gamma = pp_gamma,
                      alpha = pp_alpha,
                      mu = pp_mu,
                      N = nrow(md),
                      Z = Z_clin,
                      level = 0.9, 
                      plot = T, 
                      data = md) 
pl + 
  xlim(NA, 300) +
  ggtitle('Posterior predictive checks \nfit to clinical data showing 90% CI')



#-----Run Stan Horseshoe-------#
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
nChain <- 4
stanfile <- 'ptcmbc/nmchs.stan'
ph_fit <- stan(stanfile,
               data = gen_stan_data(data = md, Eset = brcaES, formula ='~ stage + nodes' ),
               init = gen_inits(M_clinical = 6, M_gene = 17213),
               iter = 1000,
               thin = 1,  
               cores = min(nChain, parallel::detectCores()),
               seed = 7327,
               chains = nChain)

