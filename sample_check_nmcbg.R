library(Rlab)
library(tidyverse)
library(survival)
set_theme

## ----sim-data-function-----------------------------------------------
fun_sim_data <- function(alpha, mu, N, beta, gamma, X,...) {
  
  assertthat::assert_that(length(gamma) == length(beta) - 1);
  riski <- X %*% beta
  #cure probabilities
  cure <- exp(-exp(riski))

  data <- data.frame(surv_months = rweibull(n = N, alpha, exp(-(mu + X[,-1] %*% gamma) /alpha)),
                     censor_months = rexp(n = N, rate = 1/10),
                     stringsAsFactors = F, 
                     cure_prob = cure,
                     random_prob = runif(N)
  ) %>%
    dplyr::mutate(
      surv_months = ifelse(random_prob < cure_prob, 99999, surv_months), # cured subjects
      dfs_status = ifelse(surv_months < censor_months,
                          'Recurred/Progressed', 'Disease Free'
      ),
      dfs_months = ifelse(surv_months < censor_months,
                          surv_months, censor_months
      )
    ) %>% cbind(test_X)
  
  
  return(data)
}

test_alpha <- 1
test_mu <- -5
test_N = 1000

z2 <- matrix(stats::rnorm(200000), ncol = 200) #simulate gene vars
colnames(z2) <- paste('gene', 1:200, sep = "_")
cli_1 <- Rlab::rbern(1000, prob = 0.5) #simulate clinical vars
cli_2 <- stats::runif(1000, -1, 2)
test_X <- as.matrix(cbind(rep(1, test_N), cli_1, cli_2, z2)) #joint vars

gene_beta <- runif(200, min = -.0001, max = .0001 ) #Coefficients for genes
gene_effect <- sample(1:200, 10)
gene_true_effect <- sample(c(seq(from = -3, to = -1.5, by = 0.5), seq(from = 1.5, to = 4, by = 0.5)) , 10, replace = TRUE)
gene_beta[gene_effect] <- gene_true_effect
gene_gamma <- runif(200, min = -.0001, max = .0001 ) #Coefficients for genes
gene_effect <- sample(1:200, 10)
gene_true_effect <- sample(c(seq(from = -3, to = -1.5, by = 0.5), seq(from = 1.5, to = 5, by = 0.5)) , 10, replace = TRUE)
gene_gamma[gene_effect] <- gene_true_effect

bcli_1 = 0 #Coefficients for clinical vars
bcli_2 = -2
bint = 3 #Coefficient intercept
test_beta <- c(bint, bcli_1, bcli_2, gene_beta)
gcli_1 = 1
gcli_2 = 0
test_gamma = c(gcli_1, gcli_2, gene_gamma)

test_alpha <- 5
test_mu <- -15
## ----sim-data--------------------------------------------------------
simd <- fun_sim_data(
  alpha = test_alpha,
  mu = test_mu,
  N = test_N,
  X = test_X,
  beta = test_beta,
  gamma = test_gamma
)
#Time to Event Distribution#
simd %>% ggplot(aes(x = dfs_months,
                    group = dfs_status,
                    colour = dfs_status,
                    fill = dfs_status)) + geom_density(alpha = 0.5)

mle.surv <- survfit(Surv(dfs_months, dfs_recurred) ~ 1,
                    data = simd %>%
                      mutate(dfs_recurred = (dfs_status == 'Recurred/Progressed')))
ggplot2::autoplot(mle.surv, conf.int = F) +
  ggtitle('KM survival for GGM Cohort')
colnames(simd[7:207]) <- paste('var', colnames(simd[7:207]), sep = "_")

## MLE fit
library(miCoPTCM)
fm <- formula(Surv(dfs_months, I(dfs_status == 'Recurred/Progressed')) ~ z1 + z2)
varCov <- matrix(nrow=3,ncol=3,0)
resMY <- PTCMestimBF(fm, simd , varCov = varCov, init=rnorm(3))
resMY
summary(resMY)


############################################################################
# This formula admits tag classical when work
gen_stan_data <- function(data, clinical_formula = as.formula(~ 1), genomic_formula = ...) {
  
  data <- data %>% mutate(
    dfs_progression = (dfs_status == 'Recurred/Progressed')
  )
  
  if (!inherits(clinical_formula, 'clinical_formula'))
    clinical_formula <- as.formula(clinical_formula)
  X_clin <- data %>% 
    model.matrix(clinical_formula, data = . )
  M_clinical <- ncol(X_clin)
  
  X_gene <- data %>% 
    select(contains('gene')) %>%
    as.matrix()
  M_genomic <- ncol(X_gene)
  
  stan_data <- list(
    N = nrow(data),
    yobs = data$dfs_months,
    v_i = as.numeric(data$dfs_progression),
    M_clinical = M_clinical,
    M_genomic = M_genomic,
    X_clin = array(X_clin, dim = c(nrow(data), M_clinical)),
    X_gene = array(X_gene, dim = c(nrow(data), M_genomic))
  )
}

into_data <- gen_stan_data(simd, clinical_formula = '~cli_1 + cli_2', genomic_formula )
into_data %>% glimpse

#--- Set Inits ----#
gen_inits <- function(M_clinical, M_genomic){
  function()
    list(
      alpha_raw = mean(0.01*rnorm(1000)),
      tau_mu = mean(abs(rcauchy(1000, 0, 2))),
      mu = mean(rnorm(1000, 0, 10)),
      
      tau_s_clin_raw = 0.1*abs(rnorm(1)),
      tau_clin_raw = array(abs(rnorm(M_clinical)), dim = c(M_clinical)),
      tau_s1_gene_raw = 0.1*abs(rnorm(1)),
      tau_s2_gene_raw = 0.1*abs(rnorm(1)),
      tau_gene_raw = abs(rnorm(M_genomic)),
      tau1_gene_raw = abs(rnorm(M_genomic)),
      tau2_gene_raw = abs(rnorm(M_genomic)),
      beta_clin_raw = array(rnorm(M_clinical), dim = c(M_clinical)),
      beta_gene_raw = rnorm(M_genomic),
      
      tau_gam_s_clin_raw = 0.1*abs(rnorm(1)),
      tau_gam_clin_raw = array(abs(rnorm(M_clinical-1)), dim = c(M_clinical-1)),
      tau_gam_s1_gene_raw = 0.1*abs(rnorm(1)),
      tau_gam_s2_gene_raw = 0.1*abs(rnorm(1)),
      tau_gam_gene_raw = abs(rnorm(M_genomic)),
      tau_gam1_gene_raw = abs(rnorm(M_genomic)),
      tau_gam2_gene_raw = abs(rnorm(M_genomic)),
      gamma_clin_raw = array(rnorm(M_clinical-1), dim = c(M_clinical-1)),
      gamma_gene_raw = rnorm(M_genomic)
    )
}

#-----Run Stan-------#
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
nChain <- 4
stanfile <- 'ptcmbc/nmcm.hs.bg.stan'
sim_fit2 <- stan(stanfile,
                data = gen_stan_data(simd, clinical_formula = '~cli_1 + cli_2', genomic_formula),
                init = gen_inits(M_clinical = 3, M_genomic = 200),
                iter = 1000,
                thin = 1,  #applying a thin of 10 for avoiding the autocorrelation in baseline hazard
                cores = min(nChain, parallel::detectCores()),
                seed = 7327,
                chains = nChain,
                #  control = list(adapt_delta = 0.95),
                pars = c("alpha", "mu",  "lp__", "beta_gene", "beta_clin"))

print(sim_fit,   pars = c( "alpha", "mu",  "lp__"))
if (interactive())
  shinystan::launch_shinystan(sim_fit)


# rstan::stan_rdump(ls(into_data), file = "checking.data.R",
#                   envir = list2env(into_data))
# 
# rstan::stan_rdump(ls(inits), file = "checking.inits.R",
#                   envir = list2env(inits))

#---- Using WAIC ----#
library(loo)
loo_sim <- loo::loo(extract_log_lik(sim_fit, parameter_name = "lp__"))
