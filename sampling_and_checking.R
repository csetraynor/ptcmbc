library(Rlab)
library(tidyverse)
library(survival)

alpha_raw <- 0.2
tau_al <- 0.5
log_alpha <- alpha_raw * tau_al
test_alpha <- exp(log_alpha)
print(test_alpha)
test_mu <- -3

z1 <- Rlab::rbern(100, prob = 0.5)
z2 <- stats::runif(100, 2.1, 3.1)
test_X <- as.matrix(cbind(z1, z2))
#Survival months simulation
b0 <- 2
b1 <- 1
b2 <- -1
test_beta <- c(b0, b1, b2)

test_N = 1000

## ----sim-data-function-----------------------------------------------
fun_sim_data <- function(alpha, mu, N, beta , X) {
  
  intercept <- beta[1] 
  riski <- X %*% beta[-1]
  #cure probabilities
  cure <- exp(-exp(intercept + riski))
  
  
  
  data <- data.frame(surv_months = rweibull(n = N, alpha, exp(-(mu)/alpha)),
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
## ----sim-data--------------------------------------------------------

simd <- fun_sim_data(
  alpha = test_alpha,
  mu = test_mu,
  N = test_N,
  X = test_X,
  beta = test_beta
)
#Time to Event Distribution#
simd %>% ggplot(aes(x = dfs_months,
                    group = dfs_status,
                    colour = dfs_status,
                    fill = dfs_status)) + geom_density(alpha = 0.5)
autoplot(survival::survfit(Surv(dfs_months, I(dfs_status == 'Recurred/Progressed')) ~ 1,data = simd), conf.int = F)



## MLE fit

library(miCoPTCM)
fm <- formula(Surv(dfs_months, I(dfs_status == 'Recurred/Progressed')) ~ z1 + z2)
varCov <- matrix(nrow=3,ncol=3,0)
resMY <- PTCMestimBF(fm, simd , varCov = varCov, init=rnorm(3))
resMY
summary(resMY)


############################################################################
gen_stan_data <- function(data, formula = as.formula(~ 1)) {
  
  if (!inherits(formula, 'formula'))
    formula <- as.formula(formula)
  
  data <- data %>% mutate(
    dfs_progression = (dfs_status == 'Recurred/Progressed')
  )

  X <- data %>% 
    model.matrix(formula, data = . )
  
  M <- ncol(X)
  
  stan_data <- list(
    N = nrow(data),
    yobs = data$dfs_months,
    v_i = as.numeric(data$dfs_progression),
    M_clinical = M,
    X_clin = array(X, dim = c(nrow(data), M))
  )
}

into_data <- gen_stan_data(simd, '~ z1 + z2')
into_data %>% glimpse


rstan::stan_rdump(ls(into_data), file = "checking.data.R",
                  envir = list2env(into_data))


gen_inits <- function(M){
    list(
      alpha_raw = 0.01*rnorm(1),
      mu = rnorm(1),
      tau_s_clin_raw = 0.1*abs(rnorm(1)),
      tau_clin_raw = array(abs(rnorm(M)), dim = c(M)),
      beta_clin_raw = array(rnorm(M), dim = c(M))
    )
}

inits <- gen_inits(M = 3)
rstan::stan_rdump(ls(inits), file = "checking.inits.R",
                  envir = list2env(inits))
