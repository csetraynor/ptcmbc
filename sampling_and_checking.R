library(Rlab)
library(tidyverse)
library(survival)

# alpha_raw <- 0.2
# tau_al <- 0.5
# log_alpha <- alpha_raw * tau_al
# alpha <- exp(log_alpha)
# print(test_alpha)

## ----sim-data-function-----------------------------------------------
fun_sim_data <- function(alpha, mu, N, beta, gamma , icept, X, a) {
  
  b0 <- icept
  risk <- X %*% beta #prognostic index
  f <- exp(b0 + risk) #link function mixed model
  th <- exp(X %*% gamma) #link function time promotion 
  
  cr <- (1 - a * (f / (1 + a * f) )   )^(round(1/a, 5)) #cure rate
  cure = ifelse(cr > runif(N), T, F) #F "uncured" T "cured"

  spop<- runif(N - sum(cure), min = cr[cure == F], max = 1) #set Surv pop to a unif rv
  
  term1 <- log(1 - ( (1 - (spop)^(a) ) / ( a * f[cure == F] / (1 + a * f[cure == F]))) ) #inverse function
  dfs <- ( - term1 / exp( -(mu + th[cure == F])/alpha) ) ^ (1 / alpha)
  
  data_uncured <- data_frame(dfs_time = dfs) %>% cbind(X[cure == F,])
  data_cured <- data_frame(dfs_time = 9999) %>% cbind(X[cure == T,])
  
  data <- rbind(data_cured, data_uncured) %>% mutate(
    censor_time = rexp(N, 1)) %>%
  dplyr::mutate(
      d = ifelse(dfs_time < censor_time, 'Recurred/Progressed', 'Disease Free')) %>%
  dplyr::mutate( 
    dfs_obs = ifelse(dfs_time < censor_time, dfs_time, censor_time))
  
  return(data)
}
## ----sim-data--------------------------------------------------------

test_X <- matrix(rnorm(2000), ncol = 2)
colnames(test_X) <- c("z1", "z2")
#Survival months simulation
test_icept <- 2
b1 <- 1
b2 <- -1
test_beta <- c(b1, b2)
g1 <- 2
g2 <- -2
test_gamma <- c(g1, g2)
test_N = 1000
test_a <- 0.5
test_alpha <- 1
test_mu <- 1

simd <- fun_sim_data(
  alpha = test_alpha,
  mu = test_mu,
  N = test_N,
  X = test_X,
  beta = test_beta,
  gamma =  test_gamma,
  a = test_a,
  icept = test_icept
)
#Time to Event Distribution#
simd %>% ggplot(aes(x = dfs_obs,
                    group = d,
                    colour = d,
                    fill = d)) + geom_density(alpha = 0.5)
require(ggfortify)
autoplot(survival::survfit(Surv(dfs_obs, I(d == 'Recurred/Progressed')) ~ 1,data = simd), conf.int = F)

## MLE fit

library(miCoPTCM)
fm <- formula(Surv(dfs_months, I(dfs_status == 'Recurred/Progressed')) ~ z1 + z2)
varCov <- matrix(nrow=3,ncol=3,0)
resMY <- PTCMestimBF(fm, simd , varCov = varCov, init=rnorm(3))
resMY
summary(resMY)







