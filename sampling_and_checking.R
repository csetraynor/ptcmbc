library(Rlab)
library(tidyverse)
library(survival)

# alpha_raw <- 0.2
# tau_al <- 0.5
# log_alpha <- alpha_raw * tau_al
# alpha <- exp(log_alpha)
# print(test_alpha)

## ----sim-data-function-----------------------------------------------
fun_sim_data <- function(alpha, mu, N, beta , icept, X, a) {
  
  b0 <- icept
  risk <- X %*% beta #prognostic index
  f <- exp(b0 + risk) #link function mixed model
  th <- exp(risk) #link function time promotion 
  
  cr <- (1 + a * f)^-2 #cure rate
  cure = ifelse(cr > runif(N), T, F) #F "uncured" T "cured"

  spop<- runif(N, min = cr, max = 1) #set Surv pop to a unif rv
  
  term1 <- log(1 - (1 - (spop)^(1/2) ) / ( a * f / (1 + a * f))) #inverse function
  term2 <- ( - term1 / (mu + th) ) ^ (1 / alpha)
  
  data <- data.frame(dfs_time = term2,
                     cure_i = cure, #F "uncured" T "cured"
                     censor_time = rexp(N, 1),
                     stringsAsFactors = F)  %>%
  dplyr::mutate(
    dfs_time = ifelse( cure_i == T, 10e10, dfs_time )) %>%  #Inf time for cure
  dplyr::mutate(
      d = ifelse(dfs_time < censor_time, 'Recurred/Progressed', 'Disease Free')) %>%
  dplyr::mutate( 
    dfs_obs = ifelse(dfs_time < censor_time, dfs_time, censor_time)) %>%
    cbind(test_X)
  
  return(data)
}
## ----sim-data--------------------------------------------------------

z1 <- Rlab::rbern(1000, prob = 0.5)
z2 <- stats::runif(1000, 2.1, 3.1)
test_X <- as.matrix(cbind(z1, z2))
#Survival months simulation
test_icept <- 2
b1 <- 1
b2 <- -1
test_beta <- c(b1, b2)
test_N = 1000
test_a <- 0.7
test_alpha <- 1
test_mu <- 0.5

simd <- fun_sim_data(
  alpha = test_alpha,
  mu = test_mu,
  N = test_N,
  X = test_X,
  beta = test_beta,
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







