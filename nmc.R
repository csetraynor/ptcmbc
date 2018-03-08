library(Rlab)
library(tidyverse)
library(survival)

# alpha_raw <- 0.2
# tau_al <- 0.5
# log_alpha <- alpha_raw * tau_al
# alpha <- exp(log_alpha)
# print(test_alpha)


num.samples <- 1000
U <- runif(num.samples, f, 1)
X <- 50 * (- log (1 - (log(U) / log(f))) ) ^(1/alpha) ; plot(X)

## ----sim-data-function-----------------------------------------------
fun_sim_data <- function(alpha, mu, N, b, icept, Z) {

  b0 <- icept
  risk <- Z %*% beta #prognostic index
  f <- exp(b0 + risk) / (1 + exp(b0 + risk) ) #link function mixed model
  uncured = ifelse( f > runif(N), T, F) #T: "uncured" F: "cured"
  
  U <- runif(N, min = f, max = 1) #set Surv pop to a unif rv
  sigma <- exp( - (mu) / alpha)
  #inverse function
  X <-  sigma * (- log (1 - (log(U) / log(f))) ) ^(1/alpha)  
  
  data <- data.frame( dfs_months = X,
                      censor_months = rexp(n = N, rate = 1/100),
                      uncured_status = uncured,
                      lp = f,
                      stringsAsFactors = F) %>% 
    dplyr::mutate(
      dfs_months = ifelse(uncured_status == F,  10e10, dfs_months))%>%
    dplyr::mutate(
      dfs_status = ifelse(dfs_months < censor_months, 'Recurred/Progressed', 'Disease Free')) %>%
    dplyr::mutate( 
      dfs_obs = ifelse(dfs_months < censor_months, dfs_months, censor_months)) %>%
    cbind(Z)
  
  return(data)
}
## ----sim-data--------------------------------------------------------

test_Z <- matrix(rnorm(2000), ncol = 2)
colnames(test_Z) <- c("z1", "z2")
#Survival months simulation
test_icept <- 1
b1 <- 1
b2 <- -1
test_beta <- c(b1, b2)
g1 <- 2
g2 <- -2
test_gamma <- c(g1, g2)
test_N = 1000
test_mu <- -10
test_alpha <- 1.8

simd <- fun_sim_data(
  alpha = test_alpha,
  mu = test_mu,
  N = test_N,
  Z = test_Z,
  b = test_beta,
  icept = test_icept
)
#Time to Event Distribution#
simd %>% ggplot(aes(x = dfs_obs,
                    group = dfs_status,
                    colour = dfs_status,
                    fill = dfs_status)) + geom_density(alpha = 0.5)
require(ggfortify)
autoplot(survival::survfit(Surv(dfs_obs, I(dfs_status == 'Recurred/Progressed')) ~ 1,data = simd), conf.int = F)
