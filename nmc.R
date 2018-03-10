library(Rlab)
library(tidyverse)
library(survival)

# alpha_raw <- 0.2
# tau_al <- 0.5
# log_alpha <- alpha_raw * tau_al
# alpha <- exp(log_alpha)
# print(test_alpha)

## ----sim-data-function-----------------------------------------------
fun_sim_data <- function(alpha, mu, N, b, g, Z) {
  b <- as.vector(b %>% unlist); g <- as.vector(g %>% unlist)
  alpha <- as.numeric(alpha %>% unlist); mu <- as.numeric(mu %>% unlist)

  assertthat::assert_that(length(g) == length(b) -1)
  
  lp <- Z %*% b #prognostic index
  f <- exp(lp) / (1 + exp(lp) ) #link function mixed model

        #p <- rbeta(N,1,1) 
  #uncured = ifelse( p > f, T, F) #T: "uncured" F: "cured"
  sigma <- exp( - (mu + Z[,-1] %*% g) / alpha)
  U <- runif(N, min = f, max = 1) #for an inverse transform
  X =  sigma * (- log(1 - (log(U) / log(f)) ) )^(1/alpha) 
  
  data <- data.frame( dfs_months = X,
                      censor_months = runif(N, 0, 500),
                      #dfs_status = uncured,
                      lp = f,
                      stringsAsFactors = F) %>% 
   # dplyr::mutate(
    #  dfs_months = ifelse(dfs_status == F,  10e10, dfs_months))%>%
    dplyr::mutate(
      dfs_status = ifelse(dfs_months < censor_months, 'Recurred/Progressed', 'Disease Free')) %>%
    dplyr::mutate( 
      dfs_obs = ifelse(dfs_months < censor_months, dfs_months, censor_months)) %>%
    cbind(Z)
  #Assert that non NaNs produced
  assertthat::assert_that(sum(is.nan(data$dfs_obs)) == 0)
  
  return(data)
}
## ----sim-data--------------------------------------------------------

test_Z <- matrix(rnorm(2000), ncol = 2)
colnames(test_Z) <- c("z1", "z2")
#Survival months simulation
test_icept <- -1
b1 <- 4
b2 <- 6
test_beta <- c(b1, b2);beta= test_beta

test_N = 1000
test_mu <- -8
test_alpha <- 1

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

