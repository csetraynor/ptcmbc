library(caret)
library(dplyr)
library(rstan)

#--- Cure Model 1 ------#
gen_stan_data <- function(data, ind, formula = as.formula(~ 1)) {
  
  if (!inherits(formula, 'formula'))
    formula <- as.formula(formula)
  
  data <- data %>% mutate(
    dfs_progression = (dfs_status == 'Recurred/Progressed')
  )
  
  data_t = data[-ind,]  #training
  data_h = data[ind,] #holdout

  stan_data <- list(
    N_t = nrow(data_t),
    yobs_t = as.numeric(data_t$dfs_months),
    v_t = as.numeric(data_t$dfs_progression),

    N_h = nrow(data_h),                        #holdout
    yobs_h = as.numeric(data_h$dfs_months),    #holdout
    v_h = as.numeric(data_h$dfs_progression)  #holdout

  )
}

gen_stan_data(md, ind = folds[[1]],
              formula =  ~ 1) %>%glimpse
#---Update inits----#
gen_inits <- function(M){
  function()
    list(
      alpha_raw = 0.01*rnorm(1),
      mu = rnorm(1, 0, 10),

      beta0 = rnorm(1)
    )
}

################################################################
## Create data partition K folds for Cross Validation
K = 10
folds <- caret::createFolds(md$dfs_status, k = K, list = TRUE)

##Run Stan
stan_file <- "ptcmbc/nullcv.stan"

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
nChain <- 4
for( i in 1:K){
  stanfit <- rstan::stan(stan_file,
                         data = gen_stan_data(md, ind = folds[[1]],
                                              formula =  ~ 1),
                         cores = min(nChain, parallel::detectCores()),
                         chains = nChain,
                         iter = 2000,
                         init = gen_inits(M = 0)
  )
  save(stanfit, file = paste("NullModel","CVFold", i, sep="_"))
}
