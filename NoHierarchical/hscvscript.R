library(caret)
library(dplyr)
library(rstan)

#--- Cure Model 1 ------#
gen_stan_data <- function(data, Eset, ind, formula = as.formula(~ 1)) {
  
  if (!inherits(formula, 'formula'))
    formula <- as.formula(formula)
  
  data <- data %>% mutate(
    dfs_progression = (dfs_status == 'Recurred/Progressed')
  )
  
  Z <- data %>% 
    model.matrix(formula, data = . )
  prop <- apply(Z[,-1], 2, sum) / apply(Z[,-1], 2, length) 
  Z[,-1] <- sweep(Z[,-1], 2, prop) #centering
  M <- ncol(Z) - 1
  Z_gene <- Eset
  M_gene <- ncol(Z_gene)
  
  data_t = data[-ind,]  #training
  data_h = data[ind,] #holdout
  Z_t = Z[-ind,-1] #training
  Z_h = Z[ind,-1] #holdout
  icept_t = Z[-ind,1] #training
  icept_h = Z[ind, 1] #holdout

  
  stan_data <- list(
    N_t = nrow(data_t),
    yobs_t = as.numeric(data_t$dfs_months),
    v_t = as.numeric(data_t$dfs_progression),
    M_clinical = M,
    M_g = M_gene,
    Z_clin_t = array(Z_t, dim = c(nrow(data_t), M)),
    icept_t = as.numeric(icept_t),
    Z_gene_t = Z_gene[-ind,],
    N_h = nrow(data_h),                        #holdout
    yobs_h = as.numeric(data_h$dfs_months),    #holdout
    v_h = as.numeric(data_h$dfs_progression),  #holdout
    Z_clin_h = array(Z_h, dim = c(nrow(data_h), M)),
    icept_h = as.numeric(icept_h),
    Z_gene_h = Z_gene[ind,],
    
    nu_global = 100, #1 gives Half-Normal Prior
    nu_local =1, #1 corresponds to Horseshoe prior
    scale_global = 1e-4 #corresponds to formulae
  )
}

gen_stan_data(md, Eset = brcaES, ind = folds[[1]],
              formula =  '~ stage + nodes + erandpr ') %>%glimpse
#---Update inits----#
gen_inits <- function(M, M_g){
  function()
    list(
      alpha_raw = 0.01*rnorm(1),
      mu = rnorm(1, 0, 10),
      beta_clin = rnorm(M),
      beta0 = rnorm(1),
      
      tau1_global = 0.1*abs(rnorm(1)),
      tau2_global = 0.1*abs(rnorm(1)),
      tau1_local = abs(rnorm(M_g)),
      tau2_local = abs(rnorm(M_g)),
      beta_g_raw = rnorm(M_g)
    )
}

################################################################
## Create data partition K folds for Cross Validation
K = 10
folds <- caret::createFolds(md$dfs_status, k = K, list = TRUE)

##Run Stan
stan_file <- "ptcmbc/hscv.stan"

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
nChain <- 4
for( i in 1:K){
  stanfit <- rstan::stan(stan_file,
                         data = gen_stan_data(data = md, Eset = brcaES, ind = folds[[i]],
                                              formula =  '~ stage + nodes + erandpr'),
                         cores = min(nChain, parallel::detectCores()),
                         chains = nChain,
                         iter = 2000,
                         init = gen_inits(M = 6, M_g = 17213)
  )
  save(stanfit, folds[[i]], file = paste("HorseShoeModel","CVFold", i, sep="_"))
}


