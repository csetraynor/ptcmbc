library(caret)
library(dplyr)
library(rstan)

#--- Cure Model 1 ------#
gen_stan_data <- function(data, eset,  formula = as.formula(~ 1), varoi = NULL) {
  
  if (!inherits(formula, 'formula'))
    formula <- as.formula(formula)
  
  data <- data %>% mutate(
    dfs_progression = (dfs_status == 'Recurred/Progressed')
  )
  
  Z <- data %>% 
    model.matrix(formula, data = . )
  prop <- apply(Z[,-1], 2, sum) / apply(Z[,-1], 2, length) 
  Z[,-1] <- sweep(Z[,-1], 2, prop) #centering
  if("(Intercept)" %in% colnames(Z)){
    Z = Z[,-1] #deselect intercept 
  }
  M <- ncol(Z)
  
  if(!is.null(varoi)){
    var <- data %>% select(varoi) 
    
    varm <- var %>% 
      model.matrix(as.formula(paste( "~" ,varoi, sep = " ")), data = . )
    prop <- apply(varm, 2, sum) / apply(varm, 2, length) 
    varm <- sweep(varm, 2, prop) #centering
    if("(Intercept)" %in% colnames(varm)){
      varm = varm[,-1] #deselect intercept 
    }
  }
  
  Z_gene <- eset
  M_gene <- ncol(Z_gene)
  
  stan_data <- list(
    N = nrow(data),
    yobs = as.numeric(data$dfs_months),
    v = as.numeric(data$dfs_progression),
    M = M,
    Z_c = array(Z, dim = c(nrow(data), M)),
    J = dplyr::n_distinct(data$type),
    type = data$type,
    Z_g = Z_gene,
    M_g = M_gene,
    nu_global = 1, #100 gives Half-Normal Prior, 1 Half-Cauchy
    nu_local =1, #1 corresponds to Horseshoe prior
    scale_global = 1e-4 #corresponds to formulae
  )
}
into_data <- gen_stan_data(md,
                           formula =  '~ stage + nodes + erandpr', eset= brcaES) 
rstan::stan_rdump(ls(into_data), file = "checking.data.R",
                  envir = list2env(into_data))
gen_stan_data(md,
              formula =  '~ stage + nodes + erandpr', eset= brcaES) %>%glimpse

#---Update inits----#
gen_inits <- function(J,M,M_g){
 # function()
    list(
      alpha_raw = 0.01*rnorm(1),
      
      beta_c = rnorm(M),
      
      beta0 = rnorm(J),
      loc0 = rnorm(1),
      s0 = abs(rnorm(1)),
      
      tau1_global = 0.1*abs(rnorm(1)),
      tau2_global = 0.1*abs(rnorm(1)),
      tau1_local = abs(rnorm(M_g)),
      tau2_local = abs(rnorm(M_g)),
      beta_g_raw = rnorm(M_g)
      
    )
}
inits <- gen_inits(M = 6, J = 4, M_g = 17213)
rstan::stan_rdump(ls(inits), file = "checking.init.R",
                  envir = list2env(inits))
################################################################
##Run Stan
stan_file <- "ptcmbc/Hierarchical/hsml.stan"

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
nChain <- 4
memory.size(500000000000)

stanfit <- rstan::stan(stan_file,
                       data =gen_stan_data(md,
                                           formula =  '~ stage + nodes + erandpr',
                                           eset= brcaES),
                       cores = min(nChain, parallel::detectCores()),
                       chains = nChain,
                       iter = 2000,
                       init =  gen_inits(M = 6, J = 4, M_g = 17213))


save(stanfit, file = "clinicalFit.Rdata")


