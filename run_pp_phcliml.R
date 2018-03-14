library(caret)
library(dplyr)
library(rstan)

#--- Cure Model 1 ------#
gen_stan_data <- function(data, formula = as.formula(~ 1), varoi = NULL) {
  
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
  
  stan_data <- list(
    N = nrow(data),
    yobs = as.numeric(data$dfs_months),
    v = as.numeric(data$dfs_progression),
    M = M,
    Z_c = array(Z, dim = c(nrow(data), M)),
    J = dplyr::n_distinct(data$type),
    type = data$type,
    varm = as.numeric(varm)
  )
}

gen_stan_data(md,
              formula =  '~ stage + nodes', varoi = "erandpr") %>%glimpse

into_data <- gen_stan_data(md,
                           formula =  '~ stage + nodes', varoi = "erandpr")
rstan::stan_rdump(ls(into_data), file = "checking.data.R",
                  envir = list2env(into_data))
#---Update inits----#
gen_inits <- function(J,M){
  # function()
  list(
    alpha_raw = 0.01*rnorm(1),
    
    mu_mu = rnorm(1),
    mu_J = rnorm(1),
    s_mu = abs(rnorm(J)),
    
    gI_mu = rnorm(1),
    gI_J = rnorm(1),
    s_gI = abs(rnorm(J)),
    
    beta_c = rnorm(M),
    
    beta0 = rnorm(J),
    loc0 = rnorm(1),
    s0 = abs(rnorm(1)),
    
    betaE = rnorm(J),
    locE = rnorm(1),
    sE = abs(rnorm(1)),
    
    tau_s_g_raw = 0.1*abs(rnorm(1)),
    tau_g_raw = abs(rnorm(M)),
    gamma_g_raw = rnorm(M)
    
  )
}
inits <- gen_inits(M = 5, J = 4)
rstan::stan_rdump(ls(inits), file = "checking.init.R",
                  envir = list2env(inits))
################################################################
##Run Stan
stan_file <- "ptcmbc/Hierarchical/nmclinphml.stan"

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
nChain <- 1

stanfit <- rstan::stan(stan_file,
                       data = gen_stan_data(md,
                                            formula =  '~ stage + nodes',
                                            varoi = "erandpr") ,
                       cores = min(nChain, parallel::detectCores()),
                       chains = nChain,
                       iter = 1000,
                       init = gen_inits(M = 5, J = 4))

save(stanfit,  file = "CureModel")

if (interactive())
  shinystan::launch_shinystan(stanfit)

