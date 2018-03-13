##### Run Stan Clinical Model
#setwd('c:/RFactory/ptcm_project/ptcmbc')
#load(data)
############################################################################
gen_stan_data <- function(data, formula = as.formula(~ 1)) {
  
  if (!inherits(formula, 'formula'))
    formula <- as.formula(formula)
  
  data <- data %>% mutate(
    dfs_progression = (dfs_status == 'Recurred/Progressed')
  )

  stan_data <- list(
    N = nrow(data),
    yobs = as.numeric(data$dfs_months),
    v = as.numeric(data$dfs_progression)
  )
}
#---NULL MODEL Set initial Values ---#
gen_inits <- function(){
  function()
    list(
      alpha_raw = 0.01*rnorm(1),
      mu = rnorm(1, 0, 10)
    )
}
############################################################################
#--- Cure Model 1 ------#
gen_stan_data <- function(data, formula = as.formula(~ 1)) {
  
  if (!inherits(formula, 'formula'))
    formula <- as.formula(formula)
  
  data <- data %>% mutate(
    dfs_progression = (dfs_status == 'Recurred/Progressed')
  )
  
  Z <- data %>% 
    model.matrix(formula, data = . )
  prop <- apply(Z_clin[,-1], 2, sum) / apply(Z_clin[,-1], 2, length) 
  Z_clin[,-1] <- sweep(Z_clin[,-1], 2, prop)

  M <- ncol(Z)
  
  stan_data <- list(
    N = nrow(data),
    yobs = as.numeric(data$dfs_months),
    v = as.numeric(data$dfs_progression),
    M_clinical = M,
    Z_clin = array(Z, dim = c(nrow(data), M))
  )
}

into_data <- gen_stan_data(md, '~ stage + nodes + erandpr')
into_data %>% glimpse
rstan::stan_rdump(ls(into_data), file = "checking.data.R",
                  envir = list2env(into_data))
rm(into_data)

inits <- gen_inits(M = 6)
rstan::stan_rdump(ls(inits), file = "checking.inits.R",
                  envir = list2env(inits))


#---CLINICAL MODEL MO Update initial Values ---#
gen_inits <- function(M){
  function()
    list(
      alpha_raw = 0.01*rnorm(1),
      mu = rnorm(1, 0, 10),
      
      tau_s_cb_raw = abs(rnorm(1)),
      tau_cb_raw = array(abs(rnorm(M)), dim = c(M)),
      beta_clin_raw = array(rnorm(M), dim = c(M))
    )
}

#---CLINICAL MODEL with PH Update initial Values ---#
gen_inits <- function(M){
  function()
  list(
    alpha_raw = 0.01*rnorm(1),
    
    mu = rnorm(1, 0, 10),
    
    tau_s_cb_raw = 0.1*abs(rnorm(1)),
    tau_cb_raw = array(abs(rnorm(M)), dim = c(M)),
    beta_clin_raw = array(rnorm(M), dim = c(M)),
    
    tau_s_cg_raw = 0.1*abs(rnorm(1)),
    tau_cg_raw = array(abs(rnorm(M-1)), dim = c(M-1)),
    gamma_clin_raw = array(rnorm(M-1), dim = c(M-1))
  )
}

#Update gen data for Genomic vars
############################################################################
gen_stan_data <- function(data, Eset, formula = as.formula(~ 1), ... ) {
  
  if (!inherits(formula, 'formula'))
    formula <- as.formula(formula)
  
  data <- data %>% mutate(
    dfs_progression = (dfs_status == 'Recurred/Progressed')
  )
  
  Z_clin <- data %>% 
    model.matrix(formula, data = . )
  prop <- apply(Z_clin[,-1], 2, sum) / apply(Z_clin[,-1], 2, length) 
  Z_clin[,-1] <- sweep(Z_clin[,-1], 2, prop)
  
  M_clin <- ncol(Z_clin)
  
  Z_gene <- Eset

  M_gene <- ncol(Z_gene)
  
  stan_data <- list(
    N = nrow(data),
    yobs = as.numeric(data$dfs_months),
    v = as.numeric(data$dfs_progression),
    M_c = M_clin,
    M_g = M_gene,
    Z_clin = Z_clin,
    Z_gene = Z_gene,
    nu_global = 100, #1 gives Half-Normal Prior
    nu_local =1, #1 corresponds to Horseshoe prior
    scale_global = 0.0013 #corresponds to formulae
  )
}

into_data <- gen_stan_data(data = md, Eset = brcaES, formula ='~ stage + nodes + erandpr' )
into_data %>% glimpse
rstan::stan_rdump(ls(into_data), file = "checking.data.R",
                  envir = list2env(into_data))
rm(into_data)
#---Horseshoe Update initial Values ---#
gen_inits <- function(M_c, M_g){
  function()
  list(
    alpha_raw = 0.01*rnorm(1),
    mu = rnorm(1, 0, 10),
    
    tau_s_cb_raw = 0.1*abs(rnorm(1)),
    tau_cb_raw = abs(rnorm(M_c)),
    beta_c_raw = rnorm(M_c),
    
    tau1_global = 0.1*abs(rnorm(1)),
    tau2_global = 0.1*abs(rnorm(1)),
    tau1_local = abs(rnorm(M_g)),
    tau2_local = abs(rnorm(M_g)),
    beta_g_raw = rnorm(M_g)
  )
}
inits <- gen_inits(M_c = 6, M_g = 1722)
rstan::stan_rdump(ls(inits), file = "checking.inits.R",
                  envir = list2env(inits))
