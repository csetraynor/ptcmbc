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
  
  Z <- data %>% 
    model.matrix(formula, data = . )
  M <- ncol(Z)
  
  stan_data <- list(
    N = nrow(data),
    yobs = as.numeric(data$dfs_months),
    v = as.numeric(data$dfs_progression),
    M_clinical = M,
    Z_clin = array(Z, dim = c(nrow(data), M))
  )
}

into_data <- gen_stan_data(md, '~ stage + nodes')
into_data %>% glimpse


rstan::stan_rdump(ls(into_data), file = "checking.data.R",
                  envir = list2env(into_data))
#---NULL MODEL Set initial Values ---#
gen_inits <- function(M){
  function()
  list(
    alpha_raw = 0.01*rnorm(1),
    mu = rnorm(1, 0, 10),
    beta_clin_raw = array(rnorm(M), dim = c(M))
  )
}
inits <- gen_inits(M = 6)
rstan::stan_rdump(ls(inits), file = "checking.inits.R",
                  envir = list2env(inits))


#---CLINICAL MODEL MO Update initial Values ---#
gen_inits <- function(M){
  #function()
    list(
      alpha_raw = 0.01*rnorm(1),
      mu = rnorm(1, 0, 10),
      
      tau_s_cb_raw = 0.1*abs(rnorm(1)),
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
  attr(Z_clin, "dimnames") <- NULL
  attr(Z_clin, "assign") <- NULL
  attr(Z_clin, "contrasts") <- NULL
  M_clin <- ncol(Z_clin)
  
  Z_gene <- t(exprs(Eset))
  attr(Z_gene, "dimnames") <- NULL
  M_gene <- ncol(Z_gene)
  
  Z <- cbind(Z_clin, Z_gene)
  M <- ncol(Z)
  
  
  stan_data <- list(
    N = nrow(data),
    yobs = as.numeric(data$dfs_months),
    v = as.numeric(data$dfs_progression),
    M = M,
    Z = Z,
    nu_global = 1, #1 gives Half-Cauchy Prior
    nu_local =1, #1 corresponds to Horseshoe prior
    scale_global = 0.0002 #corresponds to formulae
  )
}

into_data <- gen_stan_data(data = md, Eset = brcaES, formula ='~ stage + nodes' )
into_data %>% glimpse
rm(into_data)
#---Horseshoe Update initial Values ---#
gen_inits <- function(M){
  function()
  list(
    alpha_raw = 0.01*rnorm(1),
    mu = rnorm(1, 0, 10),
    beta0 = rnorm(1),
    
    tau1_global = 0.1*abs(rnorm(1)),
    tau2_global = 0.1*abs(rnorm(1)),
    tau1_local = abs(rnorm(M)),
    tau2_local = abs(rnorm(M)),
    beta_g_raw = rnorm(M)
  )
}

