library(caret)
library(dplyr)
library(rstan)

#--- Cure Model 1 ------#
gen_stan_data <- function(data,  formula = as.formula(~ 1)) {
  
  if (!inherits(formula, 'formula'))
    formula <- as.formula(formula)
  
  data <- data %>% mutate(
    dfs_progression = (dfs_status == 'Recurred/Progressed')
  )

  
  stan_data <- list(
    N = nrow(data),
    yobs = as.numeric(data$dfs_months),
    v = as.numeric(data$dfs_progression),
    type = data$type,
    J = dplyr::n_distinct(data$type)
  )
}

gen_stan_data(md, 
              formula =  ~ 1) %>%glimpse
#---Update inits----#
gen_inits <- function(J){
  function()
    list(
      alpha_raw = 0.01*rnorm(1),
      
      beta0 = rnorm(J),
      loc0 = rnorm(1),
      s0 = abs(rnorm(1))
    )
}

################################################################

##Run Stan
stan_file <- "ptcmbc/Hierarchical/nmcnullml.stan"

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
nChain <- 4

stanfit <- rstan::stan(stan_file,
                         data = gen_stan_data(md,
                                              formula =  ~ 1),
                         cores = min(nChain, parallel::detectCores()),
                         chains = nChain,
                         iter = 2000,
                         init = gen_inits(J = 4))
save(stanfit, file = "NullModel.Rdata")

