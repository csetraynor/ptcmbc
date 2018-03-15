#--- Gen Stan Data ---#
gen_stan_data <- function(data) {
  observed_data <- data %>%
    dplyr::filter(os_status == 'DECEASED')
  Jobs = observed_data$cohort
  
  censored_data <- data %>%
    dplyr::filter(os_status != 'DECEASED')
  Jcen = censored_data$cohort
  
  assertthat::assert_that(n_distinct(Jobs) == n_distinct(Jcen))

  stan_data <- list(
    Nobs = nrow(observed_data),
    Ncen = nrow(censored_data),
    yobs = observed_data$os_months,
    ycen = censored_data$os_months,
    J = n_distinct(data$cohort),
    Jobs = as.numeric(Jobs),
    Jcen = as.numeric(Jcen)
  )
}
into_data <- gen_stan_data(md)
rstan::stan_rdump(ls(into_data), file = "checking.data.R",
                  envir = list2env(into_data))

gen_inits <- function(J) {
  list(
    alpha_raw = 0.01*rnorm(1),
    theta = rnorm(J),
    mu_theta = rnorm(1),
    tau_theta = abs(rnorm(1))
  )
}
inits <- gen_inits(J = 6)
rstan::stan_rdump(ls(inits), file = "checking.init.R",
                  envir = list2env(inits))
