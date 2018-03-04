gc()
#----Load libraries---#
library(dplyr)
library(readr)
library(tidyverse)
library(survival)
library(rstan)
library(assertthat)
library(cgdsr)
theme_set(theme_bw())
suppressMessages(library(VIM))
suppressMessages(library(caret))
library(mice)
library(parallel)
library(synapseClient)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source("https://bioconductor.org/biocLite.R")
suppressMessages(library(Biobase))

library(globaltest) #global test
gt.options(trace=FALSE, max.print=45)
library(matrixStats)
suppressMessages(library(Hmisc))
library(affy)  #For gene expression set
library(MSnbase)
#---Load Data----#
dfs.data <- read_csv("Complete_METABRIC_Clinical_Features_Data.txt")
glimpse(dfs.data)
#---Quick Data Preprocess ---#
dfs.data %>% glimpse
dfs.data %>%  filter(is.na(status) | status == "" | is.na(time) | time <= 0) %>%
  glimpse
#Missing data combination plot
dfs.data %>% 
  VIM::aggr(prop = FALSE, combined = TRUE, numbers = TRUE, sortVars = TRUE, sortCombs = TRUE)
# Nottingham Prognostic Index Calculator
calculate_npi <- function(x, stage, grade, size){
    x <- stage + grade + (size * 0.2)
}
dfs.data <- dfs.data %>% mutate(
  npi = if_else(is.na(npi), calculate_npi(stage = stage, grade = grade, size = size), npi))

#--- Distribution event times ---#
dfs.data %>% mutate_at(vars(status), funs(as.factor)) %>%
  ggplot(aes(x = time, group = status, colour = status,  fill = status)) + geom_density(alpha = 0.5)
mle.surv <- survfit(Surv(time, status) ~ 1, data = dfs.data)
require(ggfortify)
ggplot2::autoplot(mle.surv, conf.int = F) +  ggtitle('KM disease free survival for BRCA Cohort')

#--- Prepare data for Stan ---#
## This function will take a formula object as input. 
## Genetic formula is not possible to use due model.matrix too big.
## Order data for observed time then,
## num id is the numerical id for the arranged observations.
## X clin are clinical covariates of interest in the model.
## M clin are the number of cov.
## v is the censoring indicator and "y" observed time

gen_stan_data <- function(data, formula = as.formula(~1) ) {
  if(!inherits(formula, 'formula'))
    formula <- as.formula(formula)
  
  data <- data %>% 
    dplyr::arrange(time) %>%
    mutate(num_id = seq(n()))
  
  X_clin <- data %>%
    model.matrix(formula, data = .)
  M_clinical <- ncol(X_clin)
  
  if (M_clinical > 1){
    if("(Intercept)" %in% colnames(X_clin))
      X_clin <- array(X_clin[,-1], 
                      dim = c(nrow(data), M_clinical - 1)) }
  M_clinical <- ncol(X_clin)
  
  stan_data <- list(
    N =  nrow(data),
    yobs =  data$time,
    num_id = data$num_id,
    v_i = data$status,
    M_clinical =  M_clinical,
    X_clin =  array(X_clin, dim = c(nrow(data), M_clinical))
  )
}

clinical_data <- gen_stan_data(data = dfs.data %>% filter(!is.na(npi)), formula = '~ npi')

rstan::stan_rdump(ls(clinical_data), file = "npi.data.R",
                  envir = list2env(clinical_data))
survreg(Surv(time,status)~npi, data = dfs.data, dist="w")
