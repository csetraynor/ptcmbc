library(dplyr)
library(readr)
library(caret)
library(mice)
library(survival)
#memory.size(400000)
#download data

sampledat <- read_tsv("brca_data/brca_tcga_pub2015/data_clinical_sample.txt", skip = 4)
patientdat <- read_tsv("brca_data/brca_tcga_pub2015/data_clinical_patient.txt", skip = 4)
gendat <- read_tsv("brca_data/brca_tcga_pub2015/data_expression_median.txt", col_names = TRUE)
# Get an easier code to read
names(sampledat) <- tolower(names(sampledat))
names(patientdat) <- tolower(names(patientdat))
md <- left_join(patientdat, sampledat) #joint clinical data
rm(patientdat, sampledat)
md <- md %>% filter(sample_id %in% intersect(sample_id, colnames(gendat))) 
#---- Data Cleaning ----#
#convert missig values into NA
convert_blank_to_na <- function(x) {
  if(!purrr::is_character(x)){
    warning('input variate not character - return original')
    return(x)
  } else {
    ifelse(x == "[Not Available]" , NA, x)
  }
}
md <- md %>% dplyr::mutate_all(funs(convert_blank_to_na))
glimpse(md)
md %>%
  VIM::aggr(prop = FALSE, combined = TRUE, numbers = TRUE, sortVars = TRUE, sortCombs = TRUE, plot = FALSE, only.miss = FALSE)
md %>% filter(is.na(dfs_status) | dfs_months <= 0) %>% select(dfs_status, dfs_months) %>% glimpse

md <- md %>% filter(!is.na(dfs_status), dfs_months > 0) %>% #remove NA observation
  select(
  patient_id, sample_id, cancer_type_detailed, age, history_other_malignancy, ajcc_pathologic_tumor_stage, ajcc_nodes_pathologic_pn, histological_diagnosis, surgical_procedure_first, er_status_by_ihc, pr_status_by_ihc, dfs_status, dfs_months, her2_ihc_score
) %>% mutate_at(vars(c("age", "dfs_months")), funs(as.numeric))
#Time to Event Distribution#
md %>% ggplot(aes(x = dfs_months,
                  group = dfs_status,
                  colour = dfs_status,
                  fill = dfs_status)) + geom_density(alpha = 0.5)
require(ggfortify)
autoplot(survival::survfit(Surv(dfs_months, I(dfs_status == 'Recurred/Progressed')) ~ 1 ,data = md), conf.int = F)
#Imputation#
md %>%
  VIM::aggr(prop = FALSE, combined = TRUE, numbers = TRUE, sortVars = TRUE, sortCombs = TRUE)
md$dfs_months <- as.numeric(md$dfs_months)
md$ajcc_pathologic_tumor_stage <- as.factor(md$ajcc_pathologic_tumor_stage)


VIM::marginplot(md[c("ajcc_pathologic_tumor_stage","dfs_months")])
table(md$ajcc_pathologic_tumor_stage, useNA = "always")
#Following NPI and AJCC guidelines
md$stage <- NA
md$stage[grepl("I$|IA$|IB$",md$ajcc_pathologic_tumor_stage )] <- "1"
md$stage[grepl("II$|IIA$|IIB$",md$ajcc_pathologic_tumor_stage )] <- "2"
md$stage[grepl("IIIA$|IIIA$|IIIC$",md$ajcc_pathologic_tumor_stage )] <- "3"
md$stage[grepl("IV$|X$",md$ajcc_pathologic_tumor_stage )] <- "4"
#using imputation by Bayesian poly regression
tmp <- as.factor(md$stage)
tmp <- mice::mice.impute.polyreg(y = tmp, 
                                 ry = !is.na(tmp),
                            x = model.matrix(~ ajcc_nodes_pathologic_pn  + dfs_status,
                                                  data = md)[,-1],
                                 wy = array(TRUE, dim = length(tmp)))
md$stage[is.na(md$stage)] <- tmp[is.na(md$stage)]
remove(tmp)

md$nodes <- NA
md$nodes[grepl("0",md$ajcc_nodes_pathologic_pn)] <- "1"
md$nodes[grepl("1",md$ajcc_nodes_pathologic_pn)] <- "2"
md$nodes[grepl("2|3|X",md$ajcc_nodes_pathologic_pn)] <- "3"


#--- Gene matrix preprocess ----- #
library(Biobase)
glimpse(gendat)
gene_names <- gendat %>% select(Hugo_Symbol)  #grap gene names
gendat <- gendat %>% select(intersect(colnames(gendat), md$sample_id))# get intersection with clinical and gene values
sample_names <- colnames(gendat) # get sample names
sum(is.na(gendat))

#Convert to expression set
md <- as.data.frame(md %>% filter(sample_id %in% sample_names) %>% slice(match(sample_names, sample_id))) ; row.names(md) <- md$sample_id#requires classical data frame
x <-  as.matrix(gendat) ;colnames(x) <- rownames(md); 
gene_names <- as.data.frame(gene_names);  rownames(gene_names) <- gene_names %>% unlist
brcaES <- Biobase::ExpressionSet(x,
                                  phenoData = as(md, "AnnotatedDataFrame"),
                                 featureData = as(gene_names, "AnnotatedDataFrame"))
rm(gendat)
gene_names <- gene_names %>% unlist

##Perform empirical Bayes to find differential gene expressions
library(limma)
fit <- limma::lmFit(brcaES)
fit <- limma::eBayes(fit)
### x is the input data. This function replaces the top 'perc' percent
### with the value 'rp'. 
perc=0.1
qnt = quantile(fit$F.p.value, 1-perc)
w = which(fit$F.p.value >= qnt)
subset = (brcaES)[w]

subset.top = function(x, perc, rp)
{
  qnt = quantile(x,1-perc)
  w = which(x >= qnt) 
  x[w] = rp
  return(x)
}
# Treat relative to a 10% fold-change
tfit <- treat(fit,lfc=log2(1.1))
topTreat(tfit,coef=2)

#Imputation with the min value
rm(x)
require(MSnbase)
brcaMSN <- MSnbase::as.MSnSet.ExpressionSet(brcaES)
brcaMSN <- MSnbase::impute(brcaMSN, method = "min")
Biobase::exprs(brcaES) <- MSnbase::exprs(brcaMSN)
rm(brcaMSN)
sum(is.na(Biobase::exprs(brcaES)))

#Obtain the Z scores
preProcgm <-  caret::preProcess(t(exprs(brcaES)), method = c("center", "scale")) 
brcaES <- predict(preProcgm, t(exprs(brcaES))) 
rm(preProcgm)



