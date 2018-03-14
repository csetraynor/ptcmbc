library(dplyr)
library(readr)
library(caret)
library(mice)
library(survival)
memory.size(400000000)
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
    ifelse(x == "[Not Available]" | x == "[Discrepancy]", NA, x)
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
                            x = model.matrix(~ ajcc_nodes_pathologic_pn  + dfs_status + cancer_type_detailed +dfs_months,
                                                  data = md)[,-1],
                                 wy = array(TRUE, dim = length(tmp)))
md$stage[is.na(md$stage)] <- tmp[is.na(md$stage)]
remove(tmp)

md$nodes <- NA
md$nodes[grepl("0",md$ajcc_nodes_pathologic_pn)] <- "1"
md$nodes[grepl("1",md$ajcc_nodes_pathologic_pn)] <- "2"
md$nodes[grepl("2|3|X",md$ajcc_nodes_pathologic_pn)] <- "3"

#--- Impute er_status and pr_status
md$erandpr <- "Negative"
md$erandpr[md$er_status_by_ihc == "Positive" & md$pr_status_by_ihc == "Positive"] <- "Positive"
md$erandpr[is.na(md$er_status_by_ihc) | md$pr_status_by_ihc == "Indeterminate"] <- NA


tmp <- as.factor(md$erandpr)
VIM::marginplot(md[c("erandpr","dfs_months")])
tmp <- mice::mice.impute.polyreg(y = tmp, 
                                 ry = !is.na(tmp),
                                 x = model.matrix(~dfs_status+ 
                                                    cancer_type_detailed + dfs_months ,
                                                  data = md)[,-1],
                                 wy = array(TRUE, dim = length(tmp)))
md$erandpr[is.na(md$erandpr)] <- tmp[is.na(md$erandpr)]
remove(tmp)

#--- Impute HER2 score
md$her2_ihc_score[is.na(md$her2_ihc_score)] <- "NA"

#--- Cancer type
md$type <- as.numeric(as.factor(md$cancer_type_detailed))

#----Download Metabric Data

#--- Gene matrix preprocess ----- #
source("https://bioconductor.org/biocLite.R")
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
assertthat::assert_that(all(md$patient_id == brcaES$patient_id))
rm(list = c("gendat", "x"))
gene_names <- gene_names %>% unlist
#Imputation using MSnbase
require(MSnbase)
brcaMSN <- MSnbase::as.MSnSet.ExpressionSet(brcaES)
brcaMSN <- MSnbase::impute(brcaMSN, method = "knn")
Biobase::exprs(brcaES) <- MSnbase::exprs(brcaMSN)
rm(brcaMSN)
sum(is.na(Biobase::exprs(brcaES)))

##Perform empirical Bayes to find differential gene expressions
library(limma)
fit <- limma::eBayes(limma::lmFit(brcaES))
volcanoplot(fit)
# toptable(fit)
# # 
rm(fit)
#Center and scale 
preProcgm <-  caret::preProcess(t(exprs(brcaES)), method = c("center", "scale")) 
brcaES <- predict(preProcgm, t(exprs(brcaES))) 
rm(preProcgm)

# # ### x is the input data. This function replaces the top 'perc' percent
# # ### with the value 'rp'. 
# # 
# subset.top = function(x, perc, eset)
# {
#   qnt = quantile(x$lods,1-perc)
#   w = which(fit$lods >= qnt)
#   return(eset[w])
# }
# brcaEStest = subset.top(x = fit, perc = .2 , eset = brcaES)
# volcanoplot(fit)
# fit_gene_names = rownames(brcaEStest)
# rm(fit)


