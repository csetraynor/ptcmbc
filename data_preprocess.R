library(dplyr)
library(readr)
library(caret)
library(mice)
library(survival)
memory.size(400000000)
#download data
sampletcga <- read_tsv("brca_data/brca_tcga_pub2015/data_clinical_sample.txt", skip = 4)
patienttcga <- read_tsv("brca_data/brca_tcga_pub2015/data_clinical_patient.txt", skip = 4)
gentcga <- read_tsv("brca_data/brca_tcga_pub2015/data_expression_median.txt", col_names = TRUE)

samplemeta <- read_tsv("brca_data/brca_metabric/data_clinical_sample.txt", skip = 4)
patientmeta <- read_tsv("brca_data/brca_metabric/data_clinical_patient.txt", skip = 4)
genmeta <- read_tsv("brca_data/brca_metabric/data_expression.txt", col_names = TRUE)

# Get an easier code to read
names(sampletcga) <- tolower(names(sampletcga))
names(patienttcga) <- tolower(names(patienttcga))
names(samplemeta) <- tolower(names(samplemeta))
names(patientmeta) <- tolower(names(patientmeta))

#joint medical data
mdtcga <- left_join(patienttcga, sampletcga) 
mdmeta <- left_join(patientmeta,samplemeta)

mdmeta <- mdmeta %>% 
  select(patient_id, sample_id, histological_subtype,
         age_at_diagnosis, tumor_stage,
         er_ihc, pr_status, her2_status, 
         os_months, os_status, 
         inferred_menopausal_state,
         chemotherapy,  radio_therapy,
         threegene, cohort) %>%
  rename(histology = histological_subtype,
         age = age_at_diagnosis,
         stage = tumor_stage,
         er = er_ihc,
         pr = pr_status,
         her2 = her2_status,
         menopause = inferred_menopausal_state,
         radio = radio_therapy,
         chemo = chemotherapy)

mdtcga <- mdtcga %>%  
  select(patient_id, sample_id, histological_diagnosis,
         age, ajcc_pathologic_tumor_stage,
         er_status_by_ihc, pr_status_by_ihc, ihc_her2, 
         os_months, os_status,
         menopause_status,  pharmaceutical_tx_adjuvant,
         radiation_treatment_adjuvant ) %>% 
  rename(histology = histological_diagnosis,
         stage = ajcc_pathologic_tumor_stage,
         er = er_status_by_ihc,
         pr = pr_status_by_ihc,
         her2 = ihc_her2,
         menopause = menopause_status,
         chemo = pharmaceutical_tx_adjuvant,
         radio = radiation_treatment_adjuvant) %>%
  mutate( threegene = NA , cohort = "6")
assertthat::assert_that(all(colnames(mdtcga) == colnames(mdmeta)))
md <- rbind(mdmeta , mdtcga)
rm(patienttcga, sampletcga, patientmeta, samplemeta, mdmeta, mdtcga)
glimpse(md)
#---- Data Cleaning ----#
#convert missig values into NA
convert_blank_to_na <- function(x) {
  if(!purrr::is_character(x)){
    warning('input variate not character - return original')
    return(x)
  } else {
    ifelse(x == " " | x == "[Not Available]" | x == "[Discrepancy]", NA, x)
  }
}
md <- md %>% dplyr::mutate_all(funs(convert_blank_to_na)) %>%
  mutate_at(vars(c("age", "os_months")), funs(as.numeric)) 

md %>% VIM::aggr(prop = FALSE, combined = TRUE, numbers = TRUE,
            sortVars = TRUE, sortCombs = TRUE, plot = TRUE, only.miss = FALSE)
md %>% 
  filter(is.na(os_status) | os_status == "" |os_months <= 0 | is.na(os_months)) %>%
  select(os_months, os_status) %>%
  glimpse
#--- Remove Missing Obs ---#
clinical_data <- md
md <- md %>% 
  filter(!is.na(os_status) , os_status != "" , os_months > 0 , !is.na(os_months))
#confirm 544 fewer observations
assertthat::assert_that(nrow(md) == (nrow(clinical_data) - 544))
#--- Distribution event times ---#
clinical_data <- clinical_data %>%
  arrange(os_months)
clinical_data %>%
  ggplot(aes(x = os_months,
             group = os_status,
             colour = os_status,
             fill = os_status)) +
  geom_density(alpha = 0.5)
mle.surv <- survfit(Surv(os_months, os_deceased) ~ cohort,
                    data = clinical_data %>%
                      mutate(os_deceased = (os_status == "DECEASED")))
require(ggfortify)
ggplot2::autoplot(mle.surv, conf.int = F) +
  ggtitle('KM estimate for BCs Cohort')

#---- Pipe ---#

glimpse(md)
md$er[md$er == "Positive" |md$er == "+" ] <- "pos"
md$er[md$er == "Negative" | md$er == "-"] <- "neg"
md$er[md$er == "Indeterminate" | is.na(md$er)] <- "Unavailable" 

md$pr[md$pr == "Positive" |md$pr == "+"] <- "pos"
md$pr[md$pr == "Negative" | md$pr == "-"] <- "neg"
md$pr[md$pr == "Indeterminate" | is.na(md$pr)] <- "Unavailable" 

md$her2[md$her2 == "Positive" |md$her2 == "+" ] <- "pos"
md$her2[md$her2 == "Negative" | md$her2 == "-"] <- "neg"
md$her2[md$her2 == "Indeterminate" | md$her2 == "Equivocal"  | is.na(md$her2) ] <- "Not Available"  

md$menopause[md$menopause == "Peri (6-12 months since last menstrual period)" |
               md$menopause == "Indeterminate (neither Pre or Postmenopausal)"] <- "xperi"
md$menopause[md$menopause == "Post (prior bilateral ovariectomy OR >12 mo since LMP with no prior hysterectomy)"] <- "post"
md$menopause[md$menopause == "Pre (<6 months since LMP AND no prior bilateral ovariectomy AND not on estrogen replacement)"] <- "pre"
md$menopause[is.na(md$menopause)] <- "Unavailable"

md$stage[grepl("I$|IA$|IB$",md$stage )] <- "1"
md$stage[grepl("II$|IIA$|IIB$",md$stage )] <- "2"
md$stage[grepl("IIIA$|IIIA$|IIIC$",md$stage )] <- "3"
md$stage[grepl("IV$",md$stage )] <- "4"
md$stage[grepl("X$", md$stage)] <- "null"
md$stage[is.na(md$stage)] <- "Unavailable"





md %>% filter(is.na(dfs_status) | dfs_months <= 0) %>% select(dfs_status, dfs_months) %>% glimpse

md <- md %>% filter(!is.na(dfs_status), dfs_months > 0) %>% #remove NA observation
  select(patient_id, sample_id, cancer_type_detailed, age, history_other_malignancy, ajcc_pathologic_tumor_stage, ajcc_nodes_pathologic_pn, histological_diagnosis, surgical_procedure_first, er_status_by_ihc, pr_status_by_ihc, dfs_status, dfs_months, ihc_her2, menopause_status) %>% mutate_at(vars(c("age", "dfs_months")), funs(as.numeric))
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
md$stage[is.na(md$stage)] <- NA
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

#--- HER2 score
md$ihc_her2[is.na(md$ihc_her2) | md$ihc_her2 == "Equivocal" | md$ihc_her2 == "Indeterminate" ] <- NA
tmp <- as.factor(md$ihc_her2)

tmp <- mice::mice.impute.polyreg(y = tmp,
                                 ry = !is.na(tmp),
                                 x = model.matrix(~dfs_status+
                                                    cancer_type_detailed + dfs_months ,
                                                  data = md)[,-1],
                                 wy = array(TRUE, dim = length(tmp)))
md$ihc_her2[is.na(md$ihc_her2)] <- tmp[is.na(md$ihc_her2)]
remove(tmp)


#--- Cancer type
cancer_type <- as.factor(md$cancer_type_detailed)
levels(cancer_type)
md$type <- as.numeric(cancer_type)

#----Download Metabric Data
md <- md %>% filter(sample_id %in% intersect(sample_id, colnames(gentcga))) 
#--- Gene matrix preprocess ----- #
source("https://bioconductor.org/biocLite.R")
library(Biobase)
glimpse(gentcga)
gene_names <- gentcga %>% select(Hugo_Symbol)  #grap gene names
gentcga <- gentcga %>% select(intersect(colnames(gentcga), md$sample_id))# get intersection with clinical and gene values
sample_names <- colnames(gentcga) # get sample names
sum(is.na(gentcga))

#Convert to expression set
md <- as.data.frame(md %>% filter(sample_id %in% sample_names) %>% slice(match(sample_names, sample_id))) ; row.names(md) <- md$sample_id#requires classical data frame
x <-  as.matrix(gentcga) ;colnames(x) <- rownames(md); 
gene_names <- as.data.frame(gene_names);  rownames(gene_names) <- gene_names %>% unlist
brcaES <- Biobase::ExpressionSet(x,
                                  phenoData = as(md, "AnnotatedDataFrame"),
                                 featureData = as(gene_names, "AnnotatedDataFrame"))
assertthat::assert_that(all(md$patient_id == brcaES$patient_id))
rm(list = c("gentcga", "x"))
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
# brcaES = subset.top(x = fit, perc = .2 , eset = brcaES)
# volcanoplot(fit)
# fit_gene_names = rownames(brcaEStest)
# rm(f
