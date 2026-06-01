################################################
#### GWAS on selected plasma proteins       ####
#### in fasting OGTT study                  ####
################################################

rm(list = ls())
options(stringsAsFactors = F)

## ============================================================
## Configuration: set UKB data paths for your environment
## ============================================================
ukb_pheno_parquet  <- "/path/to/ukb_phenotype_data.parquet"           # Path to UKB participant table (parquet)
ukb_pheno_dict     <- "/path/to/ukb_data_dictionary.txt"              # UKB data dictionary
ukb_olink_imputed_dir <- "/path/to/olink/imputed_data"                # directory with imputed Olink data files
## ============================================================

## faster data handling
require(data.table)
require(arrow)
require(tidyverse)
require(gridExtra)

##################################################
####      import covariate data variables     ####
##################################################

## Mapping table between col ids and names (for reference)
lab.main <- read.table(ukb_pheno_dict, sep = "\t", header = T)

## An example list of columns:
cl.select <- c(
  "f.eid",
  "f.21022.0.0",
  "f.31.0.0",
  "f.54.0.0",
  "f.22000.0.0",
  paste0("f.22009.0.", 1:10)
)
## names to assign
cl.names <- c("f.eid", "age", "sex", "centre", "batch", paste0("pc", 1:10))

## import data from the main release
ukb.dat <- read_parquet(ukb_pheno_parquet, col_select = cl.select)
## change names
names(ukb.dat) <- cl.names

## transform batch to binary
ukb.dat$batch <- ifelse(ukb.dat$batch < 0, 1, 0)
table(ukb.dat$batch)

##################################################
####           get imputed olink data         ####
##################################################
panel_files <- c(
  "Olink.imputed.1.txt",
  "Olink.imputed.2.txt",
  "Olink.imputed.3.txt",
  "Olink.imputed.4.txt"
)

olink_tables <- vector("list", length(panel_files))
for (i in seq_along(panel_files)) {
  olink_tables[[i]] <- read.table(
    file.path(ukb_olink_imputed_dir, panel_files[i]),
    sep = "\t",
    header = T
  )
}

olink.dat <- do.call(cbind, olink_tables)
olink.dat <- olink.dat[, !duplicated(colnames(olink.dat), fromLast = TRUE)] # remove duplicated metadata cols

## move eid column to end
olink.dat$eid <- NULL
olink.dat$eid <- olink_tables[[1]]$eid

rm(olink_tables)

## select relevant cols
prots <- c("anxa10")

olink.dat.sel <- olink.dat[, c(
  "eid",
  "month_blood",
  "time_blood",
  "fast",
  "sample_age",
  prots
)]

##################################################
####      combine ukb dat with olink data     ####
##################################################

## combine with covariate data
ukb.dat.olink <- merge(ukb.dat, olink.dat.sel, by.x = "f.eid", by.y = "eid")

## apply inverse normal transformation (??)
prots_int <- apply(ukb.dat.olink[, prots], 2, function(x) {
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
})
colnames(prots_int) <- paste0(colnames(prots_int), "_int")

ukb.dat.olink <- cbind(ukb.dat.olink, prots_int)

##################################################
####      create input for regenie            ####
##################################################

ukb.dat.olink$FID <- ukb.dat.olink$IID <- ukb.dat.olink$f.eid

## make sex binary
ukb.dat.olink$sex <- ifelse(ukb.dat.olink$sex == "Female", 1, 0)

## write covariates to file
fwrite(
  ukb.dat.olink[, c(
    "FID",
    "IID",
    "age",
    "sex",
    "centre",
    "batch",
    "fast",
    "month_blood",
    "time_blood",
    "sample_age",
    paste0("pc", 1:10)
  )],
  "input/covariates.txt",
  sep = "\t",
  row.names = F,
  quote = F,
  na = NA
)

## write phenotypes to file
fwrite(
  ukb.dat.olink[, c("FID", "IID", prots, colnames(prots_int))],
  "input/phenotypes.txt",
  sep = "\t",
  row.names = F,
  quote = F,
  na = NA
)

## create list of outcome names
fwrite(
  prots,
  "input/outcomes.txt",
  sep = "\t",
  row.names = F,
  quote = F,
  na = NA
)