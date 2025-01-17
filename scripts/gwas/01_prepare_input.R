################################################
#### GWAS on selected plasma proteins       ####
#### in fasting OGTT study                  ####
################################################

rm(list=ls())
setwd("/sc-projects/sc-proj-computational-medicine/people/Burulca/06_fasting_ogtt/04_gwas")
options(stringsAsFactors = F)
load("RData/01_prepare_input.RData")

## faster data handling
require(data.table)
require(arrow)
require(tidyverse)
require(gridExtra)

##################################################
####      import covariate data variables     ####
##################################################

## Mapping table between col ids and names (for reference)
mapping_tbl="/sc-projects/sc-proj-computational-medicine/data/UK_biobank/phenotypes/working_data/Data.dictionary.UKBB.main.dataset.45268.txt"
lab.main <- read.table(mapping_tbl, sep="\t", header=T)

## An example list of columns:
cl.select      <- c("f.eid", "f.21022.0.0", "f.31.0.0", "f.54.0.0", "f.22000.0.0", paste0("f.22009.0.", 1:10))
## names to assign
cl.names       <- c("f.eid", "age", "sex", "centre", "batch", paste0("pc", 1:10))

## import data from the main release
ukb.dat        <- read_parquet("/sc-projects/sc-proj-computational-medicine/data/UK_biobank/phenotypes/working_data/parquet_files/ukb45268.parquet",
                               col_select = cl.select)
## change names
names(ukb.dat) <- cl.names

## transform batch to binary
ukb.dat$batch  <- ifelse(ukb.dat$batch < 0, 1, 0)
table(ukb.dat$batch)

##################################################
####           get imputed olink data         ####
##################################################

setwd("/sc-projects/sc-proj-computational-medicine/data/UK_biobank/proteomics/01_Olink/01_release_explore_15k/imputed_data")

olink1 <- read.table("Olink.release.1.5k.baseline.imputed.panel.1.dataset.1.20230421.txt", sep="\t", header = T)
olink2 <- read.table("Olink.release.1.5k.baseline.imputed.panel.2.dataset.1.20230421.txt", sep="\t", header = T)
olink3 <- read.table("Olink.release.1.5k.baseline.imputed.panel.3.dataset.1.20230421.txt", sep="\t", header = T)
olink4 <- read.table("Olink.release.1.5k.baseline.imputed.panel.4.dataset.1.20230421.txt", sep="\t", header = T)

setwd("/sc-projects/sc-proj-computational-medicine/people/Burulca/06_fasting_ogtt/04_gwas")

olink.dat <- cbind(olink1,olink2,olink3,olink4)
olink.dat <- olink.dat[, !duplicated(colnames(olink.dat), fromLast = TRUE)] # remove duplicated metadata cols

## move eid col to end
olink.dat$eid <- NULL
olink.dat$eid <- olink1$eid

rm(olink1,olink2,olink3,olink4)

## select relevant cols
prots <- c("anxa10","gcg","igfbp1","fabp6","lpl","fgf21","il6","ctsv","mln","wif1","fst","rspo1","wars","spock1","spink5","agxt","il12rb1")
# 17 proteins
olink.dat.sel <- olink.dat[,c("eid",
                              "month_blood","time_blood","fast","sample_age",
                              prots)]

# --> Out of the top 10 significant proteins in both ogtts plus incretins, I selected those that are available in explore 1536
# --> plus AGXT and IL12RB1 as negative controls (not sig in ogtt but has high % above LOD)
# --> E.g. PYY and GAST are on explore expansion panels (cardiometabolic ii and neurology ii), so not in UKBB :(

##################################################
####      combine ukb dat with olink data     ####
##################################################

## combine with covariate data
ukb.dat.olink <- merge(ukb.dat, olink.dat.sel, by.x="f.eid", by.y="eid")

## apply inverse normal transformation (??)
prots_int  <- apply(ukb.dat.olink[, prots], 2, function(x){
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
})
colnames(prots_int)<-paste0(colnames(prots_int),"_int")

ukb.dat.olink <- cbind(ukb.dat.olink, prots_int)

ukb.dat.olink <- ukb.dat.olink[, c(1:19, order(names(ukb.dat.olink)[20:53])+19 ) ]

##################################################
####      draw histograms of the prot data    ####
##################################################

hist_list <- list()

# Plot histograms for the blood metadata and protein cols
for (col_name in names(ukb.dat.olink)[16:53]) {
  hist_data <- ukb.dat.olink[[col_name]]  # Extract the column data
  hist_plot <- ggplot(data.frame(x = hist_data), aes(x)) +
    geom_histogram(bins=30) +
    labs(title = col_name, x = "NPX") +
    theme_minimal() +
    theme(axis.title.y = element_blank())
  
  # Store the histogram plot in the list
  hist_list[[col_name]] <- hist_plot
}

pdf("graphics/phenotypes_histograms.pdf", width=10, height=25)
grid.arrange(grobs = hist_list, ncol = 4, nrow = 10)
dev.off()

##################################################
####      create input for regenie            ####
##################################################

ukb.dat.olink$FID  <- ukb.dat.olink$IID <- ukb.dat.olink$f.eid

## make sex binary
ukb.dat.olink$sex  <- ifelse(ukb.dat.olink$sex == "Female", 1, 0)

## write covariates to file
fwrite(ukb.dat.olink[, c("FID", "IID", "age", "sex", "centre", "batch", "fast", "month_blood", "time_blood", "sample_age", paste0("pc", 1:10))],
       "input/covariates.txt",
       sep = "\t", row.names=F, quote=F, na = NA)

## write phenotypes to file - take only non-INT protein values
fwrite(ukb.dat.olink[, c("FID", "IID", prots, colnames(prots_int))], "input/phenotypes.txt", sep = "\t", row.names=F, quote=F, na = NA)

## create list of outcome names
write.table(prots, "input/outcomes.txt", sep="\t", col.names = F, row.names = F, quote = F)

save.image("RData/01_prepare_input.RData")
