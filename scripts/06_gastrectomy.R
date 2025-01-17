setwd("/sc-projects/sc-proj-computational-medicine/people/Burulca/06_fasting_ogtt/00_clean_code/")

library(arrow)
library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)
library(stringr)

####### Mapping table for main ukb table #####
mapping_tbl="/sc-projects/sc-proj-computational-medicine/data/UK_biobank/phenotypes/working_data/Data.dictionary.UKBB.main.dataset.45268.txt"
lab.main <- read.table(mapping_tbl, sep="\t", header=T)

####### Get cols from main UKB table #####
slct_cols <- c("f.eid", "f.21022.0.0", 	"f.53.0.0", "f.31.0.0")
slct_cols_names <- c("f.eid", "age", "date_assesment", "sex")

phenotype_tbl="/sc-projects/sc-proj-computational-medicine/data/UK_biobank/phenotypes/working_data/parquet_files/ukb45268.parquet"
ukbb.dat <- arrow::read_parquet(phenotype_tbl, col_select = all_of(slct_cols)) 
names(ukbb.dat) <- slct_cols_names
ukbb.dat <- as.data.table(ukbb.dat)


####### Get Olink data and join #######
olink <- read.table("/sc-resources/ukb/data/projects/44448/qc_data/blood_assays/proteomics/olink_explore_3k/20231208/qc_data/Olink.release.3k.baseline.20231211.tsv", sep="\t", header = T)
colnames(olink)[2] <- "f.eid"
ukbb.dat <- ukbb.dat %>% left_join(olink[,c("f.eid","anxa10","gast")], by="f.eid")


################################################################################

####### Get gastrectomy participants (full or partial) #########
## (generated the tsv on dnanexus)

gastrectomy <- read.table("/sc-projects/sc-proj-computational-medicine/people/Burulca/06_fasting_ogtt/13_anxa10_gastrectomy/data/gastrectomy.tsv", sep="\t", header = T)
ukbb.dat <- ukbb.dat %>% left_join(gastrectomy, by=c("f.eid"="eid"))

ukbb.dat$operation_type <- case_when( substr(ukbb.dat$oper4, 1, 3)=="G27" ~ "Total gastrectomy",
                                      substr(ukbb.dat$oper4, 1, 3)=="G28" ~ "Partial gastrectomy",
                                      TRUE ~ NA)

ukbb.dat$opdate <- as.Date.character(ukbb.dat$opdate)

## Check if operation date was before assessment center
ukbb.dat$operated_at_baseline <- ukbb.dat$date_assesment > ukbb.dat$opdate
table(ukbb.dat$operated_at_baseline, !is.na(ukbb.dat$anxa10))
#        FALSE TRUE
# FALSE   475   66 --> operated after baseline
# TRUE    124   14 --> operated before baseline

## Filter only those with Olink data
ukb.olink <- ukbb.dat %>% filter(!is.na(anxa10)|!is.na(gast)) # n=80 operated, 50,313 total with olink data

ukb.olink$operated_at_baseline <- case_when(ukb.olink$operated_at_baseline ~ "Gastrectomy before baseline",
                                            !ukb.olink$operated_at_baseline ~ "Gastrectomy after baseline",
                                            TRUE ~ "No gastrectomy")
ukb.olink$operated_at_baseline <- factor(ukb.olink$operated_at_baseline, 
                                       labels = str_wrap(levels(factor(ukb.olink$operated_at_baseline)), width = 10))

# Filter the data to include only the operated group for geom_jitter()
operated_data <- ukb.olink %>% filter(!is.na(operation_type))

pdf("graphics/06_gastrectomy_gast.pdf", height = 6, width = 8)
ggplot(ukb.olink, aes(x = factor(operated_at_baseline), y = gast)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(data = operated_data, aes(color = factor(operation_type)), width = 0.2) + 
  theme_minimal() + labs(color = "OPCS4") +
  theme(axis.title.x = element_blank(), text = element_text(size = 16))
dev.off()

pdf("graphics/06_gastrectomy_anxa.pdf", height = 6, width = 8)
ggplot(ukb.olink, aes(x = factor(operated_at_baseline), y = anxa10)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(data = operated_data, aes(color = factor(operation_type)), width = 0.2) + 
  theme_minimal() + labs(color = "OPCS4") +
  theme(axis.title.x = element_blank(), text = element_text(size = 16))
dev.off()

################################################################################