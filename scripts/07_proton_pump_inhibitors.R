setwd("/sc-projects/sc-proj-computational-medicine/people/Burulca/06_fasting_ogtt/16_proton_pump_inhibitors/")

require(data.table)
#require(arrow)
require(tidyverse)
require(dplyr)
require(lmerTest)
require(readxl)
require(foreach)
require(doMC)
require(ggrepel)
require(rms)

load("16.RData")

####> Does rs2990223 increase ANXA10 plasma levels through
####> gastric disease --> PPI intake --> increased plasma level of stomach proteins

### Test models individually

## trans pQTL effect on anxa10 lvl - lower anxa10
anxa_res <- lm(paste0("ANXA10 ~ rs2990223 + age + sex + ",  paste(paste0("pc", 1:10),collapse=" + "), " + rcs(fast.0, c(.5,3,20)) + rcs(time_blood.num, c(9,14,20)) + rcs(month_blood, c(1,6,12)) + rcs(sample_age, c(11,12.7,15))"), data = filt_df)
summary(anxa_res)$coefficients["rs2990223",c("Pr(>|t|)", "Std. Error", "Estimate")]
## ppi effect - higher anxa10
ppi_res <- lm(paste0("ANXA10 ~ A02BC + age + sex + ", paste(paste0("pc", 1:10),collapse=" + ")), data = filt_df)
summary(ppi_res)$coefficients["A02BC1",c("Pr(>|t|)", "Std. Error", "Estimate")]
## h2 blocker effect - not significant
h2_res <- lm(paste0("ANXA10 ~ A02BA + age + sex + ",  paste(paste0("pc", 1:10),collapse=" + ")), data = filt_df)
summary(h2_res)$coefficients["A02BA1",c("Pr(>|t|)", "Std. Error", "Estimate")]


### Now filter out people taking PPI or H2 blocker
filt_df_no_ppi_h2 <- filt_df %>% filter(acid_med != 1) # n=34,351
anxa_res <- lm(paste0("ANXA10 ~ rs2990223 + age + sex + ",  paste(paste0("pc", 1:10),collapse=" + "), " + rcs(fast.0, c(.5,3,20)) + rcs(time_blood.num, c(9,14,20)) + rcs(month_blood, c(1,6,12)) + rcs(sample_age, c(11,12.7,15))"), data = filt_df_no_ppi)
summary(anxa_res)$coefficients["rs2990223",c("Pr(>|t|)", "Std. Error", "Estimate")]
#>> signal still there

### Now filter out people taking PPI only
filt_df_no_ppi <- filt_df %>% filter(A02BC != 1) # n=34,869
anxa_res <- lm(paste0("ANXA10 ~ rs2990223 + age + sex + ",  paste(paste0("pc", 1:10),collapse=" + "), " + rcs(fast.0, c(.5,3,20)) + rcs(time_blood.num, c(9,14,20)) + rcs(month_blood, c(1,6,12)) + rcs(sample_age, c(11,12.7,15))"), data = filt_df_no_ppi)
summary(anxa_res)$coefficients["rs2990223",c("Pr(>|t|)", "Std. Error", "Estimate")]

### Now filter people taking PPI only
filt_df_ppi <- filt_df %>% filter(A02BC == 1) # n=3488
anxa_res <- lm(paste0("ANXA10 ~ rs2990223 + age + sex + ",  paste(paste0("pc", 1:10),collapse=" + "), " + rcs(fast.0, c(.5,3,20)) + rcs(time_blood.num, c(9,14,20)) + rcs(month_blood, c(1,6,12)) + rcs(sample_age, c(11,12.7,15))"), data = filt_df_ppi)
summary(anxa_res)$coefficients["rs2990223",c("Pr(>|t|)", "Std. Error", "Estimate")]


################################################################################

#### Do all proteins

registerDoMC(cores = 10) # Adjust the number of cores as needed

## Get effect of pqtl on all proteins, excl PPI
res_rs2990223 <- foreach(a = assays, .combine = rbind) %dopar% {
  model <- lm(paste0(a," ~ age + sex + ",  paste(paste0("pc", 1:10),collapse=" + "),
                     " + rs2990223 + rs11589479 + rs4390169 + rcs(fast.0, c(.5,3,20)) + rcs(time_blood.num, c(9,14,20)) + rcs(month_blood, c(1,6,12)) + rcs(sample_age, c(11,12.7,15))"),
              data = filt_df_no_ppi)
  summary_model <- summary(model)
  res = summary_model$coefficients["rs2990223", c("Pr(>|t|)","Std. Error", "Estimate")]
  c(Assay = a, res)
}

# Convert results to a data frame
res_rs2990223_df <- data.frame(res_rs2990223)
colnames(res_rs2990223_df) <- c("Assay","pval","se","Estimate")
res_rs2990223_df$pval <- as.numeric(res_rs2990223_df$pval)
res_rs2990223_df$Estimate <- as.numeric(res_rs2990223_df$Estimate)
res_rs2990223_df$se <- as.numeric(res_rs2990223_df$se)
res_rs2990223_df$fdr <- p.adjust(res_rs2990223_df$pval,method = "BH")

################################################################################

## Get effect of ppis on all proteins, excl PPI and H2 blockers
res_rs2990223 <- foreach(a = assays, .combine = rbind) %dopar% {
  model <- lm(paste0(a," ~ age + sex + ",  paste(paste0("pc", 1:10),collapse=" + "),
                     " + rs2990223 + rs11589479 + rs4390169 + rcs(fast.0, c(.5,3,20)) + rcs(time_blood.num, c(9,14,20)) + rcs(month_blood, c(1,6,12)) + rcs(sample_age, c(11,12.7,15))"),
              data = filt_df_no_ppi)
  summary_model <- summary(model)
  res = summary_model$coefficients["rs2990223", c("Pr(>|t|)","Std. Error", "Estimate")]
  c(Assay = a, res)
}

# Convert results to a data frame
res_rs2990223_df <- data.frame(res_rs2990223)
colnames(res_rs2990223_df) <- c("Assay","pval","se","Estimate")
res_rs2990223_df$pval <- as.numeric(res_rs2990223_df$pval)
res_rs2990223_df$Estimate <- as.numeric(res_rs2990223_df$Estimate)
res_rs2990223_df$se <- as.numeric(res_rs2990223_df$se)
res_rs2990223_df$fdr <- p.adjust(res_rs2990223_df$pval,method = "BH")



## Compare to estimates including people taking PPI
comp.df <- res_rs2990223_df %>% inner_join(pqtl_adj_df, by="Assay", suffix = c(".excl_ppi",".incl_ppi"))

ggplot(comp.df, aes(x=fdr.incl_ppi, y=fdr.excl_ppi)) + geom_point() +
  geom_text_repel(aes(label=Assay), size=2, max.overlaps = 20) +
  theme_minimal()

ggplot(comp.df, aes(x=-log10(fdr.incl_ppi), y=-log10(fdr.excl_ppi))) + geom_point() +
  geom_text_repel(aes(label=Assay), size=2) +
  theme_minimal()

pdf("graphics/incl_vs_excl_ppi_users.pdf", width = 10, height = 10)
ggplot(comp.df, aes(x=Estimate.incl_ppi, y=Estimate.excl_ppi)) + geom_point() +
  #geom_errorbar(aes(xmin = Estimate.incl_ppi - se.incl_ppi, xmax = Estimate.incl_ppi + se.incl_ppi), width = 0.1, alpha=0.1) +  # Add error bars
  #geom_errorbar(aes(ymin = Estimate.excl_ppi - se.excl_ppi, ymax = Estimate.excl_ppi + se.excl_ppi), width = 0.1, alpha=0.1) +  # Add error bars
  geom_text_repel(aes(label=Assay), size=3, max.overlaps = 20) +
  xlab("rs2990223 effect") + ylab("rs2990223 effect excl. PPI intake") +
  theme_minimal()
dev.off()

## Get effect of pqtl on all proteins, excl PPI and H2 blockers
res_rs2990223 <- foreach(a = assays, .combine = rbind) %dopar% {
  model <- lm(paste0(a," ~ age + sex + ",  paste(paste0("pc", 1:10),collapse=" + "),
                     " + rs2990223 + rs11589479 + rs4390169 + rcs(fast.0, c(.5,3,20)) + rcs(time_blood.num, c(9,14,20)) + rcs(month_blood, c(1,6,12)) + rcs(sample_age, c(11,12.7,15))"),
              data = filt_df_no_ppi)
  summary_model <- summary(model)
  res = summary_model$coefficients["rs2990223", c("Pr(>|t|)","Std. Error", "Estimate")]
  c(Assay = a, res)
}