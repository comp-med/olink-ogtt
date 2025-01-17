###################################################
#### Analysis of Olink fasting OGTT study      ####
#### clinical data                             ####
#### Burulca Uluvar 10.07.2024                 ####
###################################################

require(readxl)
require(data.table)
require(ggplot2)
require(gridExtra)
require(lmerTest)
require(dplyr)
require(tidyr)
require(stringr)
require(RColorBrewer)

################## Read clinical data ##########################################

## Excel with glucose, insulin, ffa, BUT
ogtt.dat <- data.frame(read_excel("data/JJ_OGTT_Fast_Round2_Glucose_FFA_BUT_(2019-7-23)_cleaned.xlsx"))
ogtt.dat_long <- ogtt.dat %>% pivot_longer(cols = starts_with("X"), names_to = "t.point", values_to = "measurement")
ogtt.dat_long$t.point <- as.numeric(gsub("X", "", ogtt.dat_long$t.point))
ogtt.dat_long$biomarker <- gsub("ffa", "free fatty acids", ogtt.dat_long$biomarker)

## create time also as factor
ogtt.dat_long$ogtt.session <- as.factor(ogtt.dat_long$ogtt.session)

################## Write summary to file #######################################

ogtt.summ <-ogtt.summ <- ogtt.dat_long %>% group_by(ogtt.session,t.point,biomarker) %>%
  summarise(mean=mean(measurement,na.rm=T),
            sd=sd(measurement,na.rm=T))

write.table(ogtt.summ, "output/02_clinical_mean_and_sd.txt", sep="\t", row.names = F)

################## Convert units ###############################################

ogtt.dat_long2 <- ogtt.dat_long

ogtt.dat_long2$measurement[ogtt.dat_long2$biomarker == "insulin"] <- 
  ogtt.dat_long2$measurement[ogtt.dat_long2$biomarker == "insulin"] / 6

ogtt.dat_long2$unit[ogtt.dat_long2$biomarker == "insulin"] <- "mIU/L"

ogtt.dat_long2$measurement[ogtt.dat_long2$biomarker == "glucose"] <- 
  ogtt.dat_long2$measurement[ogtt.dat_long2$biomarker == "glucose"] * 18.018

ogtt.dat_long2$unit[ogtt.dat_long2$biomarker == "glucose"] <- "mg/dl"

################## Plot conf intervals for each biomarker ######################

# biomarker order for the plot
biomarkers = c("insulin", "glucose", "free fatty acids", "3-hydroxybutyrate")
exclude <- c("FP28") # participants to exclude - FP21 no proteomic data in ogtt2, fp28 no protemic data

ci_plots <- list()

for (b in biomarkers) {
  
  mean_data <- ogtt.dat_long2 %>% filter(biomarker == b,
                                         !is.na(measurement),
                                         !participant %in% exclude,
                                         t.point != 180) %>%
    group_by(t.point, ogtt.session, unit) %>%
    summarize(
      mean_value = mean(measurement, na.rm = TRUE),
      se = qt(0.975, df = sum(!is.na(measurement)) - 1) * sd(measurement, na.rm = TRUE) / sqrt(sum(!is.na(measurement))),
      n = sum(!is.na(measurement))
    )
  # insulin is in picomolar
  #unit <- ifelse(b=="insulin","uU_", "mM")
  
  # Determine if the plot should display the legend
  show_legend <- b == biomarkers[2]
  
  unit = unique(mean_data$unit)
  
  p <- ggplot(mean_data, aes(x = t.point, y = mean_value, col=ogtt.session)) +
    geom_line() +
    geom_point() + 
    geom_errorbar(aes(ymin=mean_value-se, ymax=mean_value+se), width=0) +
    scale_color_manual(values=c("#3288BD","#D53E4F")) +
    theme_minimal() +
    labs(x = "Time [min]", col="OGTT session") + ggtitle(paste0( str_to_title(b), " [", unit, "]" )) +
    theme(panel.grid.minor = element_blank(),
          legend.position = ifelse(show_legend, "bottom", "none"),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 16),
          plot.tag.position = "top",
          text = element_text(size = 16)) +
    scale_x_continuous(breaks = unique(mean_data$t.point))
  
  ci_plots[[b]] <- p
}

#svg("graphics/02_ogtt_biomarkers_ci.svg", width = 3, height = 6)
svg("graphics/02_ogtt_biomarkers_ci_units.svg", width = 3, height = 6)
grid.arrange(grobs = ci_plots[1:2], ncol=1, heights = c(1, 1.2)) # Adjust last height slightly for legend space
dev.off()

################## Compute HOMA-IR and QUICKI during fasting ###################

## Participant data - Fasting gluc & insulin & BMI etc
fasting.dat <- data.frame(read_excel("data/Participant_data_fasting_study_20220718.xlsx"))
fasting.dat <- fasting.dat %>% dplyr::filter(!participant%in%c("FP21","FP28")) # FP28 not in ogtt

## convert units
fasting.dat$insulin_uU_ml <- fasting.dat$insulin / 6
fasting.dat$glucose_mg_dl <- fasting.dat$glucose * 18.018

## compute indices
fasting.dat$homa_ir <- fasting.dat$insulin_uU_ml * fasting.dat$glucose / 22.5
fasting.dat$homa_b <- (20*fasting.dat$insulin_uU_ml) / (fasting.dat$glucose-3.5)
fasting.dat$quicki <-  1 / ( log10(fasting.dat$insulin_uU_ml) + log10(fasting.dat$glucose_mg_dl) )

fasting.dat <- fasting.dat %>% filter(t.point != 10)

## Plot
svg("graphics/02_homa_ir_fasting.svg", width = 8, height = 4)
ggplot(fasting.dat, aes(x = factor(t.point), y = homa_ir, group = t.point))+
  geom_boxplot() + theme_minimal() + theme(text = element_text(size=20)) +
  labs(x="Day",y="HOMA-IR")
dev.off()

svg("graphics/02_quicki_fasting.svg", width = 8, height = 4)
ggplot(fasting.dat, aes(x = factor(t.point), y = quicki, group = t.point))+
  geom_boxplot() + theme_minimal() + theme(text = element_text(size=20)) +
  labs(x="Day",y="QUICKI")
dev.off()

################## Prepare df to compute insulin indices from OGTT #############

## Re-read the data
## Participant data - OGTT gluc & insulin
ogtt.dat <- subset(ogtt.dat, !participant%in%c("FP28","FP21") & biomarker %in% c("glucose","insulin")) ## exclude participant w/out ogtt

## Pivot wider to calculate indices per row
names(ogtt.dat)[2:8] <- sub("^X", "", names(ogtt.dat)[2:8])
ogtt.dat <- ogtt.dat[,-11] %>% pivot_wider(id_cols=c("participant","ogtt.session"), names_from = "biomarker", values_from=2:8, names_glue = "{.name}{.value}")
colnames(ogtt.dat) <- sub(".*_", "", colnames(ogtt.dat))

## Insulin is in pmol/L, glucose is in mmol/L in the dfs till now
## Add insulin as μIU/ml and glucose as mg/dL - This is needed for formulas for indices later

insulin_cols <- colnames(ogtt.dat)[startsWith(colnames(ogtt.dat), "insulin")]
for (col in insulin_cols) {
  new_col_name <- paste0(col, "_uU_ml")
  ogtt.dat[[new_col_name]] <- ogtt.dat[[col]] / 6
}
gluc_cols <- colnames(ogtt.dat)[startsWith(colnames(ogtt.dat), "glucose")]
for (col in gluc_cols) {
  new_col_name <- paste0(col, "_mg_dl")
  ogtt.dat[[new_col_name]] <- ogtt.dat[[col]] * 18.018
}
# add the mean ins & gluc
ogtt.dat$ins_uU_ml_mean <- rowMeans(ogtt.dat[, c("insulin0_uU_ml", "insulin30_uU_ml", "insulin60_uU_ml", "insulin90_uU_ml", "insulin120_uU_ml")])
ogtt.dat$gluc_mg_dl_mean <- rowMeans(ogtt.dat[, c("glucose0_mg_dl", "glucose30_mg_dl", "glucose60_mg_dl", "glucose90_mg_dl", "glucose120_mg_dl")])


################## Compute Stumvoll indices ####################################

## https://diabetesjournals.org/care/article/24/4/796/23380/Oral-Glucose-Tolerance-Test-Indexes-for-Insulin
## Insulin (Ins) measured in picomoles per liter; glucose (Gluc) measured in millimoles per liter.

# Insulin sensitivity index (μmol · kg–1 · min–1 · pmol/l)
ogtt.dat$isi_0_120 <- 0.156 - 0.0000459*ogtt.dat$insulin120 -  0.000321*ogtt.dat$insulin0 - 0.00541*ogtt.dat$glucose120
ogtt.dat$isi_0_60 <- 0.149 - 0.000467*ogtt.dat$insulin0 - 0.00466*ogtt.dat$glucose60

# Metabolic clearance rate (ml · kg–1 · min–1)
ogtt.dat$mcr_0_120 <- 13.273 - 0.00384*ogtt.dat$insulin120 - 0.0232*ogtt.dat$insulin0 - 0.463*ogtt.dat$glucose120
ogtt.dat$mcr_0_60 <- 12.464 - 0.0357*ogtt.dat$insulin0 - 0.376*ogtt.dat$glucose60

# Stumvoll 1st phase
ogtt.dat$ph1_0_30 <- 1283 + 1.829*ogtt.dat$insulin30 - 138.7*ogtt.dat$glucose30 + 3.772*ogtt.dat$insulin0
ogtt.dat$ph1_0_60 <- 1194 + 4.724*ogtt.dat$insulin0 - 117.0*ogtt.dat$glucose60 + 1.414*ogtt.dat$insulin60

# Stumvoll 2nd phase
ogtt.dat$ph2_0_30 <- 286 + 0.416*ogtt.dat$insulin30 - 25.94*ogtt.dat$glucose30 + 0.926*ogtt.dat$insulin0
ogtt.dat$ph2_0_60 <- 295 + 0.349*ogtt.dat$insulin60 - 25.72*ogtt.dat$glucose60 + 1.107*ogtt.dat$insulin0


################## AUC and incremental AUC #####################################

## AUC minutes 0 to 120
ogtt.dat$auc_gluc_0_120 <- 0.125*ogtt.dat$glucose0 + 0.25*ogtt.dat$glucose15 + 0.375*ogtt.dat$glucose30 + 0.5*ogtt.dat$glucose60 + 0.5*ogtt.dat$glucose90 + 0.25*ogtt.dat$glucose120

## Incremental AUC

# first get the timepoint until which glucose is higher than min 0
gluc_nadir <- ogtt.dat[,c(1:14)] %>%
  mutate(gluc_nadir_tpoint = case_when( ## determine which tpoint glucose not decrease below baseline
    glucose30 < glucose0 ~ "glucose15",
    glucose60 < glucose0 ~ "glucose30",
    glucose90 < glucose0 ~ "glucose60",
    glucose120 < glucose0 ~ "glucose90",
    glucose120 > glucose0 ~ "glucose120",
    TRUE ~ NA  # If none of the conditions are met, return NA
  ))

# Replace values with 0 for glucose columns after the nadir timepoint
gluc_nadir <- gluc_nadir %>%
  mutate(across(starts_with("glucose"), ~ifelse(as.numeric(sub("glucose", "", cur_column())) > as.numeric(sub("glucose", "", gluc_nadir_tpoint)), 0, .)))

gluc_nadir$hours_to_nadir <- case_when( 
  gluc_nadir$gluc_nadir_tpoint == "glucose30" ~ 0.5,
  gluc_nadir$gluc_nadir_tpoint == "glucose60" ~ 1,
  gluc_nadir$gluc_nadir_tpoint == "glucose120" ~ 2,
  TRUE ~ 2  # If none of the conditions are met, return NA
)

# incremental AUC minutes 0 to 120
ogtt.dat$iauc_gluc_0_120 <- 0.125*gluc_nadir$glucose0 +
  0.25*gluc_nadir$glucose15 +
  0.375*gluc_nadir$glucose30 +
  0.5*gluc_nadir$glucose60 +
  0.5*gluc_nadir$glucose90 +
  0.25*gluc_nadir$glucose120 -
  gluc_nadir$hours_to_nadir*gluc_nadir$glucose0


################## HOMA and QUICKI #############################################

## Compute HOMA-IR for both OGTTs
ogtt.dat$homa_ir <- ogtt.dat$insulin0_uU_ml*ogtt.dat$glucose0/22.5
ogtt.dat$homa_b <- (20*ogtt.dat$insulin0_uU_ml)/(ogtt.dat$glucose0-3.5)
ogtt.dat$quicki <- 1 / ( log10(ogtt.dat$insulin0_uU_ml) + log10(ogtt.dat$glucose0_mg_dl) )

################## Hepatic insulin resistance index ############################

## HIRI was calculated as: ([glucoseAUC 0–30 in mmol/L*h] × [insulinAUC 0–30 in pmol/L*h]). - Gijbels 2024
ogtt.dat$auc_ins_0_30 <- 0.125*ogtt.dat$insulin0 + 0.25*ogtt.dat$insulin15 + 0.125*ogtt.dat$insulin30
ogtt.dat$auc_gluc_0_30 <- 0.125*ogtt.dat$glucose0 + 0.25*ogtt.dat$glucose15 + 0.125*ogtt.dat$glucose30

## HIRI
ogtt.dat$hiri <- ogtt.dat$auc_ins_0_30 * ogtt.dat$auc_gluc_0_30 

# remove intermediate cols
ogtt.dat$auc_ins_0_30 <- NULL
ogtt.dat$auc_gluc_0_30 <- NULL


################## Other indices ###############################################

# ins fold change
ogtt.dat$ifc_0_120 <- log(ogtt.dat$insulin120 / ogtt.dat$insulin0)

# insulinogenic index
ogtt.dat$insulinogenic_idx <- (ogtt.dat$insulin30 - ogtt.dat$insulin0) / (ogtt.dat$glucose30 - ogtt.dat$glucose0)

# Matsuda index
ogtt.dat$matsuda <- 10000 / sqrt(ogtt.dat$glucose0_mg_dl * ogtt.dat$insulin0_uU_ml * ogtt.dat$gluc_mg_dl_mean * ogtt.dat$ins_uU_ml_mean) 

################## Summarise results ###########################################

## Reshape df
ogtt.indices <- ogtt.dat[,c(1,2,33:49)]
ogtt.indices <- ogtt.indices %>% pivot_longer(cols=3:ncol(ogtt.indices), names_to = "index_name") %>% pivot_wider(names_from = ogtt.session)

colnames(ogtt.indices)[3:4] <- c("OGTT1","OGTT2")
ogtt.indices$change <- (ogtt.indices$OGTT2 - ogtt.indices$OGTT1) / ogtt.indices$OGTT1

## Do ttest
ttest_res <- data.frame(matrix(nrow = 0, ncol = 5))
for (index in unique(ogtt.indices$index_name)) {
  index_data <- ogtt.indices %>% filter(index_name==index) %>% na.omit() %>% filter_all(all_vars(!is.infinite(.)))
  mean1 <- mean(index_data$OGTT1)
  mean2 <- mean(index_data$OGTT2)
  sd1 <- sd(index_data$OGTT1)
  sd2 <- sd(index_data$OGTT2)
  t <- t.test(index_data$OGTT1, index_data$OGTT2, paired = TRUE, alternative = "two.sided")
  beta <- -t$estimate
  beta_rel <- -t$estimate / mean1
  pval <- t$p.value
  t <- data.frame(index, mean1, mean2, sd1, sd2, beta, beta_rel, pval)
  ttest_res <- rbind(ttest_res,t)
}
colnames(ttest_res) <- c("index","ogtt1_mean","ogtt2_mean","ogtt1_sd","ogtt2_sd","beta","beta_rel","pval")
ttest_res <- arrange(ttest_res, pval)
write.csv(ttest_res, "output/02_insulin_indices_ttest.csv", row.names = F)

## Plot significant results as box plot

sig_indices <- ttest_res %>% filter(pval<(0.05/17)) %>% pull(index) 

ogtt.indices_long <- ogtt.indices %>% filter(index_name %in% sig_indices) %>% pivot_longer(cols=3:4,names_to = "ogtt_session", values_to = "value")

facet_titles <- c("ph1_0_60"="Stumvoll 1st phase",
                  "ph2_0_60"="Stumvoll 2nd phase",
                  "ifc_0_120"="Insulin fold change",
                  "hiri"="HIRI",
                  "iauc_gluc_0_120"="Incremental glucose AUC",
                  "matsuda"="Matsuda")
facet_titles <- as.data.frame(tibble(index_name = names(facet_titles),
                                     formatted_name = unlist(facet_titles)))

ogtt.indices_long <- ogtt.indices_long %>% left_join(facet_titles,by="index_name")

ogtt.indices_long$formatted_name <- factor(ogtt.indices_long$formatted_name, levels=c("Stumvoll 2nd phase","Stumvoll 1st phase","HIRI","Insulin fold change","Incremental glucose AUC","Matsuda"))

## Plot 2 handpicked indices for figure 1
svg("graphics/02_matsuda_and_stumvoll.svg", height = 4, width = 8)
ggplot(ogtt.indices_long[ogtt.indices_long$index_name%in%c("ph1_0_60","matsuda"),], aes(x=ogtt_session,y=value)) + 
  geom_boxplot() + geom_line(aes(group=participant), alpha=0.2) +
  theme_minimal() +
  theme(axis.title = element_blank(), text = element_text(size = 20)) +
  facet_wrap(~formatted_name, scales = "free_y", ncol = 2, labeller = labeller(idx = facet_titles))
dev.off()

## Plot the rest
svg("graphics/02_rest_of_indices.svg", height = 4, width = 16)
ggplot(ogtt.indices_long[ogtt.indices_long$index_name%in%c("ph2_0_60","hiri","ifc_0_120","iauc_gluc_0_120"),], aes(x=ogtt_session,y=value)) + 
  geom_boxplot() + geom_line(aes(group=participant), alpha=0.2) +
  theme_minimal() +
  theme(axis.title = element_blank(), text = element_text(size = 18)) +
  facet_wrap(~formatted_name, scales = "free_y", ncol = 4, labeller = labeller(idx = facet_titles))
dev.off()




