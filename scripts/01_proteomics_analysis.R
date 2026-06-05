###################################################
#### Analysis of Olink fasting OGTT study      ####
#### proteomics data                           ####
#### Burulca Uluvar 10.07.2024                 ####
###################################################

## packages needed
require(data.table)
require(lmerTest)
require(doMC)
require(colorspace)
require(basicPlotteR)
require(gprofiler2)
require(ggplot2)
require(dplyr)
require(gridExtra)
require(tidyr)
require(stringr)
require(ggpubr)
require(RColorBrewer)
require(ggrepel)
require(gghighlight)
require(scales)

## Color palette
color_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

########################################
####  Import QCed data set          ####
########################################

ogtt.dat <- read.table(
  "data/ogtt_data.txt",
  header = T,
  sep = "\t"
)

## convert t.point to factor
ogtt.dat$t.factor <- factor(
  ogtt.dat$t.point,
  levels = c("0", "15", "30", "60", "120")
)

## Label mapping for proteins
prot.label <- read.table(
  "data/protein_labels.txt",
  header = T,
  sep = "\t"
)

##  Convert to long format
oid_cols <- grep("^OID", names(ogtt.dat), value = TRUE)

ogtt.dat_long <- ogtt.dat %>%
  pivot_longer(cols = all_of(oid_cols), names_to = "OlinkID", values_to = "NPX")
ogtt.dat_long <- na.omit(ogtt.dat_long)

# Merge with prot.label to get Panel information
ogtt.dat_long <- ogtt.dat_long %>%
  left_join(prot.label[, c("OlinkID", "Panel")], by = c("OlinkID" = "OlinkID"))

# ogtt session as text for plotting
ogtt.dat_long$ogtt.session.txt <- ifelse(
  ogtt.dat_long$ogtt.session == 1,
  "OGTT 1",
  "OGTT 2"
)

##################################################
####  Compute median npx            ####
##################################################

panel_median <- ogtt.dat_long %>%
  group_by(ogtt.session.txt, t.point, Panel, participant) %>%
  summarize(
    median_npx = median(NPX, na.rm = TRUE),
    mean_npx = mean(NPX, na.rm = T)
  )

## Plot median NPX distribution by timepoint
ggplot(
  panel_median,
  aes(x = as.factor(t.point), y = median_npx, group = t.point)
) +
  geom_boxplot() +
  facet_wrap(~ogtt.session.txt) +
  ylab("Median NPX") +
  xlab("Timepoint [Min]") +
  theme_minimal() +
  theme(text = element_text(size = 14))

## Reshape
panel_median <- panel_median %>%
  pivot_wider(
    id_cols = c("ogtt.session.txt", "t.point", "participant"),
    names_from = "Panel",
    values_from = median_npx
  )

##########################################
####  Normalize NPX values          ####
##########################################

## create common scale across proteins, use day 0 minute 0 as a reference for both OGTTs
ogtt.norm <- lapply(prot.label$OlinkID, function(x) {
  ## get the relevant data
  tmp <- ogtt.dat[which(ogtt.dat$ogtt.session == 1 & ogtt.dat$t.point == 0), x]
  return(data.frame(
    OlinkID = x,
    mean.value = mean(tmp, na.rm = TRUE),
    sd.value = sd(tmp, na.rm = TRUE)
  ))
})
## combine again
ogtt.norm <- do.call(rbind, ogtt.norm)

## apply to the data
for (j in 1:nrow(ogtt.norm)) {
  ## N.B.: careful missing values
  ii <- which(!is.na(ogtt.dat[, ogtt.norm$OlinkID[j]]))
  ## scale
  ogtt.dat[ii, ogtt.norm$OlinkID[j]] <- (ogtt.dat[ii, ogtt.norm$OlinkID[j]] -
    ogtt.norm$mean.value[j]) /
    ogtt.norm$sd.value[j]
}

## separate to 2 tables
ogtt.dat.1 <- ogtt.dat %>% filter(ogtt.session == 1)
ogtt.dat.2 <- ogtt.dat %>% filter(ogtt.session == 2)

##########################################
####  Linear model                  ####
##########################################

## Add panel median to the data
ogtt.dat.1 <- merge(
  ogtt.dat.1,
  subset(panel_median, ogtt.session.txt == "OGTT 1"),
  by = c("participant", "t.point")
)
ogtt.dat.2 <- merge(
  ogtt.dat.2,
  subset(panel_median, ogtt.session.txt == "OGTT 2"),
  by = c("participant", "t.point")
)

## convert t.point to factor
ogtt.dat.1$t.factor <- factor(
  ogtt.dat.1$t.point,
  levels = c("0", "15", "30", "60", "120")
)
ogtt.dat.2$t.factor <- factor(
  ogtt.dat.2$t.point,
  levels = c("0", "15", "30", "60", "120")
)

## Load function to run linear model with mixed effects
source("scripts/functions/mixed_effect_regression.R")

## run in parallel
registerDoMC(10)

## adjust for the sample dilution by adding median npx as covariate
res.ogtt.1.linear <- mclapply(
  unique(prot.label$Panel),
  function(x) {
    ## run lmer for panel x, adjusting for median NPX
    tmp <- mixed.anova(
      ogtt.dat.1,
      "t.factor",
      subset(prot.label, Panel == x)$OlinkID,
      paste("+ ", x, "+ (1|participant)")
    )
    ## rename adjustment
    names(tmp) <- gsub(x, "median", names(tmp))
    return(tmp)
  },
  mc.cores = 5
)

res.ogtt.2.linear <- mclapply(
  unique(prot.label$Panel),
  function(x) {
    ## run lmer for panel x, adjusting for median NPX
    tmp <- mixed.anova(
      ogtt.dat.2,
      "t.factor",
      subset(prot.label, Panel == x)$OlinkID,
      paste("+ ", x, "+ (1|participant)")
    )
    ## rename adjustment
    names(tmp) <- gsub(x, "median", names(tmp))
    return(tmp)
  },
  mc.cores = 5
)

## Combine everything
res.ogtt.1.linear <- do.call(rbind, res.ogtt.1.linear)
res.ogtt.2.linear <- do.call(rbind, res.ogtt.2.linear)

## Add protein names
res.ogtt.1.linear <- merge(
  prot.label,
  res.ogtt.1.linear,
  by.x = "OlinkID",
  by.y = "outcome"
)
res.ogtt.2.linear <- merge(
  prot.label,
  res.ogtt.2.linear,
  by.x = "OlinkID",
  by.y = "outcome"
)

## FDR correction
res.ogtt.1.linear$fdr.aov <- p.adjust(
  res.ogtt.1.linear$pval.aov.t.factor,
  method = "BH"
)
res.ogtt.2.linear$fdr.aov <- p.adjust(
  res.ogtt.2.linear$pval.aov.t.factor,
  method = "BH"
)
## sort by fdr
res.ogtt.1.linear <- res.ogtt.1.linear[order(res.ogtt.1.linear$fdr.aov), ]
res.ogtt.2.linear <- res.ogtt.2.linear[order(res.ogtt.2.linear$fdr.aov), ]

##########################################
####  Sensitivity analysis          ####
##########################################

## time point analysis, excluding sample median npx as covariate
res.ogtt.1.sens <- mixed.anova(
  ogtt.dat.1,
  "t.factor",
  prot.label$OlinkID,
  "+ (1|participant)"
)
res.ogtt.2.sens <- mixed.anova(
  ogtt.dat.2,
  "t.factor",
  prot.label$OlinkID,
  "+ (1|participant)"
)

## add protein label
res.ogtt.1.sens <- merge(
  prot.label,
  res.ogtt.1.sens,
  by.x = "OlinkID",
  by.y = "outcome"
)
res.ogtt.2.sens <- merge(
  prot.label,
  res.ogtt.2.sens,
  by.x = "OlinkID",
  by.y = "outcome"
)

## how many FDR significant
res.ogtt.1.sens$fdr.aov <- p.adjust(
  res.ogtt.1.sens$pval.aov.t.factor,
  method = "BH"
)
res.ogtt.2.sens$fdr.aov <- p.adjust(
  res.ogtt.2.sens$pval.aov.t.factor,
  method = "BH"
)

###########################################################
####  Supplementary Figure 2b       ####
###########################################################

##  Plot number of increasing & decreasing proteins

res.ogtt.dfs <- list(res.ogtt.1.linear, res.ogtt.2.linear)
ogtt.dat.dfs <- list(ogtt.dat.1, ogtt.dat.2)

for (i in 1:2) {
  res.ogtt <- res.ogtt.dfs[[i]]
  ogtt.dat <- ogtt.dat.dfs[[i]]

  ## Keep only relevant cols
  selected_cols <- c("OlinkID", "Assay", "UniProt", "fdr.aov")
  pval_cols <- grep("^pval.exp|^beta", names(res.ogtt), value = TRUE)
  selected_cols <- c(selected_cols, pval_cols)
  ogtt_res <- select(res.ogtt, all_of(selected_cols))

  num_sig <- res.ogtt %>%
    filter(fdr.aov < 0.05) %>%
    select(Assay) %>%
    unique() %>%
    count() %>%
    as.numeric()

  ## Significance threshold
  sig_thr <- as.numeric(0.05 / num_sig)

  ## Count how many proteins significantly increasing or decreasing on each t.point
  beta_cols  <- grep("^beta\\.exp", names(ogtt_res), value = TRUE)
  pval_cols2 <- grep("^pval\\.exp", names(ogtt_res), value = TRUE)
  increase <- colSums(
    ogtt_res[, pval_cols2] < sig_thr &
      ogtt_res[, beta_cols]  > 0 &
      ogtt_res[, "fdr.aov"]  < 0.05,
    na.rm = TRUE
  )
  decrease <- -colSums(
    ogtt_res[, pval_cols2] < sig_thr &
      ogtt_res[, beta_cols]  < 0 &
      ogtt_res[, "fdr.aov"]  < 0.05,
    na.rm = TRUE
  )
  res <- t(rbind(increase, decrease)) %>% as.data.frame()
  assign(paste0("res", i), res)
}

inc_dec <- cbind(res1, res2) %>% as.data.frame()
colnames(inc_dec) <- c("increase.1", "decrease.1", "increase.2", "decrease.2")
inc_dec$t.point <- factor(tpoints, levels = c(15, 30, 60, 120))

color_palette <- colorRampPalette(brewer.pal(9, "YlOrRd"))(9)[c(3, 5, 7, 9)]

## ogtt1 plot
p1 <- ggplot(inc_dec, aes(t.point)) +
  geom_bar(
    aes(y = increase.1, fill = t.point),
    stat = "identity",
    position = "dodge",
    alpha = 0.8,
    color = "black"
  ) +
  geom_bar(
    aes(y = decrease.1, fill = t.point),
    stat = "identity",
    position = "dodge",
    alpha = 0.8,
    color = "black"
  ) +
  scale_fill_manual(values = color_palette) +
  ylab("No. of proteins changing during OGTT 1") +
  xlab("Minutes") +
  ylim(-45, 20) +
  geom_hline(yintercept = 0, colour = "black") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 16)
  )

## ogtt2 plot
p2 <- ggplot(inc_dec, aes(t.point)) +
  geom_bar(
    aes(y = increase.2, fill = t.point),
    stat = "identity",
    position = "dodge",
    alpha = 0.8,
    color = "black"
  ) +
  geom_bar(
    aes(y = decrease.2, fill = t.point),
    stat = "identity",
    position = "dodge",
    alpha = 0.8,
    color = "black"
  ) +
  scale_fill_manual(values = color_palette) +
  ylab("No. of proteins changing during OGTT 2") +
  xlab("Minutes") +
  ylim(-45, 20) +
  geom_hline(yintercept = 0, colour = "black") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 16)
  )

grid.arrange(p1, p2, ncol = 2, nrow = 1)

##############################################
####  Supplementary Figure 2a       ####
##############################################

## Plot ogtt1 v 2 qvalues scatterplot

#### Function to select value with max abs value
abs_max <- function(data) {
  tmp <- Filter(is.numeric, data)
  if (inherits(data, "tbl_df")) {
    tmp <- as.matrix(tmp)
  }
  tmp[cbind(1:nrow(tmp), apply(abs(tmp), 1, which.max))]
}

## add as column to result df
beta_cols1 <- res.ogtt.1.linear %>% select(starts_with("beta"))
res.ogtt.1.linear$beta.high <- abs_max(beta_cols1)

beta_cols2 <- res.ogtt.2.linear %>% select(starts_with("beta"))
res.ogtt.2.linear$beta.high <- abs_max(beta_cols2)

#### Function to add first significant tpoint
add_first_sig <- function(df) {
  # Create a new column "first_sig" initialized with NA
  df$first_sig <- NA

  # get col indexes of pval columns
  pval_cols <- grep("pval.exposure", colnames(df))

  # Iterate over each row
  for (i in 1:nrow(df)) {
    # Flag to track if any p-value is below the significance threshold
    pval_below_threshold <- FALSE

    # Iterate over each p-value column
    for (j in pval_cols) {
      # Check if the p-value is below the significance threshold
      if (df[i, j] < significance_threshold) {
        # Set the name of the first significant column in the "first_sig" column
        df[i, "first_sig"] <- names(df)[j]
        # Set the flag to True
        pval_below_threshold <- TRUE
        # Break the loop to record only the first significant column
        break
      }
    }

    # Check if the "first_sig" column is still NA and fdr.aov is less than 0.05
    if (!pval_below_threshold && df[i, "fdr.aov"] < 0.05) {
      # Find the smallest p-value among the pval.exposure columns (use named lookup)
      smallest_pval <- min(df[i, pval_cols], na.rm = TRUE)
      # Find the column name corresponding to the smallest p-value
      col_name <- names(df)[which(df[i, ] == smallest_pval)]
      # Update the "first_sig" column with the column name
      df[i, "first_sig"] <- col_name[1]
    }
  }

  # Remove "pval.exposure" prefix and convert to integer
  df$first_sig <- as.integer(sub("pval.exposure", "", df$first_sig))
  return(df)
}

## Number of significant proteins (computed, used as Bonferroni denominator)
n_sig_1 <- subset(res.ogtt.1.linear, fdr.aov < 0.05) %>% select(Assay) %>% unique() %>% nrow()
n_sig_2 <- subset(res.ogtt.2.linear, fdr.aov < 0.05) %>% select(Assay) %>% unique() %>% nrow()

## Significance threshold for first sig tpoint ogtt 1
significance_threshold <- 0.05 / n_sig_1
res.ogtt.1.linear <- add_first_sig(res.ogtt.1.linear)

## Significance threshold ogtt 2
significance_threshold <- 0.05 / n_sig_2
res.ogtt.2.linear <- add_first_sig(res.ogtt.2.linear)

## Combine ogtt1 and 2
ogtt1 <- res.ogtt.1.linear[, c(
  "OlinkID",
  "Assay",
  "pval.aov.t.factor",
  "fdr.aov",
  "beta.high",
  "first_sig"
)]
ogtt2 <- res.ogtt.2.linear[, c(
  "OlinkID",
  "Assay",
  "pval.aov.t.factor",
  "fdr.aov",
  "beta.high",
  "first_sig"
)]
ogtt.both <- ogtt1 %>%
  left_join(ogtt2, by = c("OlinkID", "Assay"), suffix = c(".1", ".2"))

# point shape circle if betas have different direction & both significant (none)
# triangle if at least one significant increase
# down triangle if at least one significant decrease
ogtt.both$beta_comb <- case_when(
  (ogtt.both$beta.high.1 * ogtt.both$beta.high.2 < 0 &
    ogtt.both$fdr.aov.1 < .05 &
    ogtt.both$fdr.aov.2 < .05) ~ 23,
  (ogtt.both$beta.high.1 > 0 & ogtt.both$fdr.aov.1 < .05) |
    (ogtt.both$beta.high.2 > 0 & ogtt.both$fdr.aov.2 < .05) ~ 24,
  (ogtt.both$beta.high.1 < 0 & ogtt.both$fdr.aov.1 < .05) |
    (ogtt.both$beta.high.2 < 0 & ogtt.both$fdr.aov.2 < .05) ~ 25,
  TRUE ~ 21
)

# pick the beta with the higher significance to plot
ogtt.both <- ogtt.both %>%
  mutate(beta_sig = ifelse(fdr.aov.1 < fdr.aov.2, beta.high.1, beta.high.2))

# pick the first sig that is earlier
ogtt.both$first.sig_comb <- pmin(
  ogtt.both$first_sig.1,
  ogtt.both$first_sig.2,
  na.rm = T
)

# pick the smaller fdr.aov
ogtt.both$fdr.aov.smaller <- pmin(
  ogtt.both$fdr.aov.1,
  ogtt.both$fdr.aov.2,
  na.rm = T
)

# deduplicate multiple olink ids per assay (e.g. IL6)
ogtt.both <- ogtt.both %>%
  group_by(Assay) %>%
  slice(which.min(fdr.aov.smaller)) %>%
  ungroup()

top.pvals <- ogtt.both %>%
  filter(
    fdr.aov.1 < 1e-3 | fdr.aov.2 < 1e-3 | (fdr.aov.1 < .05 & fdr.aov.2 < .05)
  )
mypal <- colorRampPalette(brewer.pal(9, "YlOrRd"))

# Display GCG as GLP1
ogtt.both$Assay <- gsub("GCG", "GLP-1", ogtt.both$Assay)

# Plot scatterplot of qvals
ggplot(
  ogtt.both,
  aes(
    x = -log10(fdr.aov.1),
    y = -log10(fdr.aov.2),
    fill = as.factor(first.sig_comb),
    size = abs(beta_sig)
  )
) +
  geom_point(
    shape = ogtt.both$beta_comb,
    alpha = ifelse(
      ogtt.both$fdr.aov.1 < .05 | ogtt.both$fdr.aov.2 < .05,
      1,
      0.2
    )
  ) +
  scale_fill_manual(
    name = "Time of first change",
    values = mypal(9)[c(3, 5, 7, 9)],
    breaks = c("15", "30", "60", "120")
  ) +
  labs(x = "OGTT 1 -log10(qval)", y = "OGTT 2 -log10(qval)") +
  geom_text_repel(data = top.pvals, aes(label = Assay), size = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "orange") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "orange") +
  scale_size_continuous(name = "Max. change") + #, labels = c("1", "4"), breaks=c(1, 4)) +
  guides(
    fill = guide_legend(override.aes = list(shape = 22, size = 6)),
    size = guide_legend(override.aes = list(shape = c(25, 24), size = c(1, 4)))
  ) +
  theme_minimal() +
  theme(legend.position = c(0.9, 0.25), axis.text = element_text(size = 10))

#######################################################
####  Figure 2a: betas scatterplot  ####
#######################################################

top.betas <- ogtt.both %>%
  filter(fdr.aov.smaller < 0.05 & (abs(beta.high.1) > 1 | abs(beta.high.2) > 1))

ggplot(
  ogtt.both,
  aes(x = beta.high.1, y = beta.high.2, fill = as.factor(first.sig_comb))
) +
  geom_point(
    shape = 21,
    size = 2,
    alpha = ifelse(ogtt.both$fdr.aov.smaller < 0.05, 1, 0.05)
  ) +
  labs(x = "Max. change during OGTT1", y = "Max. change during OGTT2") +
  geom_text_repel(data = top.betas, aes(label = Assay), size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  scale_fill_manual(
    name = "Time of first change",
    values = mypal(9)[c(3, 5, 7, 9)],
    breaks = c("15", "30", "60", "120")
  ) +
  guides(fill = guide_legend(override.aes = list(size = 6))) +
  theme_classic() +
  theme(panel.grid = element_blank(), legend.position = c(0.85, 0.15), axis.text = element_text(size = 10))

##########################################
####  Differential response analysis####
##########################################

## import new function
source("scripts/functions/mixed_effect_regression_interaction.R")

registerDoMC(10)

## linear model with interaction
res.ogtt.diff <- mclapply(
  unique(prot.label$Panel),
  function(x) {
    ## run lmer for panel x, adjusting for median NPX
    tmp <- mixed.regression.v2(
      rbind(ogtt.dat.1, ogtt.dat.2),
      subset(prot.label, Panel == x)$OlinkID,
      paste("ogtt.session*t.factor + ", x, " + (1|participant)"),
      "inter",
      log.y = F,
      scale = F
    )
    ## rename adjustment
    names(tmp) <- gsub(x, "median", names(tmp))
    return(tmp)
  },
  mc.cores = 5
)

res.ogtt.diff <- do.call(rbind, res.ogtt.diff)

## add protein label
res.ogtt.diff <- merge(
  prot.label,
  res.ogtt.diff,
  by.x = "OlinkID",
  by.y = "feat"
)

##########################################################
####  Sensitivity analyses          ####
##########################################################

####  Sensitivity analysis 1        ####

##(no median NPX adjustment)

res.ogtt.diff.sens <- mixed.regression.v2(
  rbind(ogtt.dat.1, ogtt.dat.2),
  prot.label$OlinkID,
  "ogtt.session*t.factor + (1|participant)",
  "inter",
  log.y = F,
  scale = F
)
## add protein label
res.ogtt.diff.sens <- merge(
  prot.label,
  res.ogtt.diff.sens,
  by.x = "OlinkID",
  by.y = "feat"
)
res.ogtt.diff.sens <- res.ogtt.diff.sens %>%
  arrange(fdr.anova_ogtt.session.t.factor)
diff.prots.sens <- subset(
  res.ogtt.diff.sens,
  fdr.anova_ogtt.session.t.factor < .2
)$Assay

####  Sensitivity analysis 2        ####

## (normalization separately on both days)

## Read the non-normalised data again
ogtt.dat.cp <- read.table(
  "data/ogtt_data.txt",
  header = T,
  sep = "\t"
)
ogtt.dat.cp$t.factor <- factor(
  ogtt.dat.cp$t.point,
  levels = c("0", "15", "30", "60", "120")
)

## Normalise OGTT 1
ogtt.dat.cp.1 <- subset(ogtt.dat.cp, ogtt.session == 1)

## create common scale across proteins, use day 0 minute 0 as a reference for both OGTTs
ogtt.norm.1 <- lapply(prot.label$OlinkID, function(x) {
  ## get the relevant data
  tmp <- ogtt.dat.cp.1[which(ogtt.dat.cp.1$t.point == 0), x]
  return(data.frame(
    OlinkID = x,
    mean.value = mean(tmp, na.rm = TRUE),
    sd.value = sd(tmp, na.rm = TRUE)
  ))
})
## combine again
ogtt.norm.1 <- do.call(rbind, ogtt.norm.1)

## apply to the data
for (j in 1:nrow(ogtt.norm.1)) {
  ## N.B.: careful missing values
  ii <- which(!is.na(ogtt.dat.cp.1[, ogtt.norm.1$OlinkID[j]]))
  ## scale
  ogtt.dat.cp.1[ii, ogtt.norm.1$OlinkID[j]] <- (ogtt.dat.cp.1[
    ii,
    ogtt.norm.1$OlinkID[j]
  ] -
    ogtt.norm.1$mean.value[j]) /
    ogtt.norm.1$sd.value[j]
}
## check whether it worked
sd(ogtt.dat.cp.1[which(ogtt.dat.cp.1$t.point == 0), "OID30150"], na.rm = T)

## Normalise OGTT 2
ogtt.dat.cp.2 <- subset(ogtt.dat.cp, ogtt.session == 2)

## create common scale across proteins, use day 0 minute 0 as a reference for both OGTTs
ogtt.norm.2 <- lapply(prot.label$OlinkID, function(x) {
  ## get the relevant data
  tmp <- ogtt.dat.cp.2[which(ogtt.dat.cp.2$t.point == 0), x]
  return(data.frame(
    OlinkID = x,
    mean.value = mean(tmp, na.rm = TRUE),
    sd.value = sd(tmp, na.rm = TRUE)
  ))
})
## combine again
ogtt.norm.2 <- do.call(rbind, ogtt.norm.2)

## apply to the data
for (j in 1:nrow(ogtt.norm.2)) {
  ## N.B.: careful missing values
  ii <- which(!is.na(ogtt.dat.cp.2[, ogtt.norm.2$OlinkID[j]]))
  ## scale
  ogtt.dat.cp.2[ii, ogtt.norm.2$OlinkID[j]] <- (ogtt.dat.cp.2[
    ii,
    ogtt.norm.2$OlinkID[j]
  ] -
    ogtt.norm.2$mean.value[j]) /
    ogtt.norm.2$sd.value[j]
}

res.ogtt.diff.sens2 <- mixed.regression.v2(
  rbind(ogtt.dat.cp.1, ogtt.dat.cp.2),
  prot.label$OlinkID,
  "ogtt.session*t.factor + (1|participant)",
  "inter",
  log.y = F,
  scale = F
)

## add protein label
res.ogtt.diff.sens2 <- merge(
  prot.label,
  res.ogtt.diff.sens2,
  by.x = "OlinkID",
  by.y = "feat"
)
res.ogtt.diff.sens2 <- res.ogtt.diff.sens2 %>%
  arrange(fdr.anova_ogtt.session.t.factor)

diff.prots.sens2 <- subset(
  res.ogtt.diff.sens2,
  fdr.anova_ogtt.session.t.factor < .2
)$Assay

###>> MLN, CDSN, WARS, PYY are consistently significant

##########################################
####  Group proteins by consistency ####
##########################################

## Add column indicating whether protein changes consistently or only in one
ogtt.both$consistency <- case_when(
  ogtt.both$Assay %in%
    c("WARS", "MLN", "CDSN", "PYY") ~ "Differential response",
  ogtt.both$fdr.aov.1 < 0.05 &
    ogtt.both$fdr.aov.2 < 0.05 ~ "Consistent in both",
  ogtt.both$fdr.aov.1 < 0.05 &
    ogtt.both$pval.aov.t.factor.2 > 0.05 ~ "Only OGTT 1",
  ogtt.both$fdr.aov.2 < 0.05 &
    ogtt.both$pval.aov.t.factor.1 > 0.05 ~ "Only OGTT 2",
  TRUE ~ "Not significant"
)

##########################################
####  Plotting functions            ####
##########################################

## Prepare df for plotting
ogtt.dat.norm_long <- rbind(ogtt.dat.1, ogtt.dat.2) %>%
  pivot_longer(
    cols = all_of(oid_cols),
    names_to = "OlinkID",
    values_to = "NPX"
  ) %>%
  left_join(prot.label, by = c("OlinkID" = "OlinkID"))

ogtt.dat.norm_long <- ogtt.dat.norm_long %>%
  select(
    "participant",
    "t.point",
    "ogtt.session",
    "t.factor",
    "Assay",
    "NPX"
  ) %>% # OlinkID
  group_by(t.point, ogtt.session, t.factor, Assay) %>%
  summarise(
    mean = mean(NPX, na.rm = T),
    se = sd(NPX, na.rm = TRUE) / sqrt(n()),
    sqrt_n = sqrt(n())
  )

ogtt.dat.norm_long <- na.omit(ogtt.dat.norm_long)
ogtt.dat.norm_long$ogtt.session <- as.factor(ogtt.dat.norm_long$ogtt.session)

## Plotting function
plot_prots <- function(prots, ncol) {
  ## Order by lowest fdr
  subset_df <- ogtt.both[ogtt.both$Assay %in% prots, ] %>%
    arrange(fdr.aov.smaller)
  subset_df$Assay <- gsub("GCG", "GLP-1", subset_df$Assay)
  prots <- subset_df$Assay

  ogtt.dat.norm_long$Assay <- gsub("GCG", "GLP-1", ogtt.dat.norm_long$Assay)

  p <- ogtt.dat.norm_long %>%
    filter(Assay %in% prots) %>%
    mutate(Assay = factor(Assay, levels = prots)) %>%
    ggplot(aes(x = t.point, y = mean, color = ogtt.session)) +
    geom_line(size = 1, position = position_dodge(width = 3)) +
    geom_point(size = 2, position = position_dodge(width = 3)) +
    geom_errorbar(
      aes(
        x = t.point,
        ymin = mean - se * qt(p = 0.975, df = 10),
        ymax = mean + se * qt(p = 0.975, df = 10)
      ),
      position = position_dodge(width = 3),
      width = 0
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    theme_classic() +
    labs(x = "Time [min]", y = "Standardized NPX") +
    theme(
      axis.title.x = element_text(hjust = 0.5),
      axis.title.y = element_text(hjust = 0.5),
      legend.position = "bottom",
      legend.title = element_text(size = 16),
      text = element_text(size = 16, hjust = 0),
      plot.margin = unit(c(0, 1, 0, 1), "cm"),
      strip.background = element_blank()
    ) +
    scale_x_continuous(breaks = c(0, 15, 30, 60, 120)) +
    scale_color_manual(
      name = "OGTT session",
      values = c("#6699CC", "#CC6677")
    ) +
    facet_wrap(~Assay, ncol = ncol) +
    guides(colour = guide_legend(nrow = 1))

  return(p)
}

## Plotting function (individual level)
plot_prots_ind <- function(prots) {
  ## Order by lowest fdr
  subset_df <- ogtt.both[ogtt.both$Assay %in% prots, ] %>%
    arrange(fdr.aov.smaller)
  subset_df$Assay <- gsub("GCG", "GLP-1", subset_df$Assay)
  prots <- subset_df$OlinkID

  ogtt.dat.ind <- rbind(ogtt.dat.1, ogtt.dat.2) %>%
    pivot_longer(
      cols = all_of(oid_cols),
      names_to = "OlinkID",
      values_to = "NPX"
    ) %>%
    left_join(prot.label, by = c("OlinkID" = "OlinkID"))

  ogtt.dat.ind$Assay <- gsub("GCG", "GLP-1", ogtt.dat.ind$Assay)

  ## Calculate mean NPX for each time point and assay
  mean_data <- ogtt.dat.ind %>%
    filter(OlinkID %in% prots) %>%
    group_by(t.point, Assay, ogtt.session) %>%
    summarise(mean_NPX = mean(NPX, na.rm = TRUE), .groups = 'drop')

  ## Plot individual lines in grey
  p <- ogtt.dat.ind %>%
    filter(OlinkID %in% prots) %>%
    ggplot(aes(x = t.point, y = NPX, group = participant)) +
    geom_line(color = "grey", size = 0.3) + # Individual lines
    geom_line(
      data = mean_data,
      aes(x = t.point, y = mean_NPX, group = 1),
      color = "black",
      size = 1.2
    ) + # Mean line
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      text = element_text(size = 12, hjust = 0),
      plot.margin = unit(c(0, 1, 0, 1), "cm"),
      panel.spacing = unit(1, "lines")
    ) +
    scale_x_continuous(breaks = c(0, 15, 30, 60, 120)) +
    scale_y_continuous(labels = scales::number_format(accuracy = 1)) + # Round y-axis
    facet_grid(cols = vars(Assay), rows = vars(ogtt.session))

  return(p)
}

##############################################################
####  Figures & Tables              ####
##############################################################

## Proteins that change early and strongly 
early_strong <- ogtt.both %>%
  filter(
    (abs(beta.high.1) > 1.5 & fdr.aov.1 < 0.05 & first_sig.1 %in% c(15, 30)) |
      (abs(beta.high.2) > 1.5 & fdr.aov.2 < 0.05 & first_sig.2 %in% c(15, 30))
  )
early_strong_prots <- early_strong %>% filter(Assay != "PYY") %>% pull(Assay) # exclude PYY (differential)

plot_prots(early_strong_prots, 2) + theme(legend.position = c(0.7, 0.1))

## Same proteins, individual level plots
plot_prots_ind(early_strong_prots)

## proteins coexpressed with anxa10
plot_prots(c("ANXA10","VSIG2","CA9","TFF2","PLA2G10","AGR2"), 3) + theme(legend.position = c(0.9,0.97))

## differentially responding proteins
plot_prots(c("MLN", "WARS", "PYY", "CDSN"), 2) + theme(legend.position = "none")

## Write linear model results to file
write.table(
  res.ogtt.1.linear,
  "output/01_res.ogtt.1.linear.txt",
  sep = "\t",
  row.names = F
)
write.table(
  res.ogtt.2.linear,
  "output/01_res.ogtt.2.linear.txt",
  sep = "\t",
  row.names = F
)

## Supplementary Table: Linear model results differential effects
write.table(
  res.ogtt.diff,
  "output/01_res.ogtt.diff.txt",
  sep = "\t",
  row.names = F
)

#############################################
####  Sex differential response analysis ####
#############################################

# add sex to the df
participant <- read_excel("data/participant_data.xlsx")
participant <- participant %>% select(participant, sex, age) %>% unique()

## OGTT 1
ogtt.dat.1 <- ogtt.dat.1 %>% left_join(participant, by = "participant")
res.ogtt1.diff.sex <- mclapply(
  unique(prot.label$Panel),
  function(x) {
    ## run lmer for panel x, adjusting for median NPX
    tmp <- mixed.regression.v2(
      ogtt.dat.1,
      subset(prot.label, Panel == x)$OlinkID,
      paste("sex*t.factor + ", x, " + (1|participant)"),
      "inter",
      log.y = F,
      scale = F
    )
    ## rename adjustment
    names(tmp) <- gsub(x, "median", names(tmp))
    return(tmp)
  },
  mc.cores = 5
)

res.ogtt1.diff.sex <- do.call(rbind, res.ogtt1.diff.sex)
## add protein label
res.ogtt1.diff.sex <- merge(
  prot.label,
  res.ogtt1.diff.sex,
  by.x = "OlinkID",
  by.y = "feat"
)

## which of them was significant overall
res.ogtt.1.linear <- fread("output/01_res.ogtt.1.linear.txt")
sig1 <- res.ogtt.1.linear %>% filter(fdr.aov < 0.05) %>% pull(Assay)
res.ogtt1.diff.sex %>%
  filter(fdr.anova_sex.t.factor < 0.2 & Assay %in% sig1) %>%
  select(Assay, fdr.anova_sex.t.factor) %>%
  arrange(fdr.anova_sex.t.factor)

## OGTT 2
ogtt.dat.2 <- ogtt.dat.2 %>% left_join(participant, by = "participant")
res.ogtt2.diff.sex <- mclapply(
  unique(prot.label$Panel),
  function(x) {
    ## run lmer for panel x, adjusting for median NPX
    tmp <- mixed.regression.v2(
      ogtt.dat.2,
      subset(prot.label, Panel == x)$OlinkID,
      paste("sex*t.factor + ", x, " + (1|participant)"),
      "inter",
      log.y = F,
      scale = F
    )
    ## rename adjustment
    names(tmp) <- gsub(x, "median", names(tmp))
    return(tmp)
  },
  mc.cores = 5
)

res.ogtt2.diff.sex <- do.call(rbind, res.ogtt2.diff.sex)

## add protein label
res.ogtt2.diff.sex <- merge(
  prot.label,
  res.ogtt2.diff.sex,
  by.x = "OlinkID",
  by.y = "feat"
)
res.ogtt2.diff.sex %>%
  filter(fdr.anova_sex.t.factor < 0.2) %>%
  select(Assay, fdr.anova_sex.t.factor) %>%
  arrange(fdr.anova_sex.t.factor)
##>>> no resutls

res.ogtt2.diff.sex %>%
  filter(Assay == "LMOD1") %>%
  select(OlinkID, fdr.anova_sex.t.factor)

##########################################
####  HPA tissue specificity          ####
##########################################

## Join HPA data for tissue specificity, data dowloaded from Human Protein Atlas
hpa_data <- read.table(
  "data/hpa_data.tsv",
  header = T
)
hpa_data <- hpa_data[, c(
  "OlinkID",
  "Assay",
  "max_ntpm_tissue",
  "max_ntpm_celltype",
  "RNA.tissue.specificity",
  "RNA.single.cell.type.specificity",
  "Subcellular.location",
  "Secretome.location"
)]
ogtt.both_hpa <- ogtt.both %>%
  left_join(hpa_data, by = "OlinkID", multiple = "all") #--> 2940 proteins

## clean up
ogtt.both_hpa$max_ntpm_tissue <- gsub(" 1", "", ogtt.both_hpa$max_ntpm_tissue)
ogtt.both_hpa$max_ntpm_tissue <- gsub(
  "choroid plexus",
  "brain",
  ogtt.both_hpa$max_ntpm_tissue
)

## consider only tissue-specific proteins for the enrichment analysis (i.e. those with "Tissue enriched" or "Tissue enhanced" specificity in HPA)
ogtt.both_hpa$max_ntpm_tissue <- ifelse(
  ogtt.both_hpa$RNA.tissue.specificity %in%
    c("Tissue enriched", "Tissue enhanced"),
  ogtt.both_hpa$max_ntpm_tissue,
  NA
)

## Check enrichment of tissue-specific proteins

## Compute no. of sig/nonsig counts dynamically for Fisher's exact test
n_sig_hpa_1    <- ogtt.both_hpa %>% filter(fdr.aov.1 < .05) %>% nrow()
n_nonsig_hpa_1 <- ogtt.both_hpa %>% filter(!fdr.aov.1 < .05) %>% nrow()
n_sig_hpa_2    <- ogtt.both_hpa %>% filter(fdr.aov.2 < .05) %>% nrow()
n_nonsig_hpa_2 <- ogtt.both_hpa %>% filter(!fdr.aov.2 < .05) %>% nrow()

#### Get numbers for contingency table
ogtt.by.tissue.1 <- ogtt.both_hpa %>%
  count(max_ntpm_tissue, sig = fdr.aov.1 < .05) %>%
  pivot_wider(names_from = sig, values_from = n) %>%
  filter(!is.na(max_ntpm_tissue))
ogtt.by.tissue.1[is.na(ogtt.by.tissue.1)] <- 0

ogtt.by.tissue.2 <- ogtt.both_hpa %>%
  count(max_ntpm_tissue, sig = fdr.aov.2 < .05) %>%
  pivot_wider(names_from = sig, values_from = n) %>%
  filter(!is.na(max_ntpm_tissue))
ogtt.by.tissue.2[is.na(ogtt.by.tissue.2)] <- 0

fisher_test_tissue <- function(ogtt_by_tissue, sig, nonsig) {
  ## Loop through tissues
  fisher_df = data.frame(matrix(vector(), 0, 3))
  for (t in ogtt_by_tissue$max_ntpm_tissue) {
    tissue <- ogtt_by_tissue %>% filter(max_ntpm_tissue == t)
    contingency_table <- matrix(
      c(
        tissue$`TRUE`,
        sig - tissue$`TRUE`,
        tissue$`FALSE`,
        nonsig - tissue$`FALSE`
      ),
      nrow = 2
    )
    f <- fisher.test(contingency_table, alternative = "greater")
    p <- f$p.value
    or <- f$estimate
    fisher_df <- rbind(fisher_df, c(t, p, or))
  }
  colnames(fisher_df) <- c("tissue", "pval", "OR")
  fisher_df$pval <- as.numeric(fisher_df$pval)
  fisher_df$OR <- as.numeric(fisher_df$OR)
  fisher_df$fdr <- p.adjust(fisher_df$pval, method = "BH")
  fisher_df$tissue <- str_to_title(fisher_df$tissue)

  return(fisher_df)
}

## apply to dfs and plot
fisher_df_1 <- fisher_test_tissue(ogtt.by.tissue.1, n_sig_hpa_1, n_nonsig_hpa_1)
fisher_df_2 <- fisher_test_tissue(ogtt.by.tissue.2, n_sig_hpa_2, n_nonsig_hpa_2)
fisher_df_1$ogtt.session = "1"
fisher_df_2$ogtt.session = "2"
fisher_df <- rbind(fisher_df_1, fisher_df_2)

## Plot enrichment results
fisher_df %>%
  ggplot(aes(
    x = -log10(fdr),
    y = reorder(tissue, -fdr),
    size = OR,
    col = ogtt.session
  )) +
  geom_point(alpha = 0.7) +
  geom_vline(
    xintercept = -log10(.05),
    linetype = "dashed",
    colour = "#CC6677"
  ) +
  labs(x = "-log10(p-adj)") +
  scale_color_manual(name = "OGTT session", values = c("#6699CC", "#CC6677")) +
  theme_minimal() +
  theme(legend.position = c(0.9, 0.5), axis.title.y = element_blank())
