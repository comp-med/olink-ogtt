###################################################
#### Analysis of Olink fasting OGTT study      ####
#### proteomics data                           ####
#### Burulca Uluvar 10.07.2024                 ####
###################################################

rm(list=ls())
options(stringsAsFactors = F)

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

load("RData/01.RData")

## Color palette
#color_palette <- colorRampPalette(brewer.pal(11, "Spectral"))(11)[-c(6, 7)]
color_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

########################################
####      import QCed data set      ####
########################################

## Import overall data
npx.data         <- read.table("data/QCed.Olink.Explore.3072.Fasting.study.20230111.txt", header=T)

## Label for proteins
prot.label       <- read.table("data/Protein.Label.Olink.Explore.3072.Fasting.study.20220620.txt", header=T)

## Separate fasting and ogtt data 
## Fasting data
fast.dat <- subset(npx.data, study.type == "Fasting") 
## OGTT data
ogtt.dat <- subset(npx.data, study.type == "OGTT")

## Set sample that did not pass QC to NA
neurII <- prot.label %>% filter(Panel=="Neurology_II") %>% pull(OlinkID) 
ogtt.dat[(ogtt.dat$participant == "FP27" & ogtt.dat$ogtt.session==1 & ogtt.dat$t.point==0), neurII] <- NA

##########################################
####  Compute median npx of samples   ####
##########################################

oid_cols <- grep("^OID", names(ogtt.dat), value = TRUE)

ogtt.dat_long <- ogtt.dat %>%
  pivot_longer(cols = all_of(oid_cols), names_to = "OlinkID", values_to = "NPX")
ogtt.dat_long <- na.omit(ogtt.dat_long)

# Merge 'ogtt.dat_long' with 'prot.label' to get Panel information
ogtt.dat_long <- ogtt.dat_long %>%
  left_join(prot.label[,c("OlinkID","Panel")], by = c("OlinkID" = "OlinkID"))
ogtt.dat_long$ogtt.session.txt <- ifelse(ogtt.dat_long$ogtt.session==1,"OGTT 1", "OGTT 2")

# Calculate median NPX per sample & panel
panel_median <- ogtt.dat_long %>%
  group_by(ogtt.session.txt, t.point, Panel, participant) %>%
  summarize(median_npx = median(NPX, na.rm = TRUE),
            mean_npx = mean(NPX, na.rm =T))

# Plot to pdf
pdf("graphics/median_npx.pdf", width = 10, height = 6)
ggplot(panel_median, aes(x=as.factor(t.point), y=median_npx, group=t.point)) + geom_boxplot() +
  facet_wrap(~ogtt.session.txt) +
  ylab("Median NPX") + xlab("Timepoint [Min]") +
  theme_minimal() + 
  theme(text = element_text(size = 14))
dev.off()

panel_median <- panel_median %>% pivot_wider(id_cols = c("ogtt.session.txt", "t.point", "participant"), names_from="Panel", values_from = median_npx)

##########################################
####       create normalization       ####
##########################################

## create common scale across proteins, use day 0 minute 0 as a reference for both OGTTs
ogtt.norm <- lapply(prot.label$OlinkID, function(x){
  ## get the relevant data
  tmp <- ogtt.dat[which(ogtt.dat$ogtt.session==1 & ogtt.dat$t.point == 0), x]
  return(data.frame(OlinkID = x, mean.value = mean(tmp, na.rm = TRUE), sd.value = sd(tmp, na.rm = TRUE)))
})
## combine again
ogtt.norm <- do.call(rbind, ogtt.norm)

## apply to the data
for(j in 1:nrow(ogtt.norm)){
  ## N.B.: careful missing values
  ii                                     <- which(!is.na(ogtt.dat[, ogtt.norm$OlinkID[j]]))
  ## scale
  ogtt.dat[ii, ogtt.norm$OlinkID[j]] <- (ogtt.dat[ii, ogtt.norm$OlinkID[j]] - ogtt.norm$mean.value[j])/ogtt.norm$sd.value[j]
}
## check whether it worked
sd(ogtt.dat[which(ogtt.dat$ogtt.session==1 & ogtt.dat$t.point == 0), "OID30150"], na.rm=T)

write.table(ogtt.dat, "output/01_ogtt.dat.norm.txt", sep="\t", row.names = F)

## separate to 2 tables
ogtt.dat.1 <- ogtt.dat %>% filter(ogtt.session ==1)
ogtt.dat.2 <- ogtt.dat %>% filter(ogtt.session ==2)

##########################################
####     Linear model                 ####
##########################################

## Add panel median to the data
ogtt.dat.1 <- merge(ogtt.dat.1, subset(panel_median, ogtt.session.txt=="OGTT 1"), by=c("participant","t.point"))
ogtt.dat.2 <- merge(ogtt.dat.2, subset(panel_median, ogtt.session.txt=="OGTT 2"), by=c("participant","t.point"))

## Load function to run linear model with mixed effects
source("scripts/functions/mixed_effect_regression.R")

## run in parallel
registerDoMC(10)

## adjust for the sample "dilution" by adding median npx as covariate
res.ogtt.1.linear <- mclapply(unique(prot.label$Panel), function(x) {
  ## run lmer for panel x, adjusting for median NPX
  tmp <- mixed.anova(ogtt.dat.1, "t.factor", subset(prot.label, Panel==x)$OlinkID, paste("+ " , x, "+ (1|participant)"))
  ## rename adjustment
  names(tmp) <- gsub(x, "median", names(tmp))
  return(tmp)
}, mc.cores=5)

res.ogtt.2.linear <- mclapply(unique(prot.label$Panel), function(x) {
  ## run lmer for panel x, adjusting for median NPX
  tmp <- mixed.anova(ogtt.dat.2, "t.factor", subset(prot.label, Panel==x)$OlinkID, paste("+ " , x, "+ (1|participant)"))
  ## rename adjustment
  names(tmp) <- gsub(x, "median", names(tmp))
  return(tmp)
}, mc.cores=5)

## Combine everything
res.ogtt.1.linear <- do.call(rbind, res.ogtt.1.linear)
res.ogtt.2.linear <- do.call(rbind, res.ogtt.2.linear)

## Add protein names
res.ogtt.1.linear <- merge(prot.label, res.ogtt.1.linear, by.x="OlinkID", by.y="outcome")
res.ogtt.2.linear <- merge(prot.label, res.ogtt.2.linear, by.x="OlinkID", by.y="outcome")

## how many FDR significant overall
res.ogtt.1.linear$fdr.aov <- p.adjust(res.ogtt.1.linear$pval.aov.t.factor, method = "BH")
res.ogtt.2.linear$fdr.aov <- p.adjust(res.ogtt.2.linear$pval.aov.t.factor, method = "BH")
nrow(subset(res.ogtt.1.linear, fdr.aov < .05)) ## n=56
nrow(subset(res.ogtt.2.linear, fdr.aov < .05)) ## n=71

## how many unique
sig1 <- subset(res.ogtt.1.linear, fdr.aov < .05) %>% select(Assay) %>% unique() ## n=52
sig2 <- subset(res.ogtt.2.linear, fdr.aov < .05) %>% select(Assay) %>% unique() ## n=71

## sort by fdr
res.ogtt.1.linear <- res.ogtt.1.linear[order(res.ogtt.1.linear$fdr.aov), ]
res.ogtt.2.linear <- res.ogtt.2.linear[order(res.ogtt.2.linear$fdr.aov), ]

## write to file
write.table(res.ogtt.1.linear, "output/01_res.ogtt.1.linear.txt", sep="\t", row.names = F)
write.table(res.ogtt.2.linear, "output/01_res.ogtt.2.linear.txt", sep="\t", row.names = F)

## read again
res.ogtt.1.linear <- read.table("output/01_res.ogtt.1.linear.txt", sep="\t", header = T)
res.ogtt.2.linear <- read.table("output/01_res.ogtt.2.linear.txt", sep="\t", header = T)


##########################################
####      Sensitivity analysis        ####
##########################################

## convert t.point to factor 
ogtt.dat.1$t.factor       <- factor(ogtt.dat.1$t.point, levels=c("0", "15", "30", "60", "120")) 
ogtt.dat.2$t.factor       <- factor(ogtt.dat.2$t.point, levels=c("0", "15", "30", "60", "120")) 

## time point analysis, excluding sample median npx as covariate
res.ogtt.1.sens         <- mixed.anova(ogtt.dat.1, "t.factor", prot.label$OlinkID, "+ (1|participant)")
res.ogtt.2.sens         <- mixed.anova(ogtt.dat.2, "t.factor", prot.label$OlinkID, "+ (1|participant)")

## add protein label
res.ogtt.1.sens         <- merge(prot.label, res.ogtt.1.sens, by.x="OlinkID", by.y="outcome")
res.ogtt.2.sens         <- merge(prot.label, res.ogtt.2.sens, by.x="OlinkID", by.y="outcome")

## how many FDR significant overall
res.ogtt.1.sens$fdr.aov <- p.adjust(res.ogtt.1.sens$pval.aov.t.factor, method = "BH")
res.ogtt.2.sens$fdr.aov <- p.adjust(res.ogtt.2.sens$pval.aov.t.factor, method = "BH")
nrow(subset(res.ogtt.1.sens, fdr.aov < .05)) ## n=51
nrow(subset(res.ogtt.2.sens, fdr.aov < .05)) ## n=208

res.ogtt.1.sens <- res.ogtt.1.sens[order(res.ogtt.1.sens$fdr.aov), ]
res.ogtt.2.sens <- res.ogtt.2.sens[order(res.ogtt.2.sens$fdr.aov), ]

## Write to file
write.table(res.ogtt.1.sens, "output/01_res.ogtt.1.sens.txt", sep="\t", row.names = F)
write.table(res.ogtt.2.sens, "output/01_res.ogtt.2.sens.txt", sep="\t", row.names = F)

## Compare pvals from both methods
compare.1 <- res.ogtt.1.sens[,c("OlinkID","Assay","Panel","fdr.aov")] %>% left_join(res.ogtt.1.linear[,c("OlinkID","fdr.aov")], by="OlinkID")
compare.2 <- res.ogtt.2.sens[,c("OlinkID","Assay","Panel","fdr.aov")] %>% left_join(res.ogtt.2.linear[,c("OlinkID","fdr.aov")], by="OlinkID")

## Plot with ggplot
top.1 <- compare.1[compare.1$fdr.aov.x<0.05|compare.1$fdr.aov.y<0.05, ]
p1 <- ggplot(compare.1, aes(x=-log10(fdr.aov.x), y=-log10(fdr.aov.y), color=Panel)) + geom_point() +
  xlab("-log10(qval) without adjustment") + ylab("-log10(qval) adjusted for median NPX") +
  theme_minimal() + ggtitle("OGTT 1") +
  geom_text_repel(data=top.1 ,aes(label = Assay), hjust = 1.2, vjust = 0.5, size=2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "#117733", linewidth=1) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#117733", linewidth=1) +
  theme(legend.position = "top") +
  scale_color_manual(values = color_palette)

top.2 <- compare.2[compare.2$fdr.aov.x<0.05|compare.2$fdr.aov.y<0.05, ]
p2 <- ggplot(compare.2, aes(x=-log10(fdr.aov.x), y=-log10(fdr.aov.y), color=Panel)) + geom_point() +
  xlab("-log10(qval) without adjustment") + ylab("-log10(qval) adjusted for median NPX") +
  theme_minimal() + ggtitle("OGTT 2") +
  geom_text_repel(data = top.2, aes(label = Assay), hjust = 1.2, vjust = 0.5, size=2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "#117733", linewidth=1) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#117733", linewidth=1) +
  theme(legend.position = "none") +
  scale_color_manual(values = color_palette)

p1b <- ggplot(compare.1, aes(x=-log10(fdr.aov.x), y=-log10(fdr.aov.y), color=Panel)) + geom_point() +
  xlab("-log10(qval) without adjustment") + ylab("-log10(qval) adjusted for median NPX") +
  theme_minimal() + ggtitle("OGTT 1 (zoom)") + xlim(0,2.5) + ylim(0,2.5) +
  geom_text_repel(data=top.1 ,aes(label = Assay), hjust = 1.2, vjust = 0.5, size=2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "#117733", linewidth=1) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#117733", linewidth=1) +
  theme(legend.position = "none") +
  scale_color_manual(values = color_palette)

p2b <- ggplot(compare.2, aes(x=-log10(fdr.aov.x), y=-log10(fdr.aov.y), color=Panel)) + geom_point() +
  xlab("-log10(qval) without adjustment") + ylab("-log10(qval) adjusted for median NPX") +
  theme_minimal() + ggtitle("OGTT 2 (zoom)") + xlim(0,2.5) + ylim(0,2.5) +
  geom_text_repel(data = top.2, aes(label = Assay), hjust = 1.2, vjust = 0.5, size=2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "#117733", linewidth=1) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#117733", linewidth=1) +
  theme(legend.position = "none") +
  scale_color_manual(values = color_palette)

pdf("graphics/qvals_comparison.pdf", width = 12, height = 12)
grid.arrange(p1,p2,p1b,p2b, nrow=2, ncol=2)
dev.off()


##########################################
####  Plot inc. and decr. proteins    ####
##########################################

res.ogtt.dfs <- list(res.ogtt.1.linear, res.ogtt.2.linear)
ogtt.dat.dfs <- list(ogtt.dat.1, ogtt.dat.2)

for (i in 1:2){
  res.ogtt <- res.ogtt.dfs[[i]]
  ogtt.dat <- ogtt.dat.dfs[[i]]
  
  #### Keep only relevant cols
  selected_cols <- c("OlinkID", "Assay", "UniProt", "fdr.aov")
  pval_cols <- grep("^pval.exp|^beta", names(res.ogtt), value = TRUE)
  selected_cols <- c(selected_cols, pval_cols)
  ogtt_res <- select(res.ogtt, all_of(selected_cols))
  
  num_sig <- res.ogtt %>% filter(fdr.aov<0.05) %>% select(Assay) %>% unique() %>% count() %>% as.numeric()
  
  # Significance threshold
  sig_thr <- as.numeric(0.05/num_sig)
  
  # Count how many proteins significantly increasing or decreasing on each t.point
  increase <- colSums(ogtt_res[, c(6, 8, 10, 12)] < sig_thr & ogtt_res[, c(5, 7, 9, 11)] > 0 & ogtt_res[,"fdr.aov"] < 0.05, na.rm = TRUE)
  decrease <- -colSums(ogtt_res[, c(6, 8, 10, 12)] < sig_thr & ogtt_res[, c(5, 7, 9, 11)] < 0 & ogtt_res[,"fdr.aov"] < 0.05, na.rm = TRUE)
  res <- t(rbind(increase,decrease)) %>% as.data.frame()
  assign(paste0("res", i), res)
}

inc_dec <- cbind(res1,res2) %>% as.data.frame()
colnames(inc_dec) <- c("increase.1","decrease.1","increase.2","decrease.2")
tpoints <- c(15,30,60,120)
inc_dec$t.point <- factor(tpoints, levels = tpoints)

color_palette <- colorRampPalette( brewer.pal( 9 , "YlOrRd" ) )(9)[c(3,5,7,9)]

## ogtt1
p1 <- ggplot(inc_dec, aes(t.point)) + 
  geom_bar(aes(y = increase.1, fill = t.point), stat = "identity", position = "dodge", alpha=0.8, color="black") +
  geom_bar(aes(y = decrease.1, fill = t.point), stat = "identity", position = "dodge", alpha=0.8, color="black") +
  scale_fill_manual(values=color_palette) +
  ylab("No. of proteins changing during OGTT 1") + xlab("Minutes") +
  ylim(-45,20) +
  geom_hline(yintercept = 0,colour = "black") +
  theme_minimal() + theme(legend.position = "none",
                          axis.text = element_text(size=12),
                          axis.title = element_text(size=16))

## ogtt2
p2 <- ggplot(inc_dec, aes(t.point)) + 
  geom_bar(aes(y = increase.2, fill = t.point), stat = "identity", position = "dodge", alpha=0.8, color="black") +
  geom_bar(aes(y = decrease.2, fill = t.point), stat = "identity", position = "dodge", alpha=0.8, color="black") +
  scale_fill_manual(values=color_palette) +
  ylab("No. of proteins changing during OGTT 2") + xlab("Minutes") +
  ylim(-45,20) +
  geom_hline(yintercept = 0,colour = "black") +
  theme_minimal() + theme(legend.position = "none",
                          axis.text = element_text(size=12),
                          axis.title = element_text(size=16))

pdf("graphics/01_no_of_inc_decr_prots.pdf", width=6, height =6)
grid.arrange(p1,p2, ncol = 2, nrow = 1)
dev.off()

##########################################
####  Plot ogtt1 v 2 scatterplot      ####
##########################################

# Function to select value with max abs value
abs_max <- function(data) {
  tmp <- Filter(is.numeric, data)
  if (inherits(data, "tbl_df")) {
    tmp <- as.matrix(tmp)
  }
  tmp[cbind(1:nrow(tmp), apply(abs(tmp), 1, which.max))]
}

# add as column to result df
beta_cols1 <- res.ogtt.1.linear %>% select(starts_with("beta"))
res.ogtt.1.linear$beta.high <- abs_max(beta_cols1)
res.ogtt.1.linear$beta.high.tpoint <- abs_max_tpoint(beta_cols1)

beta_cols2 <- res.ogtt.2.linear %>% select(starts_with("beta"))
res.ogtt.2.linear$beta.high <- abs_max(beta_cols2)

## proteins with >1.5 s.d. change
high.beta.1 <- subset(res.ogtt.1.linear, fdr.aov < .05 & abs(beta.high) >= 1.5) %>% pull(Assay) %>% unique()  ## "IGFBP1"    "ANXA10"    "VSIG2"     "IL6"       "GIP"       "CRIP2"     "GCG"       "ADCYAP1R1"
high.beta.2 <- subset(res.ogtt.2.linear, fdr.aov < .05 & abs(beta.high) >= 1.5) %>% pull(Assay) %>% unique()  ## "ANXA10"  "VSIG2"   "WARS"    "PYY"     "SPOCK1"  "GCG"     "CRIP2"   "GIP"     "SFRP1"   "PTPRN2"  "CLEC12A"

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
      # Find the smallest p-value among the pval.exposure columns
      smallest_pval <- min(df[i, c(13,16,19,22)], na.rm = TRUE)
      # Find the column name corresponding to the smallest p-value
      col_name <- names(df)[which(df[i, ] == smallest_pval)]
      # Update the "first_sig" column with the column name
      df[i, "first_sig"] <- col_name
    }
  }
  
  # Remove "pval.exposure" prefix and convert to integer
  df$first_sig <- as.integer(sub("pval.exposure", "", df$first_sig))
  return(df)
}

# Significance threshold ogtt 1
significance_threshold <- 0.05/52
res.ogtt.1.linear <- add_first_sig(res.ogtt.1.linear)

# Significance threshold ogtt 2
significance_threshold <- 0.05/71
res.ogtt.2.linear <- add_first_sig(res.ogtt.2.linear)


# Combine ogtt1 and 2
ogtt1 <-res.ogtt.1.linear[,c("OlinkID","Assay","pval.aov.t.factor","fdr.aov","beta.high","first_sig")]
ogtt2 <-res.ogtt.2.linear[,c("OlinkID","Assay","pval.aov.t.factor","fdr.aov","beta.high","first_sig")]
ogtt.both <- ogtt1 %>% left_join(ogtt2, by=c("OlinkID","Assay"), suffix=c(".1",".2"))

## nothing is in different direction & significant
ogtt.both[ogtt.both$beta.high.1 * ogtt.both$beta.high.2 < 0 &
            ogtt.both$fdr.aov.1 < 0.05 &
            ogtt.both$fdr.aov.2 < 0.05, ]

## how many are significant in one OGTT and not even nominal significance in the other
ogtt1_only <- ogtt.both %>% filter(fdr.aov.1<0.05 & pval.aov.t.factor.2>0.05)
ogtt1_only %>% filter(abs(beta.high.1)>1.5)

ogtt2_only <- ogtt.both %>% filter(fdr.aov.2<0.05 & pval.aov.t.factor.1>0.05)
ogtt2_only %>% filter(abs(beta.high.2)>1.5)

ogtt_consistent <- ogtt.both %>% filter(fdr.aov.1<0.05&pval.aov.t.factor.2<0.05 | fdr.aov.2<0.05&pval.aov.t.factor.1<0.05)
length(unique(ogtt_consistent$Assay))

# point shape circle if betas have different direction & both significant
# triangle if at least one significant increase
# down triangle if at least one significant decr
ogtt.both$beta_comb <- case_when( (ogtt.both$beta.high.1 * ogtt.both$beta.high.2 < 0 & ogtt.both$fdr.aov.1<.05 & ogtt.both$fdr.aov.2<.05 )  ~ 23,
                                  (ogtt.both$beta.high.1>0 & ogtt.both$fdr.aov.1<.05) | (ogtt.both$beta.high.2>0 & ogtt.both$fdr.aov.2<.05) ~ 24, 
                                  (ogtt.both$beta.high.1<0 & ogtt.both$fdr.aov.1<.05) | (ogtt.both$beta.high.2<0 & ogtt.both$fdr.aov.2<.05) ~ 25,
                                  TRUE ~ 21)
# pick the beta with the higher significance to plot
ogtt.both <- ogtt.both %>%
  mutate(beta_sig = ifelse(fdr.aov.1 < fdr.aov.2, beta.high.1, beta.high.2))

# pick the first sig that is earlier
ogtt.both$first.sig_comb <- pmin(ogtt.both$first_sig.1, ogtt.both$first_sig.2, na.rm = T)

# pick the smaller fdr.aov
ogtt.both$fdr.aov.smaller <- pmin(ogtt.both$fdr.aov.1, ogtt.both$fdr.aov.2, na.rm = T)

# deduplicate multiple olink ids per assay (e.g. IL6)
ogtt.both <- ogtt.both %>% group_by(Assay) %>% slice(which.min(fdr.aov.smaller)) %>% ungroup()

top.pvals <- ogtt.both %>% filter(fdr.aov.1<1e-3 | fdr.aov.2<1e-3 | (fdr.aov.1<.05 & fdr.aov.2<.05))
mypal <- colorRampPalette( brewer.pal( 9 , "YlOrRd" ) )

# Plot - I couldnt figure out the size labels
svg("graphics/01_ogtt1_vs_2_scatter.svg", width=6, height=6)

ggplot(ogtt.both, aes(x = -log10(fdr.aov.1), y = -log10(fdr.aov.2),
                      fill = as.factor(first.sig_comb), size=abs(beta_sig))) +
  geom_point(shape=ogtt.both$beta_comb,
             alpha=ifelse(ogtt.both$fdr.aov.1<.05|ogtt.both$fdr.aov.2<.05, 1, 0.2)
             ) +
  scale_fill_manual(name="Time of first change", values = mypal(9)[c(3,5,7,9)], breaks=c("15","30","60","120")) +
  labs(x = "OGTT 1 -log10(qval)", y = "OGTT 2 -log10(qval)") + 
  geom_text_repel(data = top.pvals, aes(label = Assay), size=3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "orange") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "orange") +
  scale_size_continuous(name = "Max. change") + #, labels = c("1", "4"), breaks=c(1, 4)) +
  guides(fill = guide_legend(override.aes=list(shape=22, size=6)),
         size = guide_legend(override.aes=list(shape=c(25,24), size=c(1,4)) )
         )+
  theme_minimal() +
  theme(legend.position = c(0.9,0.25),
        axis.text = element_text(size=10))

dev.off()

# Plot betas instead
top.betas <- ogtt.both %>% filter( fdr.aov.smaller<0.05 & (abs(beta.high.1)>1 | abs(beta.high.2)>1)) 

svg("graphics/01_ogtt1_vs_2_scatter_betas.svg", width=6, height=6)

ggplot(ogtt.both, aes(x = beta.high.1, y = beta.high.2)) +
  geom_point(shape=21,
             fill=ifelse(ogtt.both$fdr.aov.smaller<0.05,"darkorange","white"),
             alpha=ifelse(ogtt.both$fdr.aov.smaller<0.05, 1, 0.05)
  ) +
  labs(x = "OGTT 1 max. change", y = "OGTT 2 max. change") + 
  geom_text_repel(data = top.betas, aes(label = Assay), size=3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "orange") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "orange") +
  scale_size_continuous(name = "Max. change") + #, labels = c("1", "4"), breaks=c(1, 4)) +
  theme_minimal() +
  theme(legend.position = c(0.9,0.25),
        axis.text = element_text(size=10))

dev.off()


## Get number of proteins in each quadrant
ogtt.both %>% filter(fdr.aov.1<.05 & fdr.aov.2<.05) %>% pull(Assay) %>% unique() %>% length() # 17
ogtt.both %>% filter(fdr.aov.1<.05 & fdr.aov.2>.05) %>% pull(Assay) %>% unique() %>% length() # 35
ogtt.both %>% filter(fdr.aov.1>.05 & fdr.aov.2<.05) %>% pull(Assay) %>% unique() %>% length() # 54
ogtt.both %>% filter(fdr.aov.1<.05) %>% pull(Assay) %>% unique() %>% length() #52


## write to file (to add beta.high and first_sig)
write.table(res.ogtt.1.linear, "output/01_res.ogtt.1.linear.txt", sep="\t", row.names = F)
write.table(res.ogtt.2.linear, "output/01_res.ogtt.2.linear.txt", sep="\t", row.names = F)
write.table(ogtt.both, "output/01_ogtt.both.summary.txt", sep="\t", row.names = F)

## read from file
res.ogtt.1.linear <- read.table("output/01_res.ogtt.1.linear.txt", sep="\t", header = T)
res.ogtt.2.linear <- read.table("output/01_res.ogtt.2.linear.txt", sep="\t", header = T)
ogtt.both <- read.table("output/01_ogtt.both.summary.txt", sep="\t", header = T)

##########################################
####  Plot ogtt1 and 2 volcano plots  ####
##########################################

############# OGTT1 volcano plot 

color_palette <- colorRampPalette( brewer.pal( 9 , "YlOrRd" ) )(9)[c(4,5,7,9)]

top.pvals.1 <- ogtt.both %>% filter(fdr.aov.1<1e-3 | abs(beta.high.1)>1.5 & fdr.aov.1<.05)
v1 <- ggplot(ogtt.both, aes(x=beta.high.1, y=-log10(fdr.aov.1), col = as.factor(first_sig.1))) + 
  geom_point() +
  scale_color_manual(name="Time of first change", values = color_palette, breaks=c("15","30","60","120")) +
  labs(x = "OGTT 1 max. change", y = "OGTT 1 -log10(qval)") + 
  geom_text_repel(data = top.pvals.1, aes(label = Assay), vjust = -0.5, size=3) +
  theme_minimal() + theme(legend.position="bottom",
                          text = element_text(size=14))

############# OGTT2 volcano plot

top.pvals.2 <- ogtt.both %>% filter(fdr.aov.2<1e-3 | abs(beta.high.2)>1.5 & fdr.aov.2<.05)
v2 <- ggplot(ogtt.both, aes(x=beta.high.2, y=-log10(fdr.aov.2), col = as.factor(first_sig.2))) + 
  geom_point() +
  scale_color_manual(name="Time of first change", values = color_palette, breaks=c("15","30","60","120")) +
  labs(x = "OGTT 2 max. change", y = "OGTT 2 -log10(qval)") + 
  geom_text_repel(data = top.pvals.2, aes(label = Assay), vjust = -0.5, size=3) +
  theme_minimal() + theme(legend.position="bottom",
                          text = element_text(size=14))

pdf("graphics/ogtt1_and_2_volcano.pdf", width=12, height=8)
grid.arrange(v1,v2, ncol=2)
dev.off()

##########################################
####  Plot fasting t0 vs t7 volcano   ####
##########################################

## Read fasting results
fasting_res_raw <- read.table("/sc-projects/sc-proj-computational-medicine/people/Maik/25_Olink_fasting_study/04_analysis/data/Results.Fasting.linear.simple.20230130.txt",
                              header=TRUE)
fasting_res <- fasting_res_raw[,c("OlinkID","Assay","UniProt","fdr.aov","beta.exposure7","pval.exposure7")]

## deduplicate multiple olink ids per assay (e.g. IL6)
fasting_res <- fasting_res %>% group_by(Assay) %>% slice(which.min(fdr.aov)) %>% ungroup()

## Join ogtt to fasting results
fasting_ogtt <- fasting_res %>% left_join(ogtt.both,by="Assay")

## Assign colors
fasting_ogtt$point_color <- case_when(fasting_ogtt$fdr.aov.1<.05 & fasting_ogtt$fdr.aov.2<.05 ~ "#117733",
                                      fasting_ogtt$fdr.aov.1<.05 ~ "#6699CC",
                                      fasting_ogtt$fdr.aov.2<.05 ~ "#CC6677",
                                      TRUE ~ "#A5A5A5")
fasting_ogtt$point_color <- factor(fasting_ogtt$point_color, levels=c("#117733","#6699CC","#CC6677","#A5A5A5"))

## Pick proteins to label
top.pvals <- fasting_ogtt %>% filter( (fdr.aov.1 < .05 | fdr.aov.2 < .05) & (fdr.aov<1e-8 | abs(beta.exposure7) > 1) )

## Plot
pdf("graphics/fasting_t7_vs_t0_volcano.pdf", width=8, height=6)
ggplot(data = fasting_ogtt, aes(x = beta.exposure7, y = -log10(fdr.aov), 
                                fill = point_color,
                                shape = beta_comb,
                                size = ifelse(point_color == "#A5A5A5", 0.5, abs(beta_sig)))) +
  geom_point(data = subset(fasting_ogtt, point_color == "#A5A5A5"), stroke=0, alpha=0.5) +   # Nonsig points
  geom_point(data = subset(fasting_ogtt, point_color != "#A5A5A5")) +  # Non-grey points visually on top
  geom_text_repel(data = top.pvals, aes(x = beta.exposure7, y = -log10(fdr.aov), label = Assay, color = point_color),
                  size = 2.5,
                  fontface = "bold"
                  ) +
  scale_fill_manual(name="OGTT response",
                    labels = c("Significant in both OGTTs", "Only in OGTT1", "Only in OGTT2", "Not significant"),
                    values = c("#117733","#6699CC","#CC6677","#A5A5A5"), 
                    breaks=c("#117733","#6699CC","#CC6677","#A5A5A5")) +
  scale_shape_identity() +
  scale_color_identity() +
  guides(fill = guide_legend(override.aes=list(shape=22, size=6)),
         size = guide_legend(title="Max. change during OGTT", override.aes=list(shape=c(25,24), size=c(1,3))) ) +
  labs(x = "Fasting day 7 vs day 0 beta", y = "-log10(qval)") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#888888") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#888888") +
  theme_minimal()

dev.off()

##########################################
####  HPA tissue specificity          ####
##########################################

## Join HPA data for tissue specificity
hpa_data <- read.table("/sc-projects/sc-proj-computational-medicine/data/07_public/01_Human_Protein_Atlas/data/olink_hpa_data_2023-07-24_processed.tsv", header = T)
hpa_data <- hpa_data[,c("OlinkID","Assay","max_ntpm_tissue","max_ntpm_celltype", "RNA.tissue.specificity","RNA.single.cell.type.specificity","Subcellular.location","Secretome.location")]
ogtt.both_hpa <- ogtt.both %>% left_join(hpa_data, by="OlinkID", multiple="all") #--> 2940 proteins

## clean up
ogtt.both_hpa$max_ntpm_tissue <- gsub(" 1", "", ogtt.both_hpa$max_ntpm_tissue)
ogtt.both_hpa$max_ntpm_tissue <- gsub("choroid plexus", "brain", ogtt.both_hpa$max_ntpm_tissue)
ogtt.both_hpa$max_ntpm_tissue <- ifelse(ogtt.both_hpa$RNA.tissue.specificity %in% c("Tissue enriched","Tissue enhanced"), ogtt.both_hpa$max_ntpm_tissue, NA)
ogtt.both_hpa$max_ntpm_celltype <- ifelse(ogtt.both_hpa$RNA.single.cell.type.specificity %in% c("Cell type enriched","Cell type enhanced"), ogtt.both_hpa$max_ntpm_celltype, NA)

ogtt.both_hpa$genes_group <- case_when( ogtt.both_hpa$Assay.x %in% c("WARS","MLN","CDSN","PYY") ~ "Differential response",
                                        ogtt.both_hpa$fdr.aov.1<0.05 & ogtt.both_hpa$fdr.aov.2<0.05 ~ "Consistent in both",
                                        ogtt.both_hpa$fdr.aov.1<0.05 & ogtt.both_hpa$pval.aov.t.factor.2>0.05 ~ "Only OGTT 1",
                                        ogtt.both_hpa$fdr.aov.2<0.05 & ogtt.both_hpa$pval.aov.t.factor.1>0.05 ~ "Only OGTT 2",
                                        TRUE ~ "Not significant")

ogtt.both_hpa_long <- ogtt.both_hpa %>%
  # Select columns that end with .1 or .2 specifically
  select(Assay.x, max_ntpm_tissue, matches("pval.aov.t.factor.\\d|fdr.aov.\\d|beta.high.\\d")) %>%
  pivot_longer(
    cols = matches("pval.aov.t.factor.\\d|fdr.aov.\\d|beta.high.\\d"),
    names_to = c(".value", "ogtt.session"),
    names_pattern = "(.*)\\.(\\d)"
  )

ogtt.both_hpa_long$sig <- ogtt.both_hpa_long$fdr.aov<0.05

protein_counts <- ogtt.both_hpa_long %>%
  filter(!is.na(max_ntpm_tissue) & sig) %>%
  distinct(Assay.x, max_ntpm_tissue, ogtt.session) %>%
  count(max_ntpm_tissue, ogtt.session) %>%
  complete(max_ntpm_tissue, ogtt.session = c("1", "2"), fill = list(n = 0)) # Add all combinations with n = 0 if missing

# Calculate total proteins per tissue
total_counts <- protein_counts %>%
  group_by(max_ntpm_tissue) %>%
  summarize(total_n = sum(n), .groups = "drop") %>%
  arrange(desc(total_n))

# Reorder max_ntpm_tissue factor levels based on total counts
protein_counts <- protein_counts %>%
  mutate(max_ntpm_tissue = factor(max_ntpm_tissue, levels = total_counts$max_ntpm_tissue))

pdf("graphics/01_no_of_proteins_by_tissue.pdf", width = 6, height = 3)
ggplot(protein_counts, aes(x = max_ntpm_tissue, y = n, fill = ogtt.session)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(y = "Number of Proteins", fill = "OGTT Session" ) +
  theme_minimal() +
  scale_fill_manual(values=c("#6699CC", "#CC6677")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", axis.title.x = element_blank())
dev.off()

######### check enrichment - tissue ############################################

## How many sig and nonsig
ogtt.both_hpa %>% count(sig=fdr.aov.1<.05) #2888 nonsig, 52 sig ogtt1
ogtt.both_hpa %>% count(sig=fdr.aov.2<.05) #2869 nonsig, 71 sig ogtt2

#### Get numbers for contingency table
ogtt.by.tissue.1 <- ogtt.both_hpa %>% count(max_ntpm_tissue, sig=fdr.aov.1<.05) %>%
  pivot_wider(names_from = sig, values_from = n) %>% filter(!is.na(max_ntpm_tissue))
ogtt.by.tissue.1[is.na(ogtt.by.tissue.1)] <- 0

ogtt.by.tissue.2 <- ogtt.both_hpa %>% count(max_ntpm_tissue, sig=fdr.aov.2<.05) %>%
  pivot_wider(names_from = sig, values_from = n) %>% filter(!is.na(max_ntpm_tissue))
ogtt.by.tissue.2[is.na(ogtt.by.tissue.2)] <- 0

fisher_test_tissue <- function(ogtt_by_tissue, sig, nonsig) {
  
  ## Loop through tissues
  fisher_df = data.frame(matrix(vector(), 0, 3 ))
  for (t in ogtt_by_tissue$max_ntpm_tissue) {
    tissue <- ogtt_by_tissue %>% filter(max_ntpm_tissue==t)
    contingency_table <- matrix(c(tissue$`TRUE`, sig-tissue$`TRUE`, tissue$`FALSE`, nonsig-tissue$`FALSE`), nrow = 2)
    f <- fisher.test(contingency_table, alternative = "greater")
    p <- f$p.value
    or <- f$estimate
    fisher_df <- rbind(fisher_df, c(t,p,or))
  }
  colnames(fisher_df) <- c("tissue","pval","OR")
  fisher_df$pval <- as.numeric(fisher_df$pval)
  fisher_df$OR <- as.numeric(fisher_df$OR)
  fisher_df$fdr <- p.adjust(fisher_df$pval, method = "BH")
  fisher_df$tissue <- str_to_title(fisher_df$tissue)
  
  return(fisher_df)
}

## apply to dfs and plot
fisher_df_1 <- fisher_test_tissue(ogtt.by.tissue.1, 52, 2888)
fisher_df_2 <- fisher_test_tissue(ogtt.by.tissue.2, 71, 2869)
fisher_df_1$ogtt.session = "1"
fisher_df_2$ogtt.session = "2"
fisher_df <- rbind(fisher_df_1, fisher_df_2)

enriched <- ogtt.both_hpa %>% filter(max_ntpm_tissue %in% c("stomach","pituitary gland","intestine") & fdr.aov.smaller<.05)

## Plot enrichment to pdf
pdf("graphics/both_ogtt_tissue_enrichmt.pdf", width=6, height=6)
fisher_df %>% ggplot(aes(x=-log10(fdr),y=reorder(tissue,-fdr),size=OR, col=ogtt.session)) + 
  geom_point(alpha=0.7) + geom_vline(xintercept = -log10(.05), linetype = "dashed", colour = "#CC6677") +
  labs(x="-log10(p-adj)") +
  scale_color_manual(name="OGTT session",values=c("#6699CC", "#CC6677")) +
  theme_minimal() +
  theme(legend.position = c(0.9,0.5),
        axis.title.y = element_blank())
dev.off()

######### check enrichment - celltype ##########################################

#### Get numbers for contingency table
ogtt.by.ct.1 <- ogtt.both_hpa %>% count(max_ntpm_celltype, sig=fdr.aov.1<.05) %>%
  pivot_wider(names_from = sig, values_from = n) %>% filter(!is.na(max_ntpm_celltype))
ogtt.by.ct.1[is.na(ogtt.by.ct.1)] <- 0

ogtt.by.ct.2 <- ogtt.both_hpa %>% count(max_ntpm_celltype, sig=fdr.aov.2<.05) %>%
  pivot_wider(names_from = sig, values_from = n) %>% filter(!is.na(max_ntpm_celltype))
ogtt.by.ct.2[is.na(ogtt.by.ct.2)] <- 0

fisher_test_ct <- function(ogtt_by_ct, sig, nonsig) {
  
  ## Loop through celltypes
  fisher_df = data.frame(matrix(vector(), 0, 3 ))
  for (c in ogtt.by.celltype$max_ntpm_celltype) {
    celltype <- ogtt.by.celltype %>% filter(max_ntpm_celltype==c)
    contingency_table <- matrix(c(celltype$`TRUE`, sig-celltype$`TRUE`, celltype$`FALSE`, nonsig-celltype$`FALSE`), nrow = 2)
    f <- fisher.test(contingency_table, alternative = "greater")
    p <- f$p.value
    or <- f$estimate
    fisher_df <- rbind(fisher_df, c(c,p,or))
  }
  colnames(fisher_df) <- c("celltype","pval","OR")
  fisher_df$pval <- as.numeric(fisher_df$pval)
  fisher_df$OR <- as.numeric(fisher_df$OR)
  fisher_df$fdr <- p.adjust(fisher_df$pval, method = "BH")
  return(fisher_df)
}  

## apply to dfs and plot
fisher_df_1 <- fisher_test_ct(ogtt.by.ct.1, 52, 2888)
fisher_df_2 <- fisher_test_ct(ogtt.by.ct.2, 71, 2869)
fisher_df_1$ogtt.session = "1"
fisher_df_2$ogtt.session = "2"
fisher_df <- rbind(fisher_df_1, fisher_df_2)

enriched <- ogtt.both_hpa %>% filter(max_ntpm_celltype %in% c("Glandular and luminal cells","Gastric mucus-secreting cells","Enteroendocrine cells", "skeletal myocytes","Suprabasal keratinocytes") & fdr.aov.smaller<.05)

## Plot enrichment
pdf("graphics/both_ogtt_celltype_enrichmt.pdf", width=8, height=8)
fisher_df %>% ggplot(aes(x=-log10(fdr),y=reorder(celltype,-fdr),size=OR, col=ogtt.session)) + 
  geom_point(alpha=0.7) + geom_vline(xintercept = -log10(.05), linetype = "dashed", colour = "#CC6677") +
  labs(x="-log10(p-adj)") +
  scale_color_manual(name="OGTT session",values=c("#6699CC", "#CC6677")) +
  theme_minimal() +
  theme(legend.position = c(0.9,0.5),
        axis.title.y = element_blank(),
        text = element_text(size=10))
dev.off()


##########################################
####  differential response analysis  ####
##########################################

## import new function
source("../scripts/mixed_effect_regression_version2.R")

res.ogtt.diff <- mclapply(unique(prot.label$Panel), function(x) {
  ## run lmer for panel x, adjusting for median NPX
  tmp <- mixed.regression.v2(rbind(ogtt.dat.1, ogtt.dat.2),
                             subset(prot.label, Panel==x)$OlinkID,
                             paste("ogtt.session*t.factor + " , x, " + (1|participant)"),
                             "inter", log.y=F, scale=F)
  ## rename adjustment
  names(tmp) <- gsub(x, "median", names(tmp))
  return(tmp)
}, mc.cores=5)

res.ogtt.diff <- do.call(rbind, res.ogtt.diff)
head(res.ogtt.diff)
## add protein label
res.ogtt.diff     <- merge(prot.label, res.ogtt.diff, by.x="OlinkID", by.y="feat")
## how many findings
nrow(subset(res.ogtt.diff, fdr.anova_ogtt.session.t.factor < .2))
## (before dil.corr. = 6)
## n = 7

res.ogtt.diff <- res.ogtt.diff %>% arrange(fdr.anova_ogtt.session.t.factor)

subset(res.ogtt.diff, fdr.anova_ogtt.session.t.factor < .2)$Assay
# "MLN"    "CDSN"   "ING1"   "TIGAR"  "SCARB2" "SFRP1"  "WARS"   "PYY"   
diff.prots <- subset(res.ogtt.diff, fdr.anova_ogtt.session.t.factor < .2)$Assay

## Write to file
write.table(res.ogtt.diff, "output/01_res.ogtt.diff.txt", sep="\t", row.names = F)

##########################################
####    Plot selected proteins        ####
##########################################

## Prepare df for plotting
ogtt.dat.norm_long <- rbind(ogtt.dat.1, ogtt.dat.2) %>% 
  pivot_longer(cols = all_of(oid_cols), names_to = "OlinkID", values_to = "NPX") %>%
  left_join(prot.label, by = c("OlinkID" = "OlinkID"))

ogtt.dat.norm_long <- ogtt.dat.norm_long %>%
  select("participant","t.point","ogtt.session","t.factor","Assay","NPX") %>% # OlinkID
  group_by(t.point, ogtt.session, t.factor, Assay) %>%
  summarise(mean=mean(NPX,na.rm=T),
            se = sd(NPX, na.rm = TRUE) / sqrt(n()),
            sqrt_n = sqrt(n()))

ogtt.dat.norm_long <- na.omit(ogtt.dat.norm_long)
ogtt.dat.norm_long$ogtt.session <- as.factor(ogtt.dat.norm_long$ogtt.session)

write.table(ogtt.dat.norm_long, "output/01_ogtt.dat.norm_long.txt", sep="\t", row.names = F)

################################################################################

## Plotting function
plot_prots <- function(prots, ncol) {
  
  ## Order by lowest fdr
  subset_df <- ogtt.both[ogtt.both$Assay %in% prots, ] %>% arrange(fdr.aov.smaller)
  subset_df$Assay <- gsub("GCG","GLP-1",subset_df$Assay)
  prots <- subset_df$Assay
  
  ogtt.dat.norm_long$Assay <- gsub("GCG","GLP-1",ogtt.dat.norm_long$Assay)
  
  p <- ogtt.dat.norm_long %>% filter(Assay %in% prots) %>% 
    mutate(Assay = factor(Assay, levels = prots)) %>%
    ggplot(aes(x=t.point,y=mean, color=ogtt.session)) +
    geom_line(size=1, position=position_dodge(width=3)) + 
    geom_point(size=2, position=position_dodge(width=3)) +
    geom_errorbar(aes(x = t.point,
                      ymin=mean-se*qt(p = 0.975, df = 10),
                      ymax=mean+se*qt(p = 0.975, df = 10)),
                  position=position_dodge(width=3), width=0) +
    theme_minimal() + 
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "bottom",
          legend.title = element_text(size=16),
          text = element_text(size = 16, hjust = 0),
          plot.margin = unit(c(0,1,0,1), "cm")) +
    scale_x_continuous(breaks = c(0, 15, 30, 60, 120)) +
    scale_color_manual(name="OGTT session",values=c("#6699CC", "#CC6677")) +
    facet_wrap(~Assay, ncol=ncol) + #, scales = "free_y"
    guides(colour = guide_legend(nrow = 1))
  
  return(p)
  
}

## Plotting function (individual level)
plot_prots_ind <- function(prots) {
  
  ## Order by lowest fdr
  subset_df <- ogtt.both[ogtt.both$Assay %in% prots, ] %>% arrange(fdr.aov.smaller)
  subset_df$Assay <- gsub("GCG", "GLP-1", subset_df$Assay)
  prots <- subset_df$OlinkID
  
  ogtt.dat.ind <- rbind(ogtt.dat.1, ogtt.dat.2) %>% 
    pivot_longer(cols = all_of(oid_cols), names_to = "OlinkID", values_to = "NPX") %>%
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
    geom_line(color = "grey", size = 0.3) +  # Individual lines
    geom_line(data = mean_data, aes(x = t.point, y = mean_NPX, group = 1), 
              color = "black", size = 1.2) +  # Mean line
    theme_minimal() + 
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          text = element_text(size = 12, hjust = 0),
          plot.margin = unit(c(0, 1, 0, 1), "cm"),
          panel.spacing = unit(1, "lines")) +
    scale_x_continuous(breaks = c(0, 15, 30, 60, 120)) +
    scale_y_continuous(labels = scales::number_format(accuracy = 1)) +  # Round y-axis
    facet_grid(cols = vars(Assay), rows = vars(ogtt.session))
  
  return(p)
}


################################################################################

## No. of proteins that change strongly on ogtt 1 and 2 separately
ogtt.both %>% filter(abs(beta.high.1)>1.5 & fdr.aov.1<0.05) %>% pull(Assay) %>% unique() %>% length()
ogtt.both %>% filter(abs(beta.high.2)>1.5 & fdr.aov.2<0.05) %>% pull(Assay) %>% unique() %>% length()

## proteins that change significantly on both
ogtt.both %>% filter(fdr.aov.1<0.05 & fdr.aov.2<0.05) %>% pull(Assay)

## No. of proteins that change strongly on both ogtts
ogtt.both %>% filter(abs(beta.high.1)>1.5 & fdr.aov.1<0.05 & abs(beta.high.2)>1.5 & fdr.aov.2<0.05) %>% pull(Assay) # %>% unique() %>% length()

## Proteins that change strongly only in ogtt 1
ogtt.both %>% filter(abs(beta.high.1)>1.5 & fdr.aov.1<0.05 & pval.aov.t.factor.2>0.05) #%>% pull(Assay) %>% unique() %>% length()

## Proteins that change strongly only in ogtt 2
ogtt.both %>% filter(abs(beta.high.2)>1.5 & fdr.aov.2<0.05 & pval.aov.t.factor.1>0.05) #%>% pull(Assay) %>% unique() %>% length()

################################################################################

#### Plot sets of proteins

## Proteins significant in both, only 1 or only 2 (with beta >1.5)
sig.both <- subset(ogtt.both, (fdr.aov.1<0.05 & abs(beta.high.1)>1.5 & pval.aov.t.factor.2<0.05) |
                     (fdr.aov.2<0.05 & abs(beta.high.2)>1.5 & pval.aov.t.factor.1<0.05) ) %>%
            filter(Assay!="PYY") %>% pull(Assay) # PYY has a differential response
sig.1.only <- subset(ogtt.both, (fdr.aov.1<0.05 & pval.aov.t.factor.2>.05 & abs(beta.high.1)>1 ) ) %>% pull(Assay) #
sig.2.only <- subset(ogtt.both, (fdr.aov.2<0.05 & pval.aov.t.factor.1>.05 & abs(beta.high.2)>1) ) %>% pull(Assay) # 

svg("graphics/01_ogtt_consistent.svg", width = 6, height = 12)
plot_prots(sig.both, 2) + theme(legend.position = c(0.8, 0.1))
dev.off()

plot_prots(sig.1.only,4)
plot_prots(sig.2.only,4)

svg("graphics/01_ogtt_consistent_indvdl.svg", width = 14, height = 4)
plot_prots_ind(sig.both)
dev.off()


## interaction analysis
diff_prots <- c("WARS","PYY","MLN","CDSN")
svg("graphics/01_ogtt_diff.svg", width=6, height=6)
plot_prots(diff_prots, 2) + theme(legend.position = "none")
dev.off()

svg("graphics/01_ogtt_diff_indvdl.svg", width=7, height=4)
plot_prots_ind(diff_prots)
dev.off()

# genes coexpressed with anxa10 in GI and pancreas
coexp_prots <- c("ANXA10","TFF2","AGR2","VSIG2","PLA2G10","CA9")
svg("graphics/01_ogtt_coexp_prots.svg", width=3, height=15)
plot_prots(coexp_prots, 1)
dev.off()

##########################################################

save.image("RData/01.RData")
