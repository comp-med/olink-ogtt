
setwd("/sc-projects/sc-proj-computational-medicine/people/Burulca/06_fasting_ogtt/00_clean_code/")

library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(ggplot2)

################################################################################

#### Read data 
res.ogtt.1.linear <- read.table("output/01_res.ogtt.1.linear.txt", sep = "\t", header = T)
res.ogtt.2.linear <- read.table("output/01_res.ogtt.2.linear.txt", sep = "\t", header = T)

## Combine ogtt1 and 2
ogtt1 <-res.ogtt.1.linear[,c("OlinkID","Assay","pval.aov.t.factor","fdr.aov","beta.high","first_sig")]
ogtt2 <-res.ogtt.2.linear[,c("OlinkID","Assay","pval.aov.t.factor","fdr.aov","beta.high","first_sig")]
ogtt.both <- ogtt1 %>% left_join(ogtt2, by=c("OlinkID","Assay"), suffix=c(".1",".2"))

## Create sets of genes to test
ogtt.both <- ogtt.both %>% mutate(
    group = case_when(
      fdr.aov.1 < 0.05 & pval.aov.t.factor.2 > 0.05 ~ "Only in OGTT1",
      fdr.aov.2 < 0.05 & pval.aov.t.factor.1 > 0.05 ~ "Only in OGTT2",
      (fdr.aov.1 < 0.05 & pval.aov.t.factor.2 < 0.05) |
        (fdr.aov.2 < 0.05 & pval.aov.t.factor.1 < 0.05) ~ "Consistent in Both",
      TRUE ~ "Not Significant"  # Catch-all for non-significant cases
))

# Select row with minimum value between fdr.aov.1 and fdr.aov.2
ogtt_groups <- ogtt.both %>%
  group_by(Assay) %>%
  slice_min(pmin(fdr.aov.1, fdr.aov.2), with_ties = FALSE) %>%  
  ungroup()  %>% select(OlinkID, group) ## n=2923

olink_universe <- ogtt.both %>% pull(Assay) %>% unique() # n=2923

################################################################################

## Do heatmap using hpa data directly
hpa_data <- read.table("/sc-projects/sc-proj-computational-medicine/data/07_public/01_Human_Protein_Atlas/data/olink_hpa_data_2023-07-24_processed.tsv", header = T)

hpa_data <- hpa_data %>% dplyr::select("OlinkID","Assay","max_ntpm_tissue","RNA.tissue.specificity",starts_with("RNA.tissue.specific.nTPM."))
colnames(hpa_data) <- gsub("RNA.tissue.specific.nTPM.", "", colnames(hpa_data))
colnames(hpa_data) <- gsub(".1", "", colnames(hpa_data))

hpa_data <- hpa_data %>% mutate(brain = coalesce(brain, choroid.plexus)) %>% select(-choroid.plexus)  # combine brain and choroid plexus

hpa_data <- hpa_data %>% inner_join(ogtt_groups, by="OlinkID")
hpa_data %>% count(group)

hpa_data$OlinkID <- NULL
hpa_data[is.na(hpa_data)] <- 0
hpa_data <- unique(hpa_data)

rownames(hpa_data) <- hpa_data$Assay   # Set row names to Assay
hpa_data$Assay <- NULL

# Remove nonsig
hpa_sig <- hpa_data %>% filter(group != "Not Significant")
# Remove the Group column from the matrix for the heatmap
heatmap_data <- as.matrix(hpa_sig[, 3:37]) %>% t()
# Keep only tissues with expression among the ogtt proteins
heatmap_data <- heatmap_data[rowSums(heatmap_data) != 0,]
# log transform
heatmap_data <- log(heatmap_data + 1)

# Define custom colors for the groups
group_colors <- c("Only in OGTT1" = "#88CCEE", "Consistent in Both" = "#AA4499", "Only in OGTT2" = "#CC6677")

# Modify the annotation for groups using the defined colors
col_annotation <- columnAnnotation(
  Group = factor(hpa_sig$group),  # Assuming `group_labels` contains the group assignments
  col = list(Group = group_colors)  # Set custom colors for the groups
)

hpa_sig$group <- factor(hpa_sig$group, levels = names(group_colors))

# Create the heatmap
pdf("graphics/03_tissue_expr_ogtt_1_2_summary.pdf", height = 8, width = 20)
  Heatmap(heatmap_data, 
     name = "log(nTPM+1)", 
     clustering_distance_rows = "euclidean", 
     clustering_distance_columns = "euclidean", 
     clustering_method_rows = "complete",
     clustering_method_columns = "complete",
     col = colorRamp2(c(min(heatmap_data), max(heatmap_data)), c("white", "black")),
     rect_gp = gpar(col = "grey50", lwd = 1),
     top_annotation = col_annotation,
     column_split = hpa_sig$group,
     column_gap = unit(6, "mm"))
dev.off()

################################################################################

#### Create heatmap with only stomach, intestine and pituitary

hpa_sig <- hpa_data %>% filter(group != "Not Significant") %>% select(stomach, intestine, pituitary.gland, group)
hpa_sig <- hpa_sig[rowSums(hpa_sig[,1:3]) != 0,]
heatmap_data <- as.matrix(hpa_sig[,1:3]) %>% t() 
heatmap_data <- log(heatmap_data + 1)

# Modify the annotation for groups using the defined colors
col_annotation <- columnAnnotation(
  Group = factor(hpa_sig$group),  # Assuming `group_labels` contains the group assignments
  col = list(Group = group_colors)  # Set custom colors for the groups
)

hpa_sig$group <- factor(hpa_sig$group, levels = names(group_colors))

# Create the heatmap
pdf("graphics/03_tissue_expr_ogtt_1_2_summary_stomach_intestine.pdf", height = 4, width = 12)
Heatmap(heatmap_data, 
        name = "log(nTPM+1)", 
        clustering_distance_rows = "euclidean", 
        clustering_distance_columns = "euclidean", 
        clustering_method_rows = "complete",
        clustering_method_columns = "complete",
        col = colorRamp2(c(min(heatmap_data), max(heatmap_data)), c("white", "black")),
        rect_gp = gpar(col = "grey50", lwd = 1),
        top_annotation = col_annotation,
        column_split = hpa_sig$group,
        column_gap = unit(6, "mm"),
        show_row_dend = FALSE,  # Hide row dendrogram
        show_column_dend = FALSE)  # Hide column dendrogram
dev.off()

################################################################################

## Fisher exact test
hpa_data$max_ntpm_tissue <- gsub(" 1", "", hpa_data$max_ntpm_tissue)
hpa_data$max_ntpm_tissue <- gsub("choroid plexus", "brain", hpa_data$max_ntpm_tissue)
hpa_data$specific_tissue <- ifelse(hpa_data$RNA.tissue.specificity %in% c("Tissue enriched","Tissue enhanced"), hpa_data$max_ntpm_tissue, NA)

#### Get numbers for contingency table
ogtt.by.tissue.1 <- hpa_data %>% count(specific_tissue, sig=group=="Only in OGTT1") %>%
  pivot_wider(names_from = sig, values_from = n) %>% filter(!is.na(specific_tissue))
ogtt.by.tissue.1[is.na(ogtt.by.tissue.1)] <- 0

ogtt.by.tissue.2 <- hpa_data %>% count(specific_tissue, sig=group=="Only in OGTT2") %>%
  pivot_wider(names_from = sig, values_from = n) %>% filter(!is.na(specific_tissue))
ogtt.by.tissue.2[is.na(ogtt.by.tissue.2)] <- 0

ogtt.by.tissue.consistent <- hpa_data %>% count(specific_tissue, sig=group=="Consistent in Both") %>%
  pivot_wider(names_from = sig, values_from = n) %>% filter(!is.na(specific_tissue))
ogtt.by.tissue.consistent[is.na(ogtt.by.tissue.consistent)] <- 0

fisher_test_tissue <- function(ogtt_by_tissue, sig, nonsig) {
  
  ## Loop through tissues
  fisher_df = data.frame(matrix(vector(), 0, 3 ))
  for (t in ogtt_by_tissue$specific_tissue) {
    tissue <- ogtt_by_tissue %>% filter(specific_tissue==t)
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

## How many sig and nonsig
hpa_data %>% count(group)

## apply to dfs and plot
fisher_df_1 <- fisher_test_tissue(ogtt.by.tissue.1, 24, 2923-24)
fisher_df_2 <- fisher_test_tissue(ogtt.by.tissue.2, 35, 2923-35)
fisher_df_cons <- fisher_test_tissue(ogtt.by.tissue.consistent, 47, 2923-47)

fisher_df_1$ogtt.session = "1 only"
fisher_df_2$ogtt.session = "2 only"
fisher_df_cons$ogtt.session = "Both"
fisher_df <- rbind(fisher_df_1, fisher_df_2, fisher_df_cons)

fisher_df %>% ggplot(aes(x=-log10(fdr),y=reorder(tissue,-fdr),size=OR, col=ogtt.session)) + 
  geom_point(alpha=0.7) + geom_vline(xintercept = -log10(.05), linetype = "dashed", colour = "#CC6677") +
  labs(x="-log10(p-adj)") +
  #scale_color_manual(name="OGTT session",values=c("#6699CC", "#CC6677")) +
  theme_minimal() +
  theme(axis.title.y = element_blank())


