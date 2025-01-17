## Adapted from /sc-projects/sc-proj-computational-medicine//people/Maik/14_phecode_GWAS/10_disease_examples/01_PNLIPRP3_rosacea/

## standard starting parameters
rm(list=ls())
setwd("/sc-projects/sc-proj-computational-medicine/people/Burulca/06_fasting_ogtt/14_eqtl_coloc")
options(stringsAsFactors = F)

## --> packages needed <-- ##
library(data.table)
library(doMC)
library(dplyr)
library(tidyr)
library(zoo)
library(ggplot2)

################################################################################

## prev. steps done in scrip 02_cis_eQTL_coloc
## --> import parameters <-- ##

# # ## get regional coordinates
# pheno <- args[1]
# chr.s <- as.numeric(args[2])
# pos.l <- as.numeric(args[3])
# pos.s <- as.numeric(args[4])
# pos.e <- as.numeric(args[5])

## get regional coordinates
tmp   <- read.table("input/input.eQTL.pipeline.txt")
j     <- 1
pheno <- tmp$V1[j]
chr.s <- tmp$V2[j]
pos.l <- tmp$V3[j]
pos.s <- tmp$V4[j]
pos.e <- tmp$V5[j]

################################################################################
######## Get regenie output, LD and genes list ########

## anxa10 regenie (file with filtered region)
res           <- readRDS(paste0("input/anxa10_chr",chr.s,"_res.rds"))
res$snp.id <- 1:nrow(res) 

## read in matrix
ld         <- fread(paste("tmpdir/ld", pheno, chr.s, pos.s, pos.e, "ld", sep="."), data.table = F)

## import list of genes
tmp.genes      <- fread("/sc-projects/sc-proj-computational-medicine/people/Maik/23_GWAS_retinal_risk_states/03_eQTL_overlap/input/Genes.GRCh37.complete.txt")

## restrict to protein encoding genes for now
tmp.genes      <- subset(tmp.genes, gene_biotype %in% c("protein_coding", "processed_transcript") & chromosome_name == ifelse(chr.s == 23, "X", chr.s) & start_position >= pos.s-5e5 & end_position <= pos.e+5e5)
#tmp.genes      <- subset(tmp.genes, gene_biotype %in% c("protein_coding", "processed_transcript") )

## subset to genes no more than 500kb away from the lead variant (gene body not TSS)
tmp.genes$dist <- apply(tmp.genes[, c("start_position", "end_position")], 1, function(x) min(abs(x-pos.l)))
tmp.genes      <- subset(tmp.genes, dist <= 5e5) # 37 genes chr1, 4 genes chr 4

################################################################################
######## Get gtex data from coloc function ########

## import function to get gene stats
source("/sc-projects/sc-proj-computational-medicine/people/Burulca/06_fasting_ogtt/14_eqtl_coloc/scripts/coloc_gene_gtex_tissues.R")

## run across all genes
registerDoMC(10)
res.genes <- mclapply(tmp.genes$ensembl_gene_id, function(g){

  print(g)

  ## run coloc for specific gene
  res.tmp <- coloc.eQTL(res, chr.s, pos.s, pos.e, g, ld, r.t=T)
  res.gene <- res.tmp[[2]]
  return(res.gene)
})

saveRDS(res.genes, "data/res_genes.rds")

## until here ran on slurm

################################################################################
######## Combine the list elements into one df ########

res_genes <- readRDS("/sc-projects/sc-proj-computational-medicine/people/Burulca/06_fasting_ogtt/14_eqtl_coloc/data/res_genes.rds")

names(res_genes) <- tmp.genes$ensembl_gene_id

for (i in 1:37) {
  print(i)
  print(names(res_genes)[i])
  print(dim(res_genes[[i]]))
}

# ENSG00000251246 and ENSG00000273088 doesnt have results
res_genes <- res_genes[c(-14,-20)]

# Identify all unique column names across all data frames
all_cols <- unique(unlist(lapply(res_genes, colnames)))

# Initialize combined_df with NA
combined_df <- data.frame(matrix(NA, nrow = 0, ncol = length(all_cols)))
colnames(combined_df) <- all_cols

# Fill combined_df with data from each df in res_genes
for (df_name in names(res_genes)) {
  df <- res_genes[[df_name]]
  
  # Create df_all with all columns and fill with NA
  df_all <- data.frame(matrix(NA, nrow = nrow(df), ncol = length(all_cols)))
  colnames(df_all) <- all_cols
  
  # Fill df_all with data from df, matching by column names
  for (col in intersect(colnames(df_all), colnames(df))) {
    df_all[, col] <- df[, col]
  }
  
  # Combine df_all with combined_df
  combined_df <- rbind(combined_df, df_all)
}

################################################################################
######## Get relevant tissues ########

tissues <- c("Esophagus_Gastroesophageal_Junction", "Esophagus_Muscularis", "Liver", "Pancreas", "Stomach")
tissue_pattern <- paste(tissues, collapse = "|")

## Subset to tissues of interest
tissue_df <- combined_df[, c("phenotype_id","MarkerName", "snp.id", "CHROM", "GENPOS", "ID","ALLELE0","ALLELE1", "BETA", "SE", "LOG10P", "Allele1", "Allele2",
                               "variant_id", "pos.hg38", grep(tissue_pattern, names(combined_df), value = TRUE, ignore.case = TRUE))]
# add anxa10 pval
tissue_df$pval_nominal.anxa10 <- 10^(-tissue_df$LOG10P)
tissue_df$phenotype_id <- sub("\\..*", "", tissue_df$phenotype_id)

# Base columns
base_columns <- c("Ensembl", "MarkerName", "snp.id", "chr", "pos", "rsid", "ALLELE0","ALLELE1", "Effect.anxa10", "StdErr.anxa10", "log10p.anxa10", "Allele1", "Allele2", "variant_id", "pos.hg38")

data_types <- c("Effect", "StdErr", "pval_nominal", "N")

# Create new column names by first iterating over data types and then over tissues
new_columns <- base_columns

for (data_type in data_types) {
  for (tissue in tissues) {
    new_columns <- c(new_columns, paste0(data_type, ".", tissue))
  }
}

# renaming
names(tissue_df) <- c(new_columns,"pval_nominal.anxa10")
write.csv(tissue_df, "graphics/stacked_association_plot_data_chr1.csv",row.names = F)


################################################################################
######## Convert to long format, add LD ########

tissue_df <- read.csv("graphics/stacked_association_plot_data_chr1.csv",header=T)

tissue_df_long <- tissue_df %>%
  pivot_longer(cols = starts_with("Effect.") | starts_with("StdErr.") | starts_with("log10p.") | starts_with("pval_nominal.") | starts_with("N."),
               names_to = c(".value", "tissue"),
               names_sep = "\\.")

## Due to mistake in coloc_gene_gtex_tissues, effect of anxa was flipped (the one stored in Effect.aligned)
## Recalculate based on the unchanged BETA column
tissue_df_long <- tissue_df_long %>%
  mutate(Effect = ifelse(tissue == "anxa10" & Allele1 == ALLELE0, -Effect, Effect))

tissue_df_long <- tissue_df_long %>%
  left_join(tmp.genes[,c("ensembl_gene_id","external_gene_name")],by=c("Ensembl"="ensembl_gene_id")) %>%
  #filter(tissue!="anxa10") %>%
  filter(external_gene_name %in% c("TRIM46","MUC1","THBS3","MTX1","GBA","FAM189B"))

# Forward fill NAs in the log10p column
tissue_df_long$log10p <- na.locf(tissue_df_long$log10p, na.rm = FALSE)
tissue_df_long

## ld df
ld_rs2990223 <- ld[,1257] %>% as.data.frame()
ld_rs2990223$snp.id <- 1:nrow(ld_rs2990223)
colnames(ld_rs2990223)[1] <- "ld"
tissue_df_long <- tissue_df_long %>% left_join(ld_rs2990223,by="snp.id")

setnames(tissue_df_long, "chr", "CHROM")
setnames(tissue_df_long, "pos", "GENPOS")
colnames(tissue_df_long)

write.csv(tissue_df, "graphics/stacked_association_plot_data_chr1.csv",row.names = F)
write.csv(tissue_df_long, "graphics/stacked_association_plot_data_chr1_long.csv",row.names = F)


################################################################################

## ggplot
tissue_df_long$tissue <- gsub("Esophagus_Gastroesophageal_Junction", "GEJ", tissue_df_long$tissue)

svg("graphics/coloc/anxa10_trans_pqtl_all5.svg", width = 10, height = 6)
ggplot(tissue_df_long[tissue_df_long$tissue!="anxa10",], aes(x=GENPOS,y=-log10(pval_nominal),col=ld)) +
  geom_point(size=0.5) +
  scale_color_viridis_c(direction=-1) +
  facet_grid(rows=vars(external_gene_name),cols=vars(tissue), scales = "free_y") +
  geom_hline(yintercept = 7.3, linetype='dashed', col="#888888") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90),
        axis.title = element_blank())
dev.off()

# ################################################################################
# 
# ## create label plot
# lab.stacked         <- data.frame(name=c("anxa10", "Esophagus_Muscularis", "Esophagus_Gastroesophageal_Junction", "Liver", "Pancreas", "Stomach"),
#                                   label=c("ANXA10 Protein - Plasma", "ANXA10 - Esophagus_Muscularis", "ANXA10 - Esophagus_Gastroesophageal_Junction", "ANXA10 - Liver", "ANXA10 - Pancreas", "ANXA10 - Stomach")) 
# 
# ## import function to do so
# source("scripts/plot_stacked_locuszoom_plot.R")
# pdf("graphics/ANXA10.trans.eQTLs.pdf", width=3.15, height=6)
# # png("../graphics/PNLIPRP3.Rosacea.cis.eQTLs.png", width=8, height=8, res=900, units = "cm")
# par(mar=c(.1,1.5,.75,.5), cex.axis=.5, cex.lab=.5, bty="l", tck=-.01, mgp=c(.6,0,0), xaxs="i", lwd=.5)
# 
# layout(matrix(1:7), heights = c(rep(.3,7),.2))
# dev.off()
