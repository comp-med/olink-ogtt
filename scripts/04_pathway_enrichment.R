###################################################
#### Analysis of Olink fasting OGTT study      ####
#### proteomics data                           ####
#### pathway enrichment                        ####
###################################################

rm(list=ls())
setwd("/sc-projects/sc-proj-computational-medicine/people/Burulca/06_fasting_ogtt/00_clean_code/")
options(stringsAsFactors = F)

## Packages
library(dplyr)
library(ggplot2)
library(gprofiler2)
library(shadowtext)
library(data.table)
library(stringr)
library(ggrepel)
library(gridExtra)
library(scales)


###################################################
#### Read data                                 ####
###################################################

res.ogtt.1.linear <- read.table("data/res.ogtt.1.linear.txt", sep = "\t", header = T)
res.ogtt.2.linear <- read.table("data/res.ogtt.2.linear.txt", sep = "\t", header = T)

## Combine ogtt1 and 2
ogtt1 <-res.ogtt.1.linear[,c("OlinkID","Assay","pval.aov.t.factor","fdr.aov","beta.high","first_sig")]
ogtt2 <-res.ogtt.2.linear[,c("OlinkID","Assay","pval.aov.t.factor","fdr.aov","beta.high","first_sig")]
ogtt.both <- ogtt1 %>% left_join(ogtt2, by=c("OlinkID","Assay"), suffix=c(".1",".2"))

## Create sets of genes to test
ogtt1_only <- ogtt.both %>% filter(fdr.aov.1<0.05 & pval.aov.t.factor.2>0.05) %>% pull(Assay) %>% unique() # n=24
ogtt2_only <- ogtt.both %>% filter(fdr.aov.2<0.05 & pval.aov.t.factor.1>0.05) %>% pull(Assay) %>% unique() # n=35
ogtt1_sig <- ogtt.both %>% filter(fdr.aov.1 < .05) %>% pull(Assay) %>% unique() # n=52
ogtt2_sig <- ogtt.both %>% filter(fdr.aov.2 < .05) %>% pull(Assay) %>% unique() # n=71
ogtt_consistent <- ogtt.both %>% filter(fdr.aov.1<0.05&pval.aov.t.factor.2<0.05 | fdr.aov.2<0.05&pval.aov.t.factor.1<0.05) %>%
  pull(Assay) %>% unique() # n=48
olink_universe <- ogtt.both %>% pull(Assay) %>% unique() # n=2923


###################################################
#### Gprofiler                                 ####
###################################################

## Function to get enrichment result
test_enrichment <- function(geneset, universe, label) {
  res.enrich <- gost(query = geneset, 
                        organism = "hsapiens", ordered_query = FALSE, 
                        multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                        measure_underrepresentation = FALSE, evcodes = TRUE, 
                        user_threshold = 0.05, correction_method = "fdr", 
                        domain_scope = "annotated", custom_bg = universe,
                        numeric_ns = "", as_short_link = FALSE,
                        highlight = TRUE
                        #sources = c("KEGG", "REAC")
                     )
  
  res.enrich <- as.data.table(res.enrich$result)
  res.enrich$OR <- res.enrich$intersection_size * (res.enrich$effective_domain_size + res.enrich$term_size+res.enrich$query_size-res.enrich$intersection_size*2) /
    ( (res.enrich$query_size-res.enrich$intersection_size) * (res.enrich$term_size-res.enrich$intersection_size) )
  res.enrich$geneset <- label
  return(res.enrich)
}

## Apply to the gene sets
ogtt1_res <- test_enrichment(ogtt1_sig, olink_universe, "OGTT 1") ##>> no signif results
ogtt2_res <- test_enrichment(ogtt2_sig, olink_universe, "OGTT 2")

ogtt1_only_res <- test_enrichment(ogtt1_only, olink_universe, "OGTT 1 only")
ogtt2_only_res <- test_enrichment(ogtt2_only, olink_universe, "OGTT 2 only")

ogtt_cons_res <- test_enrichment(ogtt_consistent, olink_universe, "Consistent in both")

## Combine and write to file
res.enrich.all <- rbind(ogtt1_res, ogtt2_res, ogtt1_only_res, ogtt2_only_res, ogtt_cons_res) %>% filter(intersection_size > 2)
res.enrich.all$log10p <- -log10(res.enrich.all$p_value)
res.enrich.all$parents <-NULL
res.enrich.all$source <- factor(res.enrich.all$source, levels=c("GO:MF","GO:BP","GO:CC","KEGG","REAC","WP","HPA","CORUM","MIRNA"))
write.csv(res.enrich.all, "output/04_enriched_pathways.csv", row.names = F)


###################################################
#### Prune enrichment results                  ####
###################################################

source("scripts/functions/prune_pathways.R")

pruned_results <- res.enrich.all %>%
  group_by(geneset, source) %>% 
  do(prune_pathways(.)) %>% 
  ungroup()    


####################################################
## Plot to pdf                                  ####
####################################################

tmp <- pruned_results %>% filter(geneset %in% c("OGTT 1 only","OGTT 2 only","Consistent in both"),
                                 source %in% c("GO:MF","GO:BP","GO:CC","KEGG","REAC","WP"))

# tmp <- pruned_results %>% filter(geneset %in% c("OGTT 1","OGTT 2"),
#                                  source %in% c("GO:MF","GO:BP","GO:CC","KEGG","REAC","WP"))

hl_terms <- tmp
hl_terms$wrapped_term <- str_wrap(hl_terms$term_name, width = 30) # Adjust width as needed

gp <- ggplot(tmp, aes(x=log10p, y=OR, shape=as.factor(geneset), col=source)) +
  geom_point(size=4) +
  geom_text_repel(data = hl_terms, aes(label = wrapped_term, col=source), size=4) +
  geom_vline(xintercept = -log10(.05), linetype = "dashed", colour = "#AA4499") +
  theme_minimal() + labs(x="-log10(p-adj)", shape="Gene set") +
  theme(panel.spacing = unit(1, "lines"),
        panel.grid.minor = element_blank(),
        text = element_text(size=12),
        legend.position = c(0.8,0.8))

pdf("graphics/04_gprofiler_enrichment.pdf", width = 9, height = 9)
gp
dev.off()  

