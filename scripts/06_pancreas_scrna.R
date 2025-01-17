setwd("/sc-projects/sc-proj-computational-medicine/people/Burulca/06_fasting_ogtt/00_clean_code/")

library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(reshape2)
library(ggrepel)

## Read seurat object - related article https://doi.org/10.2337/db23-0130
## Downloaded from https://www.gaultonlab.org/pages/Islet_expression_HPAP.html 
hpap <- readRDS("/sc-projects/sc-proj-computational-medicine/people/Burulca/06_fasting_ogtt/03_single_cell/data/hpap_islets/hpap_islet_scRNAseq.rds")
# Remove spaces from the column names the meta.data
colnames(hpap@meta.data) <- gsub(" ", "_", colnames(hpap@meta.data))
# Aggregate by cell type
hpap <- SetIdent(hpap, value = hpap@meta.data$Cell_Type)
## Subset to only healthy tissue
hpap_healthy <- subset(hpap, subset = Diabetes_Status=="ND")

################################################################################

## Plot umap of cell types
svg("graphics/06_pancreas_umap.svg",width = 5,height = 5)
DimPlot(hpap, reduction = "umap", label = TRUE, label.size = 3, raster = T) +
  xlab("UMAP_1") + ylab("UMAP_2") +  NoLegend()
dev.off()

## Plot anxa10 expression on umap - highest exp cells on top

# Get expression values of anxa10
expression_data <- FetchData(hpap, vars = "ANXA10")

# Add expression values to the meta.data slot for ordering
hpap@meta.data$expression_value <- expression_data[["ANXA10"]]

# Order cells based on expression values (decreasing to plot high expression on top)
cell_order <- order(hpap@meta.data$expression_value, decreasing = F)

# Create the FeaturePlot for the current gene with ordered cells
#pdf("graphics/06_pancreas_anxa_exp_plot.pdf", height=5, width=5)
png("graphics/06_pancreas_anxa_exp_plot_compress.png", height=5, width=5, units="in", res=300)
FeaturePlot(hpap, features = "ANXA10", raster = F, cells = cell_order) +
  xlab("UMAP_1") + ylab("UMAP_2") + ggtitle("ANXA10") +
  theme(legend.position = "none")
dev.off()


## Plot dotplot with coexpressed genes
## List of genes
genes <- c("ANXA10","MUC1", "AGR2", "VSIG2", "TFF2","CA9","PLA2G10") 
#genes <- c("ANXA10","MUC1","AGR2","VSIG2","TFF2","TFF1","ATP4A","ATP4B","ATP12A","ATP1A1","ATP1A2","FXYD5")
svg("graphics/06_pancreas_dotplot.svg", height=2.5, width=4)
DotPlot(hpap_healthy, features = genes, idents=c("MUC5B+ Ductal","Acinar","Ductal")) +
  scale_color_gradient(low = "lightblue", high = "darkblue", limits = c(-1.5, 1.5), breaks=c(-1,0,1)) +
  scale_size_continuous(breaks=c(20,40,60,80), limits = c(0,100)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  guides(size = "none", col="none")
dev.off()

svg("graphics/06_pancreas_dotplot_all_cells.svg", height=5, width=6)
DotPlot(hpap, features = genes) +
  scale_color_gradient(low = "lightblue", high = "darkblue", limits = c(-1.5, 3)) +
  scale_size_continuous(breaks=c(20,40,60,80), limits = c(0,100)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
dev.off()


################################################################################

#### Coexpression analysis

rm(hpap) # clean up some memory
gc() 

#### Co-expression with anxa10 in healthy cells

# Define the function
compute_gene_correlation <- function(seurat_obj, celltype, target_gene = "ANXA10") {
  cell_subset <- subset(seurat_obj, subset = Cell_Type %in% celltype)
  #selected_cells <- WhichCells(cell_subset, expression = nFeature_RNA > 200) # filter worked before but now is not working
  selected_genes <- rownames(cell_subset)[Matrix::rowSums(cell_subset@assays$RNA@counts) > 0]
  cell_subset <- subset(cell_subset, features = selected_genes)
  count_matrix <- as.data.frame(t(as.matrix(cell_subset@assays$RNA$counts)))
  target_gene_counts <- count_matrix[, target_gene, drop = FALSE]
  
  corr_results <- data.frame()
  # Loop through each gene in the count matrix and compute the correlation with the target gene
  for (gene in colnames(count_matrix)) {
    cor_test <- cor.test(count_matrix[[gene]], target_gene_counts[[1]], method = "pearson")
    corr_results <- rbind(corr_results, data.frame(
      gene = gene, 
      cor = cor_test$estimate, 
      p_value = cor_test$p.value
    ))
  }
  return(corr_results)
}

## Compute correlation of anxa10 with other genes for relevant celltypes, and all healthy cells
acinar.corr <- compute_gene_correlation(hpap_healthy, celltype = "Acinar")
muc5.corr <- compute_gene_correlation(hpap_healthy, celltype = "MUC5B+ Ductal")

## Computing correlation for all healthy cells takes too much mem and fails
## Try doing it in batches
selected_genes <- rownames(hpap_healthy)[Matrix::rowSums(hpap_healthy@assays$RNA@counts) > 200]
cell_subset <- subset(hpap_healthy, features = selected_genes)
count_matrix <- as.data.frame(t(as.matrix(cell_subset@assays$RNA@counts)))
anxa_counts <- count_matrix[, "ANXA10", drop = FALSE]

# Initialize result dataframe
corr_results <- data.frame()
# Process genes in batches
batch_size <- 1000
i=1
gene_batches <- split(colnames(count_matrix), ceiling(seq_along(colnames(count_matrix)) / batch_size))
for (batch in gene_batches) {
  for (gene in batch) {
    cor_test <- cor.test(count_matrix[[gene]], anxa_counts[[1]], method = "pearson")
    corr_results <- rbind(corr_results, data.frame(
      gene = gene, 
      cor = cor_test$estimate, 
      p_value = cor_test$p.value
    ))
  }
  print(i)
  i=i+1
}

## join together
anxa.corr <- corr_results %>% 
  full_join(acinar.corr, by = "gene", suffix = c(".healthy", ".acinar")) %>% 
  full_join(muc5.corr, by = "gene", suffix = c("", ".muc5"))
# Rename columns for MUC5 manually 
colnames(anxa.corr)[6:7] <- c("cor.muc5","p_value.muc5")

#### Check distribution of correlations
anxa.corr.long <- anxa.corr[,c("gene","cor.healthy","cor.acinar","cor.muc5")]
anxa.corr.long <- melt(anxa.corr.long, id.vars="gene", 
                       measure.vars=c("cor.healthy","cor.acinar","cor.muc5"))
anxa.corr.long$variable <- gsub("cor.","",anxa.corr.long$variable)

corr.quantiles <- anxa.corr.long %>% group_by(variable) %>% summarise(quantile95=quantile(value, probs = 0.95, na.rm=T),
                                                                      quantile99=quantile(value, probs = 0.99, na.rm=T))

# determine 95% percentile of the correlations for each cluster (or all healthy cells)
quantiles <- corr.quantiles$quantile95
names(quantiles) <- corr.quantiles$variable

################################################################################

#### Compare with ogtt proteomics results

ogtt1 <- read.table("output/01_res.ogtt.1.linear.txt", sep="\t", header=T)
ogtt2 <- read.table("output/01_res.ogtt.2.linear.txt", sep="\t", header=T)

# keep min fdr for proteins with multiple olink ids
ogtt1 <- ogtt1 %>% group_by(Assay) %>% slice(which.min(fdr.aov)) %>% summarize(fdr.aov = first(fdr.aov))
ogtt2 <- ogtt2 %>% group_by(Assay) %>% slice(which.min(fdr.aov)) %>% summarize(fdr.aov = first(fdr.aov))

ogtt.both <- ogtt1 %>% full_join(ogtt2, by="Assay")
colnames(ogtt.both) <- c("Assay", "ogtt1.fdr.aov", "ogtt2.fdr.aov")
ogtt.both$sig <- ifelse(ogtt.both$ogtt1.fdr.aov < 0.05 | ogtt.both$ogtt2.fdr.aov < 0.05, 1, 0)

anxa.corr.ogtt <- anxa.corr %>% full_join(ogtt.both, by=c("gene"="Assay"))
anxa.corr.ogtt$coexp_healthy <- ifelse(anxa.corr.ogtt$cor.healthy > quantiles[grepl("healthy", names(quantiles))], 1, 0)
anxa.corr.ogtt$coexp_acinar <- ifelse(anxa.corr.ogtt$cor.acinar > quantiles[grepl("acinar", names(quantiles))], 1, 0)
anxa.corr.ogtt$coexp_muc5 <- ifelse(anxa.corr.ogtt$cor.muc5 > quantiles[grepl("muc5", names(quantiles))], 1, 0)

write.csv(anxa.corr.ogtt, "output/06_pancreas_anxa_coexpression_ogtt_together.csv", row.names = F)


################################################################################

#### Plot dotplot with 25 coexpressed genes in stomach (for supp fig)

anxa.corr.ogtt.pancreas <- read.csv("output/06_pancreas_anxa_coexpression_ogtt_together.csv", header = T)
anxa.corr.ogtt.stomach <- read.csv("output/05_stomach_anxa_coexpression_ogtt_together.csv", header = T)

# Calculate the percentage of cells expressing ANXA10
anxa10_expression <- FetchData(hpap_healthy, vars = "ANXA10")
anxa10_expression$cell_type <- Idents(hpap_healthy)

# Calculate percentage of cells expressing ANXA10 for each cell type
percentage_expressing <- anxa10_expression %>% group_by(cell_type) %>% summarize(percent_expressing = sum(ANXA10 > 0) / n() * 100)

# Sort cell types by percentage
sorted_cell_types <- percentage_expressing %>% arrange(desc(percent_expressing)) %>% pull(cell_type)

# Order cell types in the Seurat object
hpap_healthy <- SetIdent(hpap_healthy, value = factor(Idents(hpap_healthy), levels = rev(sorted_cell_types)))

# Pick top 25 genes in pancreas
top25_corr <- head(anxa.corr.ogtt.pancreas[order(-anxa.corr.ogtt.pancreas$cor.healthy),]$gene, 25)
# top 25 in stomach
top25_corr_names <- head(anxa.corr.ogtt.stomach[order(-anxa.corr.ogtt.stomach$healthy.cor),]$hgnc_symbol.x, 25)

svg("graphics/06_pancreas_top25_dotplot.svg", height=5, width=10)
DotPlot(hpap_healthy, features = c(top25_corr_names,"CA9","PLA2G10")) + ## these are the top25 in stomach
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  scale_size_continuous(breaks=c(20,40,60,80), limits = c(0,100)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
dev.off()

################################################################################

#### Combine with the stomach data
stomach <- anxa.corr.ogtt.stomach %>% filter(coexp_healthy==1 & length(hgnc_symbol.x)>0) %>% dplyr::select("hgnc_symbol.x","gene","healthy.cor")

pancreas <- anxa.corr.ogtt.pancreas %>% filter(coexp_healthy==1 & length(gene)>0) %>% dplyr::select("gene","cor.healthy","sig")

stomach_and_pancreas <- stomach %>% inner_join(pancreas, by=c("hgnc_symbol.x"="gene"))
plot_data <- stomach_and_pancreas %>% filter(hgnc_symbol.x!="ANXA10")

## scatterplot of correlations in stomach vs pancreas
svg("graphics/06_stomach_and_pancreas_coexp.svg", width=8, height=8)
ggplot(plot_data, aes(x= healthy.cor, y=cor.healthy, col=as.factor(sig))) +
  geom_point() +
  # Label points with sig == 1
  geom_text_repel(data = subset(plot_data, sig == 1), aes(label = hgnc_symbol.x, vjust = -0.5) ) +
  # Optionally label additional points if they don't overlap with the first
  geom_text_repel(data = subset(plot_data, is.na(sig) | sig != 1), aes(label = hgnc_symbol.x), max.overlaps = 2) +
  scale_color_manual(name="Changed in OGTT", values = c("#332288","#AA4499","#888888")) +
  theme_minimal() + xlab("Corr. with ANXA10 in stomach") + ylab("Corr. with ANXA10 in pancreas") +
  theme(legend.position=c(0.9,0.8))
dev.off()

## dotplot of ogtt genes coexp in both 
coexp_in_both <- subset(stomach_and_pancreas, sig == 1)$hgnc_symbol.x
DotPlot(hpap_healthy, features = coexp_in_both) + 
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  scale_size_continuous(breaks=c(20,40,60,80), limits = c(0,100)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 

coexp_in_both_ensembl <- subset(stomach_and_pancreas, sig == 1)$gene
DotPlot(nowicki_healthy, features = coexp_in_both_ensembl) + 
  scale_x_discrete(labels = coexp_in_both) +
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  scale_size_continuous(breaks=c(20,40,60,80), limits = c(0,100)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
