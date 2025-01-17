# publication https://doi.org/10.1158/2159-8290.cd-22-0824

# seurat object obtained from https://cellxgene.cziscience.com/collections/a18474f4-ff1e-4864-af69-270b956cee5b
# (UMAP of all data, 145,583 cells)

setwd("/sc-projects/sc-proj-computational-medicine/people/Burulca/06_fasting_ogtt/00_clean_code/")

library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(reshape2)
library(biomaRt)

################################################################################

## Read seurat object
nowicki <- readRDS("/sc-projects/sc-proj-computational-medicine/people/Burulca/06_fasting_ogtt/03_single_cell/data/nowicki/local.rds")
nowicki <- SetIdent(nowicki, value = nowicki@meta.data$Celltypes_global)

## Subset to only healthy tissue (80,381 cells)
nowicki_healthy <- subset(nowicki, subset = Tissue_in_paper %in% c("NE", "SMG", "NSCJ", "NGB", "NGC", "ND", "Ileum", "Colon", "Rectum"))


VlnPlot(object = nowicki, features = 'ENSG00000109511',idents=c("Foveolar_Intermediate","Neck-Cells","Foveolar_Differentiated","Parietal"), pt.size = 0) + ggtitle("ANXA10")
VlnPlot(object = nowicki, features = 'ENSG00000185499',idents=c("Foveolar_Intermediate","Neck-Cells","Foveolar_Differentiated","Parietal"), pt.size = 0) + ggtitle("MUC1")
VlnPlot(object = nowicki, features = 'ENSG00000215182',idents=c("Foveolar_Intermediate","Neck-Cells","Foveolar_Differentiated","Parietal"), pt.size = 0) + ggtitle("MUC5AC")


################################################################################

## Plot UMAP by cell type
svg("graphics/05_stomach_umap.svg",width = 5,height = 5)
DimPlot(nowicki, reduction = ".umap_MinDist_0.2_N_Neighbors_15", label = TRUE, label.size = 3, raster = T) +
  xlab("UMAP_1") + ylab("UMAP_2") +  NoLegend()
dev.off()

## Plot anxa10 expression on umap
#pdf("graphics/05_stomach_anxa_exp_plot.pdf", height=5, width=5, compress=TRUE)
png("graphics/05_stomach_anxa_exp_plot_compress.png", height=5, width=5, units="in", res=300)
FeaturePlot(nowicki, features = "ENSG00000109511", raster = F) +
  xlab("UMAP_1") + ylab("UMAP_2") +
  ggtitle("ANXA10") + theme(legend.position = "none")
dev.off()

## Plot dotplot with coexpressed genes (from analysis below)
genes <- c("ENSG00000109511","ENSG00000185499","ENSG00000106541","ENSG00000019102","ENSG00000160181","ENSG00000107159","ENSG00000069764")
           #"ENSG00000160182","ENSG00000105675","ENSG00000186009","ENSG00000075673",
           #"ENSG00000163399","ENSG00000174437","ENSG00000089327")
names <- c("ANXA10","MUC1", "AGR2", "VSIG2", "TFF2","CA9","PLA2G10") 
          #"TFF1","ATP4A","ATP4B","ATP12A","ATP1A1","ATP1A2","FXYD5")

svg("graphics/05_stomach_dotplot.svg", height=3, width=6)
DotPlot(nowicki_healthy, features = genes, idents=c("Foveolar_Intermediate","Neck-Cells","Foveolar_Differentiated","Parietal")) +
  scale_x_discrete(labels = names) +
  scale_color_gradient(low = "lightblue", high = "darkblue", limits = c(-1.5, 1.5), breaks=c(-1,0,1)) +
  scale_size_continuous(breaks=c(20,40,60,80), limits = c(0,100)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10))
dev.off()

################################################################################

#### Coexpression analysis

## Define the function
compute_gene_correlation <- function(seurat_obj, celltype, target_gene = "ENSG00000109511") {
  cell_subset <- subset(seurat_obj, subset = Celltypes_global %in% celltype)
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
foveolar_int.corr <- compute_gene_correlation(nowicki_healthy, celltype = "Foveolar_Intermediate")
foveolar_diff.corr <- compute_gene_correlation(nowicki_healthy, celltype = "Foveolar_Differentiated")
neck.corr <- compute_gene_correlation(nowicki_healthy, celltype = "Neck-Cells")
parietal.corr <- compute_gene_correlation(nowicki_healthy, celltype = "Parietal")
healthy.corr <- compute_gene_correlation(nowicki_healthy, celltype = unique(nowicki_healthy@meta.data$Celltypes_global))

## join together
anxa.corr <- healthy.corr %>% full_join(foveolar_diff.corr, by="gene", suffix = c("",".foveolar_diff")) %>%
  full_join(foveolar_int.corr, by="gene", suffix = c("",".foveolar_int")) %>%
  full_join(neck.corr,by="gene", suffix = c("",".neck")) %>%
  full_join(parietal.corr,by="gene", suffix = c("",".parietal"))

## add gene names
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                   values=anxa.corr$gene, mart = ensembl)
anxa.corr <- anxa.corr %>% left_join(gene_list, by=c("gene"="ensembl_gene_id"))

#### Check distribution of correlations
anxa.corr.long <- anxa.corr[,c("gene","foveolar_int.cor","foveolar_diff.cor","neck.cor","parietal.cor","healthy.cor")]
anxa.corr.long <- melt(anxa.corr.long, id.vars="gene", 
                       measure.vars=c("foveolar_int.cor","foveolar_diff.cor","neck.cor","parietal.cor","healthy.cor"))
anxa.corr.long$variable <- gsub(".cor","",anxa.corr.long$variable)

corr.quantiles <- anxa.corr.long %>% group_by(variable) %>% summarise(quantile95=quantile(value, probs = 0.95, na.rm=T),
                                                                      quantile99=quantile(value, probs = 0.99, na.rm=T))

# determine 95% percentile of the correlations for each cluster (or all healthy cells)
quantiles <- corr.quantiles$quantile95
names(quantiles) <- corr.quantiles$variable

#### Compare with ogtt proteomics results

ogtt1 <- read.table("output/01_res.ogtt.1.linear.txt", sep="\t", header=T)
ogtt2 <- read.table("output/01_res.ogtt.2.linear.txt", sep="\t", header=T)

# keep min fdr for proteins with multiple olink ids
ogtt1 <- ogtt1 %>% group_by(UniProt) %>% slice(which.min(fdr.aov)) %>% summarize(fdr.aov = first(fdr.aov))
ogtt2 <- ogtt2 %>% group_by(UniProt) %>% slice(which.min(fdr.aov)) %>% summarize(fdr.aov = first(fdr.aov))

ogtt.both <- ogtt1 %>% full_join(ogtt2, by="UniProt")
colnames(ogtt.both) <- c("UniProt", "ogtt1.fdr.aov", "ogtt2.fdr.aov")
ogtt.both$sig <- ifelse(ogtt.both$ogtt1.fdr.aov < 0.05 | ogtt.both$ogtt2.fdr.aov < 0.05, 1, 0)

gene_list <- getBM(filters= "uniprotswissprot", attributes= c("uniprotswissprot","ensembl_gene_id","hgnc_symbol"),
                   values=ogtt.both$UniProt, mart = ensembl)
ogtt.both <- ogtt.both %>% left_join(gene_list, by=c("UniProt"="uniprotswissprot"))

anxa.corr.ogtt <- anxa.corr %>% full_join(ogtt.both, by=c("gene"="ensembl_gene_id"))
anxa.corr.ogtt$coexp_foveolar_int <- ifelse(anxa.corr.ogtt$foveolar_int.cor > quantiles[grepl("foveolar_int", names(quantiles))], 1, 0)
anxa.corr.ogtt$coexp_foveolar_diff <- ifelse(anxa.corr.ogtt$foveolar_diff.cor > quantiles[grepl("foveolar_diff", names(quantiles))], 1, 0)
anxa.corr.ogtt$coexp_neck <- ifelse(anxa.corr.ogtt$neck.cor > quantiles[grepl("neck", names(quantiles))], 1, 0)
anxa.corr.ogtt$coexp_parietal <- ifelse(anxa.corr.ogtt$parietal.cor > quantiles[grepl("parietal", names(quantiles))], 1, 0)
anxa.corr.ogtt$coexp_healthy <- ifelse(anxa.corr.ogtt$healthy.cor > quantiles[grepl("healthy", names(quantiles))], 1, 0)

write.csv(anxa.corr.ogtt, "output/05_stomach_anxa_coexpression_ogtt_together.csv", row.names = F)

################################################################################

#### Plot dotplot with 25 coexpressed genes (supp fig)

anxa.corr.ogtt.stomach <- read.csv("output/05_stomach_anxa_coexpression_ogtt_together.csv", header = T)

# Calculate the percentage of cells expressing ANXA10 (to sort the plot)
anxa10_expression <- FetchData(nowicki_healthy, vars = "ENSG00000109511")
anxa10_expression$cell_type <- Idents(nowicki_healthy)

# Calculate percentage of cells expressing ANXA10 for each cell type
percentage_expressing <- anxa10_expression %>% group_by(cell_type) %>% summarize(percent_expressing = sum(ENSG00000109511 > 0) / n() * 100)

# Sort cell types by percentage
sorted_cell_types <- percentage_expressing %>% arrange(desc(percent_expressing)) %>% pull(cell_type)

# Order cell types in the Seurat object
nowicki_healthy <- SetIdent(nowicki_healthy, value = factor(Idents(nowicki_healthy), levels = rev(sorted_cell_types)))

top25_corr <- head(anxa.corr.ogtt.stomach[order(-anxa.corr.ogtt.stomach$healthy.cor),]$gene, 25)
top25_corr_names <- head(anxa.corr.ogtt.stomach[order(-anxa.corr.ogtt.stomach$healthy.cor),]$hgnc_symbol.x, 25)

svg("graphics/05_stomach_top25_dotplot.svg", height=6, width=12)
DotPlot(nowicki_healthy, features = c(top25_corr,"ENSG00000107159","ENSG00000069764")) +
  scale_x_discrete(labels = c(top25_corr_names,"CA9","PLA2G10")) +
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  scale_size_continuous(breaks=c(20,40,60,80), limits = c(0,100)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
dev.off()

DotPlot(nowicki_healthy, features = genes,
        idents=c("Foveolar_Intermediate","Neck-Cells","Foveolar_Differentiated","Parietal")) +
  scale_x_discrete(labels = names) +
  scale_color_gradient(low = "lightblue", high = "darkblue", limits = c(-1.5, 1.5), breaks=c(-1,0,1)) +
  scale_size_continuous(breaks=c(20,40,60,80), limits = c(0,100)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10))