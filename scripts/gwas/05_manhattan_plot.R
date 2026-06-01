library(data.table)
library(ggplot2)
library(ggrepel)
library(dplyr)

file_pattern <- "*anxa10_int.regenie"
files <- list.files(path = ".", pattern = "anxa10_int.regenie$")

# Read and combine data
combined_data <- rbindlist(lapply(files, function(file) {
  data <- fread(file)
  return(data)
}))

combined_data$P <- 10**-combined_data$LOG10P
filt_data <- combined_data[combined_data$LOG10P > 2, ]
filt_data <- filt_data %>% arrange(CHROM, GENPOS)

## function to plot manhattan
source("functions/plot_manhattan.R")
pdf("graphics/anxa10_manhattan.pdf", width = 10, height = 5)
plot_manhattan(
  filt_data,
  c("rs139558368", "rs2990223", "rs760077", "rs4600882", "rs6914749")
)
dev.off()
