## https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/

plot_manhattan <- function(data, snps) {
  require(ggplot2)
  require(ggrepel)

  data_cum <- data %>%
    group_by(CHROM) %>%
    summarise(max_bp = max(GENPOS)) %>%
    mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>%
    select(CHROM, bp_add)

  gwas_data <- data %>%
    inner_join(data_cum, by = "CHROM") %>%
    mutate(bp_cum = GENPOS + bp_add)

  axis_set <- gwas_data %>%
    group_by(CHROM) %>%
    summarize(center = mean(bp_cum))

  ylim <- max(gwas_data$LOG10P) + 1

  # Draw plot
  plot <- ggplot(
    gwas_data,
    aes(x = bp_cum, y = LOG10P, color = as.factor(CHROM))
  ) +
    #geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") +
    geom_point(alpha = 0.75, size = 0.5) +
    scale_x_continuous(label = axis_set$CHROM, breaks = axis_set$center) +
    scale_y_continuous(expand = c(0, 0), limits = c(2, ylim)) +
    scale_color_manual(
      values = rep(c("#332288", "#6699CC"), unique(length(axis_set$CHROM)))
    ) +
    labs(x = NULL, y = "-log10P") +
    geom_text_repel(
      data = subset(gwas_data, ID %in% snps),
      aes(label = ID),
      size = 3
    ) +
    geom_hline(yintercept = -log10(5e-8), col = "#CC6677") +
    geom_hline(yintercept = -log10(1e-5), col = "#888888") +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      text = element_text(size = 12)
    )

  return(plot)
}
