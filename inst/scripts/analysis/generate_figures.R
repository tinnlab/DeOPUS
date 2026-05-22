################################################################################
# Generate Figures for Paper
# Produces publication-ready figures for Nature Methods submission
################################################################################

library(tidyverse)
library(scales)
library(cowplot)

# Source visualization functions. Prefer the installed-package location
# (system.file resolves after install_github / R CMD INSTALL); fall back to
# the in-repo path when running from a clone.
.viz_path <- system.file("scripts/analysis/visualize_results.R",
                         package = "DeOPUS")
if (!nzchar(.viz_path)) .viz_path <- "inst/scripts/analysis/visualize_results.R"
source(.viz_path)

################################################################################
# Configuration
################################################################################

# Input/output paths
RESULTS_FILE <- "results/real_benchmark/benchmark_results_aggregated.rds"
FIGURE_DIR <- "figures"

# Create output directory
if (!dir.exists(FIGURE_DIR)) {
  dir.create(FIGURE_DIR, recursive = TRUE)
}

# Method colors (consistent across all figures)
method_colors <- c(
  "DeOPUS" = "#E64B35",
  "MuSiC" = "#00A087",
  "AutoGeneS" = "#4DBBD5",
  "CIBERSORT" = "#3C5488",
  "FARDEEP" = "#8491B4",
  "scaden" = "#7E6148",
  "AdRoit" = "#B15928"
)

################################################################################
# Load Data
################################################################################

load_benchmark_data <- function(results_file = RESULTS_FILE) {
  
  if (!file.exists(results_file)) {
    stop("Results file not found: ", results_file)
  }
  
  data <- readRDS(results_file)
  
  # Ensure method is a factor with correct order
  data$method <- factor(data$method, levels = names(method_colors))
  
  # Clean dataset names
  data <- data %>%
    mutate(dataset = gsub("\\.rds$", "", dataset))
  
  return(data)
}

################################################################################
# Figure 1: Real Dataset Benchmark Dashboard
################################################################################

generate_figure_benchmark_real <- function(data, output_file = "Figure_benchmark_real.pdf") {
  
  message("Generating Figure: Real Dataset Benchmark")
  
  # Create dashboard
  results <- create_summary_barplot_dashboard(
    data,
    method_order = names(method_colors),
    detection_threshold = 0.01,
    output_file = file.path(FIGURE_DIR, output_file),
    width = 12,
    height = 14
  )
  
  message("Saved to: ", file.path(FIGURE_DIR, output_file))
  
  return(results)
}

################################################################################
# Figure 2: Detailed Performance Comparison
################################################################################

generate_figure_detailed_comparison <- function(data, output_file = "Figure_detailed_comparison.pdf") {
  
  message("Generating Figure: Detailed Comparison")
  
  # Calculate metrics
  cor_data <- data %>%
    group_by(method, dataset) %>%
    summarise(
      cor_pearson = cor(ground_truth, predicted_value, method = "pearson", use = "complete.obs"),
      cor_spearman = cor(ground_truth, predicted_value, method = "spearman", use = "complete.obs"),
      mse = mean((ground_truth - predicted_value)^2),
      mae = mean(abs(ground_truth - predicted_value)),
      .groups = "drop"
    )
  
  # Panel A: Pearson correlation distribution (violin + box)
  panel_A <- ggplot(cor_data, aes(x = method, y = cor_pearson, fill = method)) +
    geom_violin(alpha = 0.7, width = 0.8) +
    geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_fill_manual(values = method_colors) +
    scale_y_continuous(limits = c(-1, 1)) +
    labs(
      title = "Pearson Correlation Distribution",
      x = NULL,
      y = "Pearson Correlation"
    ) +
    theme_nature() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # Panel B: Spearman correlation distribution
  panel_B <- ggplot(cor_data, aes(x = method, y = cor_spearman, fill = method)) +
    geom_violin(alpha = 0.7, width = 0.8) +
    geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_fill_manual(values = method_colors) +
    scale_y_continuous(limits = c(-1, 1)) +
    labs(
      title = "Spearman Correlation Distribution",
      x = NULL,
      y = "Spearman Correlation"
    ) +
    theme_nature() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # Panel C: MSE distribution
  panel_C <- ggplot(cor_data, aes(x = method, y = mse, fill = method)) +
    geom_violin(alpha = 0.7, width = 0.8) +
    geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
    scale_fill_manual(values = method_colors) +
    scale_y_log10() +
    labs(
      title = "Mean Squared Error Distribution",
      x = NULL,
      y = "MSE (log scale)"
    ) +
    theme_nature() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # Panel D: Scatter plot - DeOPUS vs ground truth (example)
  decopus_data <- data %>%
    filter(method == "DeOPUS")
  
  panel_D <- ggplot(decopus_data, aes(x = ground_truth, y = predicted_value)) +
    geom_point(alpha = 0.3, size = 0.5, color = method_colors["DeOPUS"]) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5) +
    scale_x_continuous(limits = c(0, 1), labels = percent_format()) +
    scale_y_continuous(limits = c(0, 1), labels = percent_format()) +
    labs(
      title = "DeOPUS: Predicted vs Ground Truth",
      x = "Ground Truth Proportion",
      y = "Predicted Proportion"
    ) +
    theme_nature() +
    coord_fixed()
  
  # Combine panels
  figure <- plot_grid(
    panel_A, panel_B,
    panel_C, panel_D,
    labels = c("A", "B", "C", "D"),
    label_size = 12,
    label_fontface = "bold",
    ncol = 2,
    align = "hv"
  )
  
  # Save
  ggsave(
    file.path(FIGURE_DIR, output_file),
    figure,
    width = 10,
    height = 10,
    dpi = 300,
    bg = "white"
  )
  
  message("Saved to: ", file.path(FIGURE_DIR, output_file))
  
  return(figure)
}

################################################################################
# Print Summary Statistics
################################################################################

print_summary_for_paper <- function(data) {
  
  message("\n", paste(rep("=", 80), collapse = ""))
  message("SUMMARY STATISTICS FOR PAPER")
  message(paste(rep("=", 80), collapse = ""), "\n")
  
  # Print correlation summary
  print_correlation_summary(data)
  
  # Print rank accuracy summary
  print_rank_accuracy_summary(data)
  
  # Additional statistics
  cor_data <- data %>%
    group_by(method, dataset) %>%
    summarise(
      cor_pearson = cor(ground_truth, predicted_value, method = "pearson", use = "complete.obs"),
      .groups = "drop"
    )
  
  # Count positive correlations
  positive_counts <- cor_data %>%
    group_by(method) %>%
    summarise(
      n_positive = sum(cor_pearson > 0, na.rm = TRUE),
      n_total = n(),
      pct_positive = n_positive / n_total * 100
    )
  
  message("\n", paste(rep("=", 80), collapse = ""))
  message("POSITIVE CORRELATION COUNTS")
  message(paste(rep("=", 80), collapse = ""), "\n")
  
  print(positive_counts)
  
  return(invisible(list(correlations = cor_data, positive_counts = positive_counts)))
}

################################################################################
# Main Execution
################################################################################

if (sys.nframe() == 0) {
  
  # Load data
  data <- load_benchmark_data()
  
  # Generate figures
  fig1 <- generate_figure_benchmark_real(data)
  fig2 <- generate_figure_detailed_comparison(data)
  
  # Print statistics
  stats <- print_summary_for_paper(data)
  
  message("\n=== Figure Generation Complete ===")
}
