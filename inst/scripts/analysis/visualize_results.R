################################################################################
# Deconvolution Benchmark Dashboard for Nature Methods
# Creates a 4-panel figure: Pearson heatmap, Spearman heatmap, 
# Rank accuracy, Detection accuracy
################################################################################

# Load required packages
library(tidyverse)
library(scales)
library(cowplot)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(gridExtra)

# Define method colors (consistent across all panels)
method_colors <- c(
  "DeOPUS" = "#E64B35",
  "MuSiC" = "#00A087",
  "AutoGeneS" = "#4DBBD5",

  "CIBERSORT" = "#3C5488",
  "FARDEEP" = "#8491B4",
  "scaden" = "#7E6148",
  "AdRoit" = "#B15928"
)

# Method order for all plots
method_order <- names(method_colors)

################################################################################
# PRINT SUMMARY STATISTICS WITH STANDARD ERRORS
################################################################################

print_correlation_summary <- function(data, method_order = names(method_colors)) {
  
  # Filter data to include only methods in method_order
  data <- data %>%
    filter(method %in% method_order)
  
  # Calculate correlations by dataset first
 cor_data <- data %>%
    group_by(method, dataset) %>%
    summarise(
      cor_pearson = cor(ground_truth, predicted_value, method = "pearson", use = "complete.obs"),
      cor_spearman = cor(ground_truth, predicted_value, method = "spearman", use = "complete.obs"),
      .groups = "drop"
    )
  
  # Calculate mean and SE across datasets for each method
  summary_stats <- cor_data %>%
    group_by(method) %>%
    summarise(
      pearson_mean = mean(cor_pearson, na.rm = TRUE),
      pearson_se = sd(cor_pearson, na.rm = TRUE) / sqrt(n()),
      spearman_mean = mean(cor_spearman, na.rm = TRUE),
      spearman_se = sd(cor_spearman, na.rm = TRUE) / sqrt(n()),
      n_datasets = n(),
      .groups = "drop"
    ) %>%
    mutate(method = factor(method, levels = method_order)) %>%
    arrange(method)
  
  # Print formatted output
  cat("\n")
  cat("=" %>% rep(80) %>% paste(collapse = ""), "\n")
  cat("CORRELATION SUMMARY (Mean ± SE across datasets)\n")
  cat("=" %>% rep(80) %>% paste(collapse = ""), "\n\n")
  
  cat(sprintf("%-12s %20s %20s %10s\n", "Method", "Pearson (r)", "Spearman (ρ)", "N datasets"))
  cat("-" %>% rep(80) %>% paste(collapse = ""), "\n")
  
  for (i in 1:nrow(summary_stats)) {
    cat(sprintf("%-12s %12.3f ± %.3f %12.3f ± %.3f %10d\n",
                as.character(summary_stats$method[i]),
                summary_stats$pearson_mean[i],
                summary_stats$pearson_se[i],
                summary_stats$spearman_mean[i],
                summary_stats$spearman_se[i],
                summary_stats$n_datasets[i]))
  }
  
  cat("-" %>% rep(80) %>% paste(collapse = ""), "\n\n")
  
  # Return the summary table invisibly
  invisible(summary_stats)
}

print_rank_accuracy_summary <- function(data, method_order = names(method_colors)) {
  
  data <- data %>%
    filter(method %in% method_order)
  
  # Calculate rank accuracy by dataset
  ranked_data <- data %>%
    group_by(method, dataset, sample) %>%
    mutate(
      rank_truth = rank(-ground_truth, ties.method = "min"),
      rank_pred = rank(-predicted_value, ties.method = "min")
    ) %>%
    ungroup()
  
  top_k_accuracy <- ranked_data %>%
    group_by(method, dataset, sample) %>%
    summarise(
      top1_correct = as.numeric(cell_type[rank_pred == 1][1] == cell_type[rank_truth == 1][1]),
      top3_overlap = sum(cell_type[rank_pred <= 3] %in% cell_type[rank_truth <= 3]) / 
                     min(3, n_distinct(cell_type)),
      .groups = "drop"
    )
  
  # Summarize by method
  summary_stats <- top_k_accuracy %>%
    group_by(method) %>%
    summarise(
      top1_mean = mean(top1_correct, na.rm = TRUE),
      top1_se = sd(top1_correct, na.rm = TRUE) / sqrt(n()),
      top3_mean = mean(top3_overlap, na.rm = TRUE),
      top3_se = sd(top3_overlap, na.rm = TRUE) / sqrt(n()),
      n_samples = n(),
      .groups = "drop"
    ) %>%
    mutate(method = factor(method, levels = method_order)) %>%
    arrange(method)
  
  # Print formatted output
  cat("\n")
  cat("=" %>% rep(80) %>% paste(collapse = ""), "\n")
  cat("RANK-BASED ACCURACY SUMMARY (Mean ± SE)\n")
  cat("=" %>% rep(80) %>% paste(collapse = ""), "\n\n")
  
  cat(sprintf("%-12s %20s %20s %10s\n", "Method", "Top-1 Accuracy", "Top-3 Overlap", "N samples"))
  cat("-" %>% rep(80) %>% paste(collapse = ""), "\n")
  
  for (i in 1:nrow(summary_stats)) {
    cat(sprintf("%-12s %11.1f%% ± %.1f%% %11.1f%% ± %.1f%% %10d\n",
                as.character(summary_stats$method[i]),
                summary_stats$top1_mean[i] * 100,
                summary_stats$top1_se[i] * 100,
                summary_stats$top3_mean[i] * 100,
                summary_stats$top3_se[i] * 100,
                summary_stats$n_samples[i]))
  }
  
  cat("-" %>% rep(80) %>% paste(collapse = ""), "\n\n")
  
  invisible(summary_stats)
}

# Wrapper function to print all summaries
print_all_summaries <- function(data, method_order = names(method_colors)) {
  cor_summary <- print_correlation_summary(data, method_order)
  rank_summary <- print_rank_accuracy_summary(data, method_order)
  
  invisible(list(correlation = cor_summary, rank_accuracy = rank_summary))
}

################################################################################
# Theme for Nature Methods (clean, publication-ready)
################################################################################

theme_nature <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      # Text elements
      text = element_text(family = "sans", color = "black"),
      axis.text = element_text(size = base_size, color = "black"),
      axis.title = element_text(size = base_size + 1, face = "bold"),
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size, hjust = 0, color = "gray40"),
      
      # Legend
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size - 1),
      legend.key.size = unit(0.4, "cm"),
      legend.background = element_blank(),
      legend.position = "right",
      
      # Panel
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      
      # Margins
      plot.margin = margin(5, 5, 5, 5, "pt"),
      
      # Strip for facets
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold")
    )
}

################################################################################
# PANEL A & B: Correlation Heatmaps by Dataset
################################################################################

calculate_correlations_by_dataset <- function(data, method_order) {
  
  # Calculate Pearson and Spearman correlations per method-dataset combination
  cor_data <- data %>%
    group_by(method, dataset) %>%
    summarise(
      cor_pearson = cor(ground_truth, predicted_value, method = "pearson", use = "complete.obs"),
      cor_spearman = cor(ground_truth, predicted_value, method = "spearman", use = "complete.obs"),
      n_samples = n_distinct(sample),
      n_cell_types = n_distinct(cell_type),
      .groups = "drop"
    )
  
  # Factor methods in specified order (only include methods present in data)
  present_methods <- intersect(method_order, unique(cor_data$method))
  cor_data$method <- factor(cor_data$method, levels = present_methods)
  
  # Clean dataset names for display
  cor_data <- cor_data %>%
    mutate(dataset_clean = gsub("\\.rds$", "", dataset))
  
  return(cor_data)
}

create_correlation_heatmap_ggplot <- function(cor_data, cor_type = "pearson", 
                                               color_low = "#E8E8E8", 
                                               color_high = "#B2182B",
                                               add_average = TRUE) {
  
  cor_col <- ifelse(cor_type == "pearson", "cor_pearson", "cor_spearman")
  title_text <- ifelse(cor_type == "pearson", "Pearson Correlation", "Spearman Correlation")
  
  plot_data <- cor_data %>%
    mutate(cor_value = .data[[cor_col]])
  
  # Add average column if requested
  if (add_average) {
    avg_data <- plot_data %>%
      group_by(method) %>%
      summarise(
        cor_value = mean(cor_value, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(dataset_clean = "Average")
    
    plot_data <- bind_rows(plot_data, avg_data)
    
    # Set factor levels for x-axis with Average at the end
    dataset_levels <- c(unique(cor_data$dataset_clean), "Average")
    plot_data$dataset_clean <- factor(plot_data$dataset_clean, levels = dataset_levels)
  }
  
  # Compute text color based on correlation value
  plot_data <- plot_data %>%
    mutate(text_color = ifelse(cor_value > 0.6, "white", "black"))
  
  # Create heatmap
  p <- ggplot(plot_data, aes(x = dataset_clean, y = method, fill = cor_value)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", cor_value), color = text_color), 
              size = 2.5, show.legend = FALSE) +
    scale_color_identity() +
    scale_fill_gradient(
      low = color_low, 
      high = color_high,
      limits = c(0, 1),
      name = "Correlation",
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = c("0.00", "0.25", "0.50", "0.75", "1.00")
    ) +
    scale_y_discrete(limits = rev(levels(plot_data$method))) +
    labs(
      title = title_text,
      x = "Dataset",
      y = "Method"
    ) +
    theme_nature() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "right",
      panel.border = element_blank()
    )
  
  # Add vertical line to separate Average column
  if (add_average) {
    n_datasets <- length(unique(cor_data$dataset_clean))
    p <- p + geom_vline(xintercept = n_datasets + 0.5, linetype = "dashed", color = "gray40", linewidth = 0.5)
  }
  
  return(p)
}

################################################################################
# PANEL C: Rank-Based Accuracy by Dataset
################################################################################

calculate_rank_accuracy_by_dataset <- function(data, method_order) {
  
  # For each sample, rank cell types by ground truth and by prediction
  ranked_data <- data %>%
    group_by(method, dataset, sample) %>%
    mutate(
      rank_truth = rank(-ground_truth, ties.method = "min"),
      rank_pred = rank(-predicted_value, ties.method = "min"),
      n_types = n()
    ) %>%
    ungroup()
  
  # Calculate top-1 and top-3 accuracy per sample
  top_k_accuracy <- ranked_data %>%
    group_by(method, dataset, sample) %>%
    summarise(
      top1_correct = as.numeric(cell_type[rank_pred == 1][1] == cell_type[rank_truth == 1][1]),
      top3_overlap = sum(cell_type[rank_pred <= 3] %in% cell_type[rank_truth <= 3]) / 
                     min(3, n_distinct(cell_type)),
      rank_cor = cor(rank_truth, rank_pred, method = "spearman"),
      .groups = "drop"
    )
  
  # Summarize by method and dataset
  summary_by_dataset <- top_k_accuracy %>%
    group_by(method, dataset) %>%
    summarise(
      top1_acc = mean(top1_correct, na.rm = TRUE),
      top3_acc = mean(top3_overlap, na.rm = TRUE),
      rank_cor = mean(rank_cor, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(dataset_clean = gsub("\\.rds$", "", dataset))
  
  # Factor methods
  present_methods <- intersect(method_order, unique(summary_by_dataset$method))
  summary_by_dataset$method <- factor(summary_by_dataset$method, levels = present_methods)
  
  return(summary_by_dataset)
}

create_rank_accuracy_heatmap <- function(rank_data, metric = "top1_acc",
                                          color_low = "#F7F7F7", 
                                          color_high = "#2166AC") {
  
  metric_labels <- c(
    "top1_acc" = "Top-1 Accuracy",
    "top3_acc" = "Top-3 Overlap",
    "rank_cor" = "Rank Correlation"
  )
  
  plot_data <- rank_data %>%
    mutate(value = .data[[metric]])
  
  p <- ggplot(plot_data, aes(x = dataset_clean, y = method, fill = value)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", value)), 
              size = 2.5, color = ifelse(plot_data$value > 0.6, "white", "black")) +
    scale_fill_gradient(
      low = color_low, 
      high = color_high,
      limits = c(0, 1),
      name = metric_labels[metric],
      breaks = c(0, 0.25, 0.5, 0.75, 1)
    ) +
    scale_y_discrete(limits = rev(levels(plot_data$method))) +
    labs(
      title = "Rank-Based Accuracy",
      subtitle = "Correctly identifying most abundant cell types",
      x = "Dataset",
      y = "Method"
    ) +
    theme_nature() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "right",
      panel.border = element_blank()
    )
  
  return(p)
}

################################################################################
# PANEL D: Detection Accuracy by Dataset
################################################################################

calculate_detection_by_dataset <- function(data, method_order, threshold = 0.01) {
  
  # Binary classification: Is the cell type present (>threshold)?
  detection_data <- data %>%
    mutate(
      true_present = ground_truth > threshold,
      pred_present = predicted_value > threshold,
      correct = true_present == pred_present,
      tp = true_present & pred_present,
      fp = !true_present & pred_present,
      tn = !true_present & !pred_present,
      fn = true_present & !pred_present
    )
  
  # Calculate metrics per method and dataset
  detection_metrics <- detection_data %>%
    group_by(method, dataset) %>%
    summarise(
      accuracy = mean(correct),
      sensitivity = sum(tp) / max(sum(tp) + sum(fn), 1),
      specificity = sum(tn) / max(sum(tn) + sum(fp), 1),
      precision = sum(tp) / max(sum(tp) + sum(fp), 1),
      f1 = 2 * precision * sensitivity / max(precision + sensitivity, 0.001),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(dataset_clean = gsub("\\.rds$", "", dataset))
  
  # Factor methods
  present_methods <- intersect(method_order, unique(detection_metrics$method))
  detection_metrics$method <- factor(detection_metrics$method, levels = present_methods)
  
  return(detection_metrics)
}

create_detection_heatmap <- function(detection_data, metric = "f1",
                                      color_low = "#FEF0D9", 
                                      color_high = "#B30000") {
  
  metric_labels <- c(
    "accuracy" = "Accuracy",
    "sensitivity" = "Sensitivity",
    "specificity" = "Specificity",
    "precision" = "Precision",
    "f1" = "F1 Score"
  )
  
  plot_data <- detection_data %>%
    mutate(value = .data[[metric]])
  
  p <- ggplot(plot_data, aes(x = dataset_clean, y = method, fill = value)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", value)), 
              size = 2.5, color = ifelse(plot_data$value > 0.6, "white", "black")) +
    scale_fill_gradient(
      low = color_low, 
      high = color_high,
      limits = c(0, 1),
      name = metric_labels[metric],
      breaks = c(0, 0.25, 0.5, 0.75, 1)
    ) +
    scale_y_discrete(limits = rev(levels(plot_data$method))) +
    labs(
      title = "Detection Accuracy",
      subtitle = "Cell type presence (>1%) vs absence",
      x = "Dataset",
      y = "Method"
    ) +
    theme_nature() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "right",
      panel.border = element_blank()
    )
  
  return(p)
}

################################################################################
# MAIN DASHBOARD FUNCTION
################################################################################

create_benchmark_dashboard <- function(allDataToPlotF, 
                                        method_order = names(method_colors),
                                        detection_threshold = 0.01,
                                        output_file = "benchmark_dashboard.pdf",
                                        width = 12, height = 10) {
  
  # Filter data to include only methods in method_order
  data <- allDataToPlotF %>%
    filter(method %in% method_order)
  
  # Calculate all metrics
  message("Calculating correlations...")
  cor_data <- calculate_correlations_by_dataset(data, method_order)
  
  message("Calculating rank accuracy...")
  rank_data <- calculate_rank_accuracy_by_dataset(data, method_order)
  

  message("Calculating detection metrics...")
  detection_data <- calculate_detection_by_dataset(data, method_order, detection_threshold)
  
  # Create individual panels
  message("Creating Panel A: Pearson correlation heatmap...")
  panel_A <- create_correlation_heatmap_ggplot(
    cor_data, 
    cor_type = "pearson",
    color_low = "#E8E8E8",
    color_high = "#B2182B"
  ) + labs(title = NULL, subtitle = NULL)
  
  message("Creating Panel B: Spearman correlation heatmap...")
  panel_B <- create_correlation_heatmap_ggplot(
    cor_data, 
    cor_type = "spearman",
    color_low = "#E8E8E8",
    color_high = "#E66100"
  ) + labs(title = NULL, subtitle = NULL)
  
  message("Creating Panel C: Rank accuracy heatmap...")
  panel_C <- create_rank_accuracy_heatmap(
    rank_data, 
    metric = "top1_acc",
    color_low = "#F7F7F7",
    color_high = "#2166AC"
  ) + labs(title = NULL, subtitle = NULL)
  
  message("Creating Panel D: Detection accuracy heatmap...")
  panel_D <- create_detection_heatmap(
    detection_data, 
    metric = "f1",
    color_low = "#FEF0D9",
    color_high = "#7F3B08"
  ) + labs(title = NULL, subtitle = NULL)
  
  # Combine panels with labels
  message("Assembling dashboard...")
  dashboard <- plot_grid(
    panel_A + theme(legend.position = "right"),
    panel_B + theme(legend.position = "right"),
    panel_C + theme(legend.position = "right"),
    panel_D + theme(legend.position = "right"),
    labels = c("A", "B", "C", "D"),
    label_size = 12,
    label_fontface = "bold",
    ncol = 2,
    nrow = 2,
    align = "hv",
    axis = "tblr"
  )
  
  # Add main title
  title <- ggdraw() + 
    draw_label(
      "Deconvolution Benchmark: Real Datasets",
      fontface = "bold",
      size = 14,
      x = 0.5,
      hjust = 0.5
    )
  
  # Add subtitle with metric descriptions
  subtitle <- ggdraw() + 
    draw_label(
      "A: Pearson correlation  |  B: Spearman correlation  |  C: Top-1 rank accuracy  |  D: F1 detection score (threshold: 1%)",
      size = 9,
      x = 0.5,
      hjust = 0.5,
      color = "gray40"
    )
  
  final_plot <- plot_grid(
    title,
    subtitle,
    dashboard,
    ncol = 1,
    rel_heights = c(0.04, 0.03, 0.93)
  )
  
  # Save to file
  message(paste("Saving to", output_file, "..."))
  ggsave(
    output_file,
    final_plot,
    width = width,
    height = height,
    dpi = 300,
    bg = "white"
  )
  
  message("Done!")
  
  # Return all components for further customization
  return(list(
    dashboard = final_plot,
    panels = list(A = panel_A, B = panel_B, C = panel_C, D = panel_D),
    data = list(
      correlations = cor_data,
      rank_accuracy = rank_data,
      detection = detection_data
    )
  ))
}

################################################################################
# ALTERNATIVE: Summary Bar Plot Version for Panels C & D
################################################################################

create_summary_barplot_dashboard <- function(allDataToPlotF,
                                              method_order = names(method_colors),
                                              detection_threshold = 0.01,
                                              output_file = "benchmark_dashboard_bars.pdf",
                                              width = 12, height = 14) {
  
  data <- allDataToPlotF %>%
    filter(method %in% method_order)
  
  # Calculate metrics
  cor_data <- calculate_correlations_by_dataset(data, method_order)
  rank_data <- calculate_rank_accuracy_by_dataset(data, method_order)
  
  # Panel A: Grouped bar chart for Pearson and Spearman correlations
  cor_summary <- cor_data %>%
    group_by(method) %>%
    summarise(
      pearson_mean = mean(cor_pearson, na.rm = TRUE),
      pearson_se = sd(cor_pearson, na.rm = TRUE) / sqrt(n()),
      spearman_mean = mean(cor_spearman, na.rm = TRUE),
      spearman_se = sd(cor_spearman, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  cor_plot_data <- cor_summary %>%
    pivot_longer(
      cols = c(pearson_mean, spearman_mean),
      names_to = "metric",
      values_to = "correlation"
    ) %>%
    mutate(
      se = ifelse(metric == "pearson_mean", pearson_se, spearman_se),
      metric = factor(
        ifelse(metric == "pearson_mean", "Pearson", "Spearman"),
        levels = c("Pearson", "Spearman")
      )
    )
  
  panel_A <- ggplot(cor_plot_data, aes(x = method, y = correlation, fill = metric)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(
      aes(ymin = pmax(0, correlation - se), ymax = pmin(1, correlation + se)),
      position = position_dodge(width = 0.8),
      width = 0.2,
      linewidth = 0.3
    ) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50", linewidth = 0.3) +
    scale_fill_manual(values = c("Pearson" = "#B2182B", "Spearman" = "#E66100")) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(title = "Correlation", x = NULL, y = "Correlation", fill = "Metric") +
    theme_nature() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top",
      legend.justification = "left"
    )
  
  # Panel B: Rank accuracy bar plot
  rank_summary <- rank_data %>%
    group_by(method) %>%
    summarise(
      top1_acc = mean(top1_acc, na.rm = TRUE),
      top1_se = sd(top1_acc, na.rm = TRUE) / sqrt(n()),
      top3_acc = mean(top3_acc, na.rm = TRUE),
      top3_se = sd(top3_acc, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  rank_plot_data <- rank_summary %>%
    pivot_longer(
      cols = c(top1_acc, top3_acc),
      names_to = "metric",
      values_to = "accuracy"
    ) %>%
    mutate(
      se = ifelse(metric == "top1_acc", top1_se, top3_se),
      metric = factor(
        ifelse(metric == "top1_acc", "Top-1 Accuracy", "Top-3 Overlap"),
        levels = c("Top-1 Accuracy", "Top-3 Overlap")
      )
    )
  
  panel_B <- ggplot(rank_plot_data, aes(x = method, y = accuracy, fill = metric)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(
      aes(ymin = pmax(0, accuracy - se), ymax = pmin(1, accuracy + se)),
      position = position_dodge(width = 0.8), 
      width = 0.2,
      linewidth = 0.3
    ) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50", linewidth = 0.3) +
    scale_fill_manual(values = c("Top-1 Accuracy" = "#2166AC", "Top-3 Overlap" = "#67A9CF")) +
    scale_y_continuous(labels = percent_format(), limits = c(0, 1), expand = c(0, 0)) +
    labs(title = "Rank-Based Accuracy", x = NULL, y = "Accuracy", fill = "Metric") +
    theme_nature() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top",
      legend.justification = "left"
    )
  
  # Panel C: Pearson correlation heatmap
  panel_C <- create_correlation_heatmap_ggplot(
    cor_data, "pearson", "#E8E8E8", "#B2182B"
  ) + labs(title = "Pearson Correlation")
  
  # Panel D: Spearman correlation heatmap
  panel_D <- create_correlation_heatmap_ggplot(
    cor_data, "spearman", "#E8E8E8", "#E66100"
  ) + labs(title = "Spearman Correlation")
  
  # Assemble dashboard
  # Row 1: A | B
  row1 <- plot_grid(
    panel_A, panel_B,
    labels = c("A", "B"),
    label_size = 12,
    label_fontface = "bold",
    ncol = 2,
    align = "h"
  )
  
  # Row 2: C (full width)
  row2 <- plot_grid(
    panel_C,
    labels = c("C"),
    label_size = 12,
    label_fontface = "bold",
    ncol = 1
  )
  
  # Row 3: D (full width)
  row3 <- plot_grid(
    panel_D,
    labels = c("D"),
    label_size = 12,
    label_fontface = "bold",
    ncol = 1
  )
  
  # Combine all rows
  final_plot <- plot_grid(
    row1,
    row2,
    row3,
    ncol = 1,
    rel_heights = c(1, 1, 1)
  )
  
  ggsave(output_file, final_plot, width = width, height = height, dpi = 300, bg = "white")
  
  return(list(dashboard = final_plot, panels = list(A = panel_A, B = panel_B, C = panel_C, D = panel_D)))
}

################################################################################
# USAGE EXAMPLE
################################################################################

# # Load your data
# # allDataToPlotF <- readRDS("your_data.rds")
# 
# # Create the heatmap-based dashboard (all 4 panels as heatmaps)
# results <- create_benchmark_dashboard(
#   allDataToPlotF,
#   method_order = names(method_colors),
#   detection_threshold = 0.01,
#   output_file = "Figure_benchmark_real.pdf",
#   width = 12,
#   height = 10
# )
# 
# # Or create the bar plot version for panels C & D
# results_bars <- create_summary_barplot_dashboard(
#   allDataToPlotF,
#   method_order = names(method_colors),
#   detection_threshold = 0.01,
#   output_file = "Figure_benchmark_real_bars.pdf",
#   width = 14,
#   height = 10
# )
# 
# # Access individual panels for further customization
# # results$panels$A  # Pearson heatmap
# # results$panels$B  # Spearman heatmap
# # results$panels$C  # Rank accuracy
# # results$panels$D  # Detection accuracy
# 
# # Access calculated metrics
# # results$data$correlations
# # results$data$rank_accuracy
# # results$data$detection
