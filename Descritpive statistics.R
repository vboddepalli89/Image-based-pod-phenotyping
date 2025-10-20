# Pod Phenotype Data Analysis in R
# ==================================
# 
# Comprehensive analysis script for pod phenotyping data including:
# 1. Data import and validation
# 2. Descriptive statistics 
# 3. Data visualization
# 4. BLUEs estimation using mixed models
# 5. Heritability calculation
# 6. GAPIT-compatible output formatting
#
# Author: Generated for pod GWAS analysis
# Date: 2025-09-04

# Required packages
required_packages <- c("readxl", "dplyr", "ggplot2", "lme4", "emmeans", 
                      "car", "reshape2", "corrplot", "VIM", "plotly", 
                      "gridExtra", "broom.mixed", "performance")

# Install missing packages
missing_packages <- required_packages[!required_packages %in% installed.packages()[,"Package"]]
if(length(missing_packages)) {
  cat("Installing missing packages:", paste(missing_packages, collapse=", "), "\n")
  install.packages(missing_packages, repos="https://cran.r-project.org")
}

# Load libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(lme4)
library(emmeans)
library(car)
library(reshape2)
library(corrplot)
library(VIM)
library(gridExtra)
library(broom.mixed)
library(performance)

# Set working directory and create output directories
setwd("C:/Users/naresh89/Desktop/pod_gwas_analysis")
dir.create("pod_analysis_results", showWarnings = FALSE)
dir.create("pod_analysis_results/plots", showWarnings = FALSE)
dir.create("pod_analysis_results/tables", showWarnings = FALSE)
dir.create("pod_analysis_results/blues", showWarnings = FALSE)
dir.create("pod_analysis_results/reports", showWarnings = FALSE)

# Data directory
data_dir <- "C:/Users/naresh89/Desktop/poddata"

cat("=== Pod Phenotype Data Analysis ===\n")
cat("Starting comprehensive analysis...\n\n")

# Function to load and clean data
load_phenotype_data <- function() {
  cat("Loading phenotype data files...\n")
  
  # Define file mappings
  file_list <- list(
    pod_length = "pod.length_2022-23.xlsx",
    pod_curvature_image = "pod.curvature.image_2022_23.xlsx", 
    pod_curvature_manual = "pod.curvature.manual_2022_23.xlsx",
    pod_width_image = "pod.width.image_2022_23.xlsx",
    seed_per_pod = "seed.per.pod_2022_23.xlsx"
  )
  
  data_list <- list()
  
  for (trait_name in names(file_list)) {
    file_path <- file.path(data_dir, file_list[[trait_name]])
    
    if (file.exists(file_path)) {
      cat("  Loading", file_list[[trait_name]], "...\n")
      
      tryCatch({
        # Read Excel file
        df <- read_excel(file_path)
        
        # Standardize column names
        names(df) <- tolower(trimws(names(df)))
        
        # Check required columns
        required_cols <- c("env", "block", "line")
        missing_cols <- required_cols[!required_cols %in% names(df)]
        
        if (length(missing_cols) > 0) {
          cat("    Warning: Missing columns", paste(missing_cols, collapse=", "), "in", file_list[[trait_name]], "\n")
          next
        }
        
        # Clean data
        df <- df %>%
          filter(!is.na(env), !is.na(block), !is.na(line)) %>%
          mutate(
            env = as.character(env),
            block = as.character(block),
            line = as.character(line)
          )
        
        data_list[[trait_name]] <- df
        cat("    Loaded", nrow(df), "observations for", trait_name, "\n")
        
      }, error = function(e) {
        cat("    Error loading", file_list[[trait_name]], ":", e$message, "\n")
      })
    } else {
      cat("  File not found:", file_list[[trait_name]], "\n")
    }
  }
  
  cat("\nSuccessfully loaded", length(data_list), "trait datasets\n\n")
  return(data_list)
}

# Function to identify trait columns
identify_trait_columns <- function(df) {
  system_cols <- c("env", "block", "line", "rep", "replicate")
  trait_cols <- names(df)[!tolower(names(df)) %in% system_cols]
  return(trait_cols)
}

# Function to calculate descriptive statistics
calculate_descriptive_stats <- function(data_list) {
  cat("Calculating descriptive statistics...\n")
  
  desc_stats_list <- list()
  
  for (trait_name in names(data_list)) {
    df <- data_list[[trait_name]]
    cat("  Processing", trait_name, "...\n")
    
    trait_cols <- identify_trait_columns(df)
    
    for (trait_col in trait_cols) {
      if (is.numeric(df[[trait_col]])) {
        # Remove missing values for analysis
        values <- df[[trait_col]][!is.na(df[[trait_col]])]
        
        if (length(values) == 0) next
        
        # Basic descriptive statistics
        desc_stats <- data.frame(
          trait = paste0(trait_name, "_", trait_col),
          n_obs = length(values),
          n_missing = sum(is.na(df[[trait_col]])),
          mean = mean(values),
          std = sd(values),
          variance = var(values),
          min = min(values),
          max = max(values),
          median = median(values),
          q25 = quantile(values, 0.25),
          q75 = quantile(values, 0.75),
          iqr = IQR(values),
          skewness = moments::skewness(values),
          kurtosis = moments::kurtosis(values),
          cv = (sd(values) / mean(values)) * 100,
          stringsAsFactors = FALSE
        )
        
        # Normality test
        if (length(values) >= 3 && length(values) <= 5000) {
          shapiro_test <- shapiro.test(values)
          desc_stats$shapiro_stat <- shapiro_test$statistic
          desc_stats$shapiro_p <- shapiro_test$p.value
          desc_stats$normal_shapiro <- shapiro_test$p.value > 0.05
        } else {
          desc_stats$shapiro_stat <- NA
          desc_stats$shapiro_p <- NA
          desc_stats$normal_shapiro <- NA
        }
        
        # Outlier detection (IQR method)
        q1 <- quantile(values, 0.25)
        q3 <- quantile(values, 0.75)
        iqr <- q3 - q1
        lower_bound <- q1 - 1.5 * iqr
        upper_bound <- q3 + 1.5 * iqr
        outliers <- values[values < lower_bound | values > upper_bound]
        desc_stats$n_outliers <- length(outliers)
        desc_stats$outlier_percentage <- (length(outliers) / length(values)) * 100
        
        # Environment and genotype counts
        desc_stats$n_environments <- length(unique(df$env))
        desc_stats$n_genotypes <- length(unique(df$line))
        desc_stats$environments <- paste(unique(df$env), collapse = ", ")
        
        # Store results
        key <- paste0(trait_name, "_", trait_col)
        desc_stats_list[[key]] <- desc_stats
      }
    }
  }
  
  cat("  Completed descriptive statistics for", length(desc_stats_list), "trait measurements\n\n")
  return(desc_stats_list)
}

# Function to create visualization plots
create_descriptive_plots <- function(data_list) {
  cat("Creating visualization plots...\n")
  
  for (trait_name in names(data_list)) {
    df <- data_list[[trait_name]]
    trait_cols <- identify_trait_columns(df)
    
    for (trait_col in trait_cols) {
      if (!is.numeric(df[[trait_col]])) next
      
      trait_key <- paste0(trait_name, "_", trait_col)
      values <- df[[trait_col]][!is.na(df[[trait_col]])]
      
      if (length(values) == 0) next
      
      # Create plots
      plot_list <- list()
      
      # 1. Histogram with normal curve
      p1 <- ggplot(data.frame(values = values), aes(x = values)) +
        geom_histogram(aes(y = after_stat(density)), bins = 30, alpha = 0.7, fill = "skyblue", color = "black") +
        stat_function(fun = dnorm, args = list(mean = mean(values), sd = sd(values)), color = "red", size = 1) +
        labs(title = "Distribution with Normal Curve", x = "Value", y = "Density") +
        theme_minimal()
      
      # 2. Q-Q plot
      p2 <- ggplot(data.frame(values = values), aes(sample = values)) +
        stat_qq() + stat_qq_line(color = "red") +
        labs(title = "Q-Q Plot (Normality Check)", x = "Theoretical Quantiles", y = "Sample Quantiles") +
        theme_minimal()
      
      # 3. Box plot by environment
      df_plot <- df[!is.na(df[[trait_col]]), c("env", trait_col)]
      names(df_plot)[2] <- "value"
      p3 <- ggplot(df_plot, aes(x = env, y = value)) +
        geom_boxplot(fill = "lightblue", alpha = 0.7) +
        labs(title = "Distribution by Environment", x = "Environment", y = "Value") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      # 4. Box plot by block
      df_plot2 <- df[!is.na(df[[trait_col]]), c("block", trait_col)]
      names(df_plot2)[2] <- "value"
      p4 <- ggplot(df_plot2, aes(x = block, y = value)) +
        geom_boxplot(fill = "lightgreen", alpha = 0.7) +
        labs(title = "Distribution by Block", x = "Block", y = "Value") +
        theme_minimal()
      
      # 5. Environment means
      env_means <- df %>%
        group_by(env) %>%
        summarise(mean_value = mean(.data[[trait_col]], na.rm = TRUE), .groups = "drop")
      
      p5 <- ggplot(env_means, aes(x = env, y = mean_value)) +
        geom_point(size = 3, color = "blue") +
        geom_text(aes(label = round(mean_value, 2)), vjust = -1) +
        labs(title = "Environment Means", x = "Environment", y = "Mean Value") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      # 6. Residual plot
      df_grouped <- df %>%
        group_by(env, line) %>%
        summarise(mean_value = mean(.data[[trait_col]], na.rm = TRUE), .groups = "drop") %>%
        filter(!is.na(mean_value))
      
      if (nrow(df_grouped) > 10) {
        overall_mean <- mean(df_grouped$mean_value)
        df_grouped$residuals <- df_grouped$mean_value - overall_mean
        
        p6 <- ggplot(df_grouped, aes(x = mean_value, y = residuals)) +
          geom_point(alpha = 0.6) +
          geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
          labs(title = "Residual Plot", x = "Fitted Values", y = "Residuals") +
          theme_minimal()
      } else {
        p6 <- ggplot() + labs(title = "Insufficient data for residual plot") + theme_minimal()
      }
      
      # Combine plots
      combined_plot <- grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3, 
                                   top = paste("Data Visualization:", trait_key))
      
      # Save plot
      plot_file <- file.path("pod_analysis_results/plots", paste0("descriptive_", trait_key, ".png"))
      ggsave(plot_file, combined_plot, width = 18, height = 12, dpi = 300)
      
      cat("  Saved plot:", plot_file, "\n")
    }
  }
  cat("\n")
}

# Function to estimate BLUEs using mixed models
estimate_blues <- function(data_list) {
  cat("Estimating BLUEs using mixed models...\n")
  
  blues_results <- list()
  heritability_results <- list()
  
  for (trait_name in names(data_list)) {
    df <- data_list[[trait_name]]
    trait_cols <- identify_trait_columns(df)
    
    for (trait_col in trait_cols) {
      if (!is.numeric(df[[trait_col]])) next
      
      trait_key <- paste0(trait_name, "_", trait_col)
      cat("  Processing", trait_key, "...\n")
      
      # Prepare data
      df_analysis <- df[, c("env", "block", "line", trait_col)] %>%
        filter(!is.na(.data[[trait_col]])) %>%
        mutate(
          env = as.factor(env),
          block = as.factor(block),
          line = as.factor(line),
          env_block = as.factor(paste0(env, "_", block))
        )
      
      names(df_analysis)[4] <- "trait_value"
      
      if (nrow(df_analysis) < 10) {
        cat("    Insufficient data for", trait_key, "\n")
        next
      }
      
      tryCatch({
        # Fit mixed model: trait_value ~ env + line + (1|env_block)
        model <- lmer(trait_value ~ env + line + (1|env_block), data = df_analysis, REML = TRUE)
        
        # Extract BLUEs using emmeans
        blues_em <- emmeans(model, ~ line)
        blues_df <- as.data.frame(blues_em) %>%
          select(line, emmean) %>%
          rename(BLUE = emmean) %>%
          arrange(line)
        
        # Calculate variance components
        var_comp <- as.data.frame(VarCorr(model))
        genotype_var <- var_comp[var_comp$grp == "line", "vcov"]
        block_var <- var_comp[var_comp$grp == "env_block", "vcov"]
        residual_var <- var_comp[var_comp$grp == "Residual", "vcov"]
        
        total_var <- genotype_var + block_var + residual_var
        heritability <- genotype_var / total_var
        
        # Model fit statistics
        model_summary <- summary(model)
        aic_val <- AIC(model)
        bic_val <- BIC(model)
        
        # Store results
        blues_results[[trait_key]] <- list(
          blues = blues_df,
          model = model,
          variance_components = list(
            genotype_var = genotype_var,
            block_var = block_var,
            residual_var = residual_var,
            total_var = total_var,
            heritability = heritability
          ),
          aic = aic_val,
          bic = bic_val,
          n_obs = nrow(df_analysis),
          n_genotypes = length(unique(df_analysis$line)),
          n_environments = length(unique(df_analysis$env))
        )
        
        heritability_results[[trait_key]] <- heritability
        
        cat("    BLUEs estimated for", nrow(blues_df), "genotypes\n")
        cat("    Heritability:", round(heritability, 3), "\n")
        
      }, error = function(e) {
        cat("    Error fitting model for", trait_key, ":", e$message, "\n")
      })
    }
  }
  
  cat("\nBLUEs estimation completed for", length(blues_results), "traits\n\n")
  return(list(blues_results = blues_results, heritability = heritability_results))
}

# Function to save results
save_results <- function(desc_stats_list, blues_results, heritability_results) {
  cat("Saving analysis results...\n")
  
  # 1. Save descriptive statistics
  desc_stats_df <- do.call(rbind, desc_stats_list)
  write.csv(desc_stats_df, "pod_analysis_results/tables/descriptive_statistics.csv", row.names = FALSE)
  cat("  Saved descriptive statistics\n")
  
  # 2. Save BLUEs for each trait
  for (trait_key in names(blues_results)) {
    blues_df <- blues_results[[trait_key]]$blues
    write.csv(blues_df, paste0("pod_analysis_results/blues/blues_", trait_key, ".csv"), row.names = FALSE)
  }
  cat("  Saved individual BLUEs files\n")
  
  # 3. Save variance components and heritability
  herit_data <- data.frame(
    trait = names(blues_results),
    genotype_variance = sapply(blues_results, function(x) x$variance_components$genotype_var),
    block_variance = sapply(blues_results, function(x) x$variance_components$block_var),
    residual_variance = sapply(blues_results, function(x) x$variance_components$residual_var),
    total_variance = sapply(blues_results, function(x) x$variance_components$total_var),
    heritability = sapply(blues_results, function(x) x$variance_components$heritability),
    n_genotypes = sapply(blues_results, function(x) x$n_genotypes),
    n_environments = sapply(blues_results, function(x) x$n_environments),
    n_observations = sapply(blues_results, function(x) x$n_obs),
    aic = sapply(blues_results, function(x) x$aic),
    bic = sapply(blues_results, function(x) x$bic),
    stringsAsFactors = FALSE
  )
  
  write.csv(herit_data, "pod_analysis_results/tables/heritability_variance_components.csv", row.names = FALSE)
  cat("  Saved heritability results\n")
  
  # 4. Create GAPIT-compatible phenotype file
  all_blues <- list()
  genotype_list <- c()
  
  for (trait_key in names(blues_results)) {
    blues_df <- blues_results[[trait_key]]$blues
    all_blues[[trait_key]] <- setNames(blues_df$BLUE, blues_df$line)
    genotype_list <- c(genotype_list, blues_df$line)
  }
  
  genotype_list <- sort(unique(genotype_list))
  gapit_df <- data.frame(Taxa = genotype_list, stringsAsFactors = FALSE)
  
  for (trait_key in names(all_blues)) {
    gapit_df[[trait_key]] <- all_blues[[trait_key]][genotype_list]
  }
  
  write.csv(gapit_df, "pod_analysis_results/blues/GAPIT_phenotypes.csv", row.names = FALSE)
  cat("  Saved GAPIT-compatible phenotypes\n")
  
  return(gapit_df)
}

