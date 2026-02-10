setwd("C:/Users/naresh89/Desktop/My projects/gwas/pod phenotyping/pod_phenotyping_manuscript/genomic prediction/weilv9/blups")

library(rrBLUP)
library(Matrix)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(stringr)

# Load Phenotypes (Use check.names=FALSE to preserve underscores)
phenotypes <- read.table("pod_blups.txt", header = TRUE, sep = "\t", check.names = FALSE)

# Load Genotypes (17k SNPs)
genotypes <- read.table("pod_numerical.txt", header = TRUE, sep = "\t", check.names = FALSE)

# 2. Data Processing
# ------------------------------------------------------------------------------
# Fix column name for merging
colnames(genotypes)[1] <- "Taxa"

# Prepare Genotype Matrix
X <- as.matrix(genotypes[, -1])
X <- 1 - X # Reverse coding if needed
rownames(X) <- genotypes$Taxa

# Filter MAF > 0.05
maf <- colMeans(X)
X_filtered <- X[, maf > 0.05]
Z <- Matrix(X_filtered, sparse = TRUE)

# Match Phenotypes to Genotypes
pheno_matched <- phenotypes[match(rownames(Z), phenotypes$Taxa), ]
traits <- setdiff(colnames(pheno_matched), "Taxa")

# 3. Run 10-Fold CV with Coincidence Calculation
# ------------------------------------------------------------------------------
set.seed(123)
K <- 10
folds <- sample(rep(1:K, length.out = nrow(Z)))

# Store results
cv_results <- list()

for (trait in traits) {
  cat("\nRunning analysis for:", trait, "\n")
  fold_metrics <- list()
  
  for (k in 1:K) {
    # Define Test/Train sets
    test_idx <- which(folds == k)
    train_idx <- setdiff(1:nrow(Z), test_idx)
    
    Z_train <- Z[train_idx, ]; Z_test <- Z[test_idx, ]
    y_train <- pheno_matched[train_idx, trait]; y_test <- pheno_matched[test_idx, trait]
    
    # Remove NAs in training
    not_na <- !is.na(y_train)
    Z_train <- Z_train[not_na, ]; y_train <- y_train[not_na]
    
    if (length(y_train) == 0) next
    
    # MAF & Variance Filter on Train
    maf_train <- colMeans(Z_train)
    valid_snps <- which(maf_train >= 0.05 & apply(Z_train, 2, var) > 0)
    Z_train <- Z_train[, valid_snps]; Z_test <- Z_test[, valid_snps]
    
    # Fit Model
    tryCatch({
      model <- mixed.solve(y = y_train, Z = Z_train)
      
      # Align Test Data
      u_names <- names(model$u)
      Z_test_aligned <- matrix(0, nrow = nrow(Z_test), ncol = length(u_names))
      colnames(Z_test_aligned) <- u_names
      common <- intersect(colnames(Z_test), u_names)
      if(length(common) > 0) Z_test_aligned[, common] <- as.matrix(Z_test[, common])
      
      # Predict
      y_pred <- drop(Z_test_aligned %*% as.vector(model$u)) + as.numeric(model$beta)
      
      # --- CALCULATION BLOCK ---
      valid <- !is.na(y_test) & !is.na(y_pred)
      y_obs_v <- y_test[valid]
      y_pred_v <- y_pred[valid]
      
      if (length(y_obs_v) > 0) {
        # 1. Accuracy
        r <- cor(y_obs_v, y_pred_v)
        
        # 2. Selection Coincidence Index (Top 10%)
        # Determine how many lines to select (e.g., Top 10 lines or Top 20%)
        n_select <- max(1, round(0.2 * length(y_obs_v))) # Select top 20% of the fold
        
        # Get indices of top performers
        top_obs <- order(y_obs_v, decreasing = TRUE)[1:n_select]
        top_pred <- order(y_pred_v, decreasing = TRUE)[1:n_select]
        
        # Count intersection
        common_lines <- length(intersect(top_obs, top_pred))
        coincidence_index <- (common_lines / n_select) * 100
        
        fold_metrics[[k]] <- data.frame(
          Trait = trait,
          Fold = k,
          Correlation = r,
          Coincidence_Index = coincidence_index,
          Selected_N = n_select,
          Common_N = common_lines
        )
      }
    }, error = function(e) cat("Error in fold", k, ":", conditionMessage(e), "\n"))
  }
  if (length(fold_metrics) > 0) cv_results[[trait]] <- do.call(rbind, fold_metrics)
}

# Compile Results
cv_df <- do.call(rbind, cv_results)
write.csv(cv_df, "rrBLUP_Full_Results_with_Coincidence.csv", row.names = FALSE)
# 4. CALCULATE GEBVs FOR ALL LINES (NEW SECTION)
# ------------------------------------------------------------------------------
cat("\nCalculating Final GEBVs for all lines using full dataset...\n")

# Initialize a dataframe with Taxa names
gebv_results <- data.frame(Taxa = rownames(Z))

for (trait in traits) {
  cat("  Fitting full model for:", trait, "\n")
  
  # Get phenotype vector
  y_full <- pheno_matched[, trait]
  
  # Filter NAs for model fitting (but we predict for everyone)
  not_na <- !is.na(y_full)
  y_fit <- y_full[not_na]
  Z_fit <- Z[not_na, ]
  
  # Fit Model on all available data
  # Note: We use the Z matrix that was already MAF filtered (>0.05) in step 2
  model_full <- mixed.solve(y = y_fit, Z = Z_fit)
  
  # Calculate GEBV for ALL lines in Z (even those with missing phenotypes)
  # Formula: GEBV = Intercept (beta) + (Genotype Matrix * Marker Effects)
  marker_effects <- model_full$u
  intercept <- as.numeric(model_full$beta)
  
  # Ensure Z matches marker effects (should be exact since we used Z subset to fit)
  gebv_values <- intercept + drop(as.matrix(Z) %*% marker_effects)
  
  # Add to results dataframe
  gebv_results[[paste0("GEBV_", trait)]] <- gebv_values
}

# Save GEBVs to CSV
write.csv(gebv_results, "GEBV_All_Lines_Predicted.csv", row.names = FALSE)
cat("GEBV calculation complete. Saved to 'GEBV_All_Lines_Predicted.csv'\n")

# 4. Summarize Data for Plotting
# ------------------------------------------------------------------------------
summary_stats <- cv_df %>%
  group_by(Trait) %>%
  summarise(
    Mean_Accuracy = mean(Correlation, na.rm = TRUE),
    SD_Accuracy = sd(Correlation, na.rm = TRUE),
    Mean_Coincidence = mean(Coincidence_Index, na.rm = TRUE),
    Mean_Common = mean(Common_N, na.rm = TRUE),
    .groups = 'drop'
  )
# Save to your current working directory
write.csv(summary_stats, "Mungbean_GP_Summary_Stats.csv", row.names = FALSE)
####plot
             trait_mapping <- data.frame(
  original = c("pod curvature_Image", "pod curvature_manual", 
               "pod length_Image", "pod length_manual", 
               "podwidthmm", 
               "seed per pod_Image", "seed per pod_manual"),
  
  display_name = c("Pod Curvature\n(Image)", "Pod Curvature\n(Manual)",
                   "Pod Length\n(Image)", "Pod Length\n(Manual)",
                   "Pod Width\n(Image)",
                   "Seeds per Pod\n(Image)", "Seeds per Pod\n(Manual)"),
  
  short_name = c("PC_Image", "PC_Manual", "PL_Image", "PL_Manual", 
                 "PW_Image", "SPP_Image", "SPP_Manual"),
  
  method = c("Image", "Manual", "Image", "Manual", "Image", "Image", "Manual")
)

# Merge mapping with results
# We use inner_join to ensure we only plot traits that successfully mapped
plot_data <- summary_stats %>%
  inner_join(trait_mapping, by = c("Trait" = "original"))

# Check if data merged correctly (If this prints 0 rows, the names are still off)
print(paste("Rows to plot:", nrow(plot_data)))

# Set Order
trait_order <- c("PC_Image", "PC_Manual", "PL_Image", "PL_Manual", "PW_Image", "SPP_Image", "SPP_Manual")
plot_data$short_name <- factor(plot_data$short_name, levels = trait_order)

# Colors
image_color <- "#2E86AB"
manual_color <- "#E63946"

# Panel A: Prediction Accuracy
panel_a <- ggplot(plot_data, aes(x = short_name, y = Mean_Accuracy, fill = method)) +
  geom_col(alpha = 0.8, width = 0.7) +
  geom_errorbar(aes(ymin = Mean_Accuracy - SD_Accuracy, ymax = Mean_Accuracy + SD_Accuracy), width = 0.2) +
  geom_text(aes(label = sprintf("%.2f", Mean_Accuracy)), vjust = -0.5, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c("Image" = image_color, "Manual" = manual_color)) +
  labs(title = "(A) Genomic Prediction Accuracies", x = "", y = "Prediction Accuracy (r)") +
  theme_minimal() +
  theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1))

# Panel B: Selection Coincidence
panel_b <- ggplot(plot_data, aes(x = short_name, y = Mean_Coincidence, fill = method)) +
  geom_col(alpha = 0.8, width = 0.7) +
  geom_text(aes(label = sprintf("%.1f", Mean_Coincidence)), vjust = -0.3, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c("Image" = image_color, "Manual" = manual_color)) +
  labs(title = "(B) Selection Coincidence Indices (Top 20%)", x = "Traits", y = "Coincidence Index (%)") +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

# Save as PNG
png("Final_Publication_Figure.png", width = 12, height = 8, units = "in", res = 300)
grid.arrange(panel_a, panel_b, nrow = 2)
dev.off()

print("Plot saved as 'Final_Publication_Figure.png'")
