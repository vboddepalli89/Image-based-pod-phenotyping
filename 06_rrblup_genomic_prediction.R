file_path <- "C:/Users/naresh89/Desktop/My projects/gwas/pod phenotyping/pod_seed_gwas_final/SVEN/seed.txt"

file.exists(file_path)  # This should return TRUE

# Read file content
file_content <- readLines(file_path)

# Compress and write gz file
gz_file_path <- paste0(file_path, ".gz")
gz_con <- gzfile(gz_file_path, "w")
writeLines(file_content, gz_con)
close(gz_con)


cat("File compressed to:", gz_file_path, "\n")

library(Matrix) # To create sparse matrices
getwd()
snp = read.table("seed.txt.gz",header=T)
dim(snp)
snp$X.Marker.
X = as.matrix(snp[,-1])
X = 1 - X
print(mean(X != 0))
mx = colMeans(X) > 0.05

# Convert to sparse format, first filetering out columns wiht MAF < 0.05
Z = Matrix(X[,mx],sparse = T)
rownames(Z) <- snp$X.Marker.
saveRDS(Z,file="zmat_sparse.rds")

###Now run the below code by modfiying the relavant areas with your file location and name
# Load required libraries
library(rrBLUP)
library(Matrix)

# Load data
Z <- readRDS("zmat_sparse.rds")
pheno <- read.csv("C:\\Users\\naresh89\\Desktop\\My projects\\gwas\\pod phenotyping\\Genomic prediction\\rrBLUP\\blups.csv _phenotype.csv")
pheno <- pheno[match(rownames(Z), pheno$Taxa), ]
stopifnot(all(rownames(Z) == pheno$Taxa))

# Traits to analyze
traits <- setdiff(colnames(pheno), "Taxa")

# 10-fold CV setup
set.seed(123)
K <- 10
folds <- sample(rep(1:K, length.out = nrow(Z)))

# Store all fold-wise results
cv_results <- list()

for (trait in traits) {
  cat("\nRunning 10-fold CV for:", trait, "\n")
  fold_metrics <- list()
  
  for (k in 1:K) {
    cat(" Fold", k, "\n")
    
    test_idx <- which(folds == k)
    train_idx <- setdiff(1:nrow(Z), test_idx)
    
    # Subset genotype and phenotype
    Z_train <- Z[train_idx, ]
    Z_test  <- Z[test_idx, ]
    y_train <- pheno[train_idx, trait]
    y_test  <- pheno[test_idx, trait]
    
    # Remove NA from y_train
    not_na <- !is.na(y_train)
    Z_train <- Z_train[not_na, ]
    y_train <- y_train[not_na]
    
    # MAF filtering
    maf <- colMeans(Z_train)
    keep_snps <- which(maf >= 0.05)
    Z_train <- Z_train[, keep_snps]
    Z_test  <- Z_test[, keep_snps]
    
    # Drop zero-variance SNPs
    non_zero_var <- apply(Z_train, 2, var) > 0
    Z_train <- Z_train[, non_zero_var]
    Z_test  <- Z_test[, non_zero_var]
    
    # Fit rrBLUP model
    model <- mixed.solve(y = y_train, Z = Z_train)
    
    # Align SNPs in Z_test to match model$u
    Z_test_full <- as.data.frame(as.matrix(Z_test))
    missing_snps <- setdiff(names(model$u), colnames(Z_test_full))
    for (snp in missing_snps) {
      Z_test_full[[snp]] <- 0
    }
    Z_test_aligned <- as.matrix(Z_test_full[, names(model$u), drop = FALSE])
    
    # Predict
    y_pred <- drop(Z_test_aligned %*% as.vector(model$u)) + as.numeric(model$beta)
    
    # Evaluate prediction
    valid <- !is.na(y_test) & !is.na(y_pred)
    y_obs <- y_test[valid]
    y_pred <- y_pred[valid]
    
    if (length(y_obs) > 0) {
      r <- cor(y_obs, y_pred)
      rsq <- r^2
      mse <- mean((y_obs - y_pred)^2)
      bias <- mean(y_pred - y_obs)
      
      fold_metrics[[k]] <- data.frame(
        Trait = trait,
        Fold = k,
        Method = "rrBLUP",
        Correlation = r,
        R_squared = rsq,
        MSE = mse,
        Bias = bias,
        N = length(y_obs)
      )
    }
  }
  
  if (length(fold_metrics) > 0) {
    cv_results[[trait]] <- do.call(rbind, fold_metrics)
  }
}

# Combine and save all fold-wise results
cv_df <- do.call(rbind, cv_results)
write.csv(cv_df, "rrBLUP_10fold_cv_foldwise_results.csv", row.names = FALSE)
cat("\n10-fold CV for rrBLUP completed. Fold-wise results saved to 'rrBLUP_10fold_cv_foldwise_results.csv'\n")

