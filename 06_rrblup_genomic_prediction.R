# Define the file path
file_path <- "C:\\Users\\naresh89\\Desktop\\My projects\\gwas\\pod phenotyping\\pod_phenotyping_manuscript\\GWAS\\SVEN\\pod.txt"
gz_file_path <- paste0(file_path, ".gz")

# Use gzfile to compress
file_content <- readLines(file_path)
gz_con <- gzfile(gz_file_path, "w")
writeLines(file_content, gz_con)
close(gz_con)

cat("File compressed to:", gz_file_path, "\n")

library(Matrix) # To create sparse matrices
getwd()
# MLScoreTraits_imputed.txt.gz is obtained from TASSEL after
# first running a kNN imputation and then a simple average
# imputation. Saved as numerical genotype. Then compressed
# using gzip to save space on disk.

# read the numeric genotype output from TASSEL
snp = read.table("pod.txt.gz",header=T)
dim(snp)
snp$X.Marker.
# The first column contains the name of the plants. Store the SNP values in X:
X = as.matrix(snp[,-1])


# Since major alleles are coded as 1 and minor alleles are 0
# reverse the coding to have more zeros (which can be ignored)
# This also does not affect the model except for a sign change of the coefficients.
X = 1 - X

# Check: the following should be much less than 0.5 for sparse format to be useful
# If you find the answer is NA, there's probably still SNPs with missing alleles values.
print(mean(X != 0))

# Now find columns with MAF > 0.05
mx = colMeans(X) > 0.05

# Convert to sparse format, first filetering out columns wiht MAF < 0.05
Z = Matrix(X[,mx],sparse = T)

# assign the row names for later reference....
# we need these to make sure the phenotype values are correctly
# matched with the genotypes.
rownames(Z) <- snp$X.Marker.
saveRDS(Z,file="zmat_sparse.rds")

##################
# Check out SNP_to_X.r if you have not already
#install.packages("bravo")
#install.packages("doParallel")
library(Matrix)
library(bravo)
library(doParallel) # For parallelizing on Linux/Unix machines

# We want to run the GWAS for each environment/phenotype parallely:
# How many cores do we want to use?
# my system has 10 cores(20 threads), but we have only four environments.
# So need only 4 cores.
# Note: detectCores() give the number of available cores.
# Note: on nova then you have to request 4 CPU cores as well.
registerDoMC(cores = 1)

# Load the BLUPs
yy = read.csv("blues.csv _phenotype.csv")
str(yy)
dim(yy)
dim(Z)
# Load the numerical SNP matrix (the output of "SNP_to_X.r")
# Note it's an RDS (single data object) file hence read using readRDS.
Z = readRDS("zmat_sparse.rds")

#check if the rows of yy and rows of Z are aligned:
all.equal(rownames(Z),yy$Taxa)

# No... seems like the first three lines of BLUPs were not genotyped.
Y = yy

# check again
all.equal(rownames(Z),Y$Taxa)
# We will use SVEN, but there's a tuning parameter called "lam". The default value of 
# lam is n/p^2 where n is the number of lines and p = number of markers.
# But in certain scenarios (e.g. when markers are highly correlated) larger values of 
# lam could be useful (see Li et al., 2020)
# Run ?sven for more info


# First try with lam=0.1... it was useful for Kevin's data on mungbeans
fits = foreach(ii = 2:ncol(Y),.packages = c("bravo","Matrix")) %dopar%
  {
    X = Z + 0
    y = Y[,ii]
    phen.name = names(Y)[ii]
    set.seed(2441139)
    fit = sven(X = X,y = y,Ntemp = 30,Miter = 200,verbose = F,lam = 0.1)
    fit$Phenotype=phen.name
    return(fit)
  }
# save the result
saveRDS(object = fits,file="GWAS-results-1.rds")

rm(fits)


################
library("bravo")

# We need the marker matrix to get the SNP names
Z = readRDS("zmat_sparse.rds")

snp.names = colnames(Z)

gwas = readRDS("GWAS-results-1.rds")
for(i in 1:length(gwas)) {  # Now loops through ALL traits
  cat("\n=== Results for:", gwas[[i]]$Phenotype, "===\n")
  print(data.frame(SNP = snp.names[gwas[[i]]$model.wam],
                   MIP = gwas[[i]]$mip.wam))
}
###results with PVE%
library(bravo)
library(Matrix)


gwas = readRDS("GWAS-results-1.rds")
for(i in 1:length(gwas)) {
  model = gwas[[i]]$model.wam
  trait = gwas[[i]]$stats$y
  Xmat = as.matrix(gwas[[i]]$stats$Xm[, model, drop=FALSE])
  
  pve = as.numeric(cor(trait,Xmat)^2)
  Summary = mip.sven(gwas[[i]],threshold = -1)[model,]
  Summary$PVE = round(pve * 100,1)
  cat("\n Phenotype =  ",gwas[[i]]$Phenotype,"\n")
  print(Summary)
}

