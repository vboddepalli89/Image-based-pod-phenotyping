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
snp = read.table("pod.txt.gz",header=T)
dim(snp)
snp$X.Marker.
X = as.matrix(snp[,-1])
X = 1 - X
print(mean(X != 0))
mx = colMeans(X) > 0.05
Z = Matrix(X[,mx],sparse = T)

rownames(Z) <- snp$X.Marker.
saveRDS(Z,file="zmat_sparse.rds")

##################
library(Matrix)
library(bravo)
library(doParallel) # For parallelizing on Linux/Unix machines

registerDoMC(cores = 1)

# Load the BLUPs
yy = read.csv("blues.csv _phenotype.csv")
str(yy)
dim(yy)
dim(Z)
Z = readRDS("zmat_sparse.rds")
all.equal(rownames(Z),yy$Taxa)
Y = yy

# check again
all.equal(rownames(Z),Y$Taxa)
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


