rm(list = ls())

library(data.table)
library(qs) # read qs files
library(arrow) # read .parquet files
library(dendextend) # add info in dendrograms
library(pheatmap)
library(rcdk)
library(fingerprint) # fp.sim.matrix()
library(ggplot2)


setwd("~/github/SingleCellPerturbations")

source("R/utils.R")
# Load data ----

filename <- "data/de_train.parquet"
de_train <- read_parquet(filename)
de_train <- setDT(de_train)

head(de_train)
dataSummary(de_train)
length(unique(de_train$sm_name))

mols <- parse.smiles(unique(de_train$SMILES))
fps <- lapply(mols, get.fingerprint, type='circular')
names(fps) <- unique(de_train$sm_name)

cat("Fingerprint for the first molecule:\n")
fps[[1]]


finger_features <- cbind(data.table("sm_name" = names(fps)),
                         fp.to.matrix(fps))
cat("Data.table with 1024 features representing the chemical structure of each compound")
head(finger_features)
dataSummary(finger_features)

# PCA for dimensionality reduction ----
pca_features_res <- prcomp(finger_features[, !"sm_name"],
                           scale. = FALSE,
                           center = TRUE)
finger_features_pca <- finger_features[, .(sm_name)]
finger_features_pca <- cbind(finger_features_pca,
                             pca_features_res$x)
head(finger_features_pca)
dataSummary(finger_features_pca)

getFig(20, 10)

cumulative_var <- cumsum(pca_features_res$sdev^2 / sum(pca_features_res$sdev^2))

plot(cumulative_var,
     xlab = "Principal Component",
     ylab = "Cumulative Proportion of Variance Explained",
     ylim = c(0, 1), type = "b", cex.axis = 1.5, cex.lab = 1.5, xaxt = "n")
breaks <- seq(0, ncol(pca_features_res$x), by = 10)
axis(1, at = breaks, labels = breaks)


abline(v = c(which.max(cumulative_var >= 0.80)), col = "red", lty = 2)
abline(h = 0.80, col = "red", lty = 2)
