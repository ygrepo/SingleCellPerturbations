rm(list = ls())
library(rJava)
library(arrow)
library(data.table)
library(qs)
library(tictoc)
library(ggplot2)
library(glmnet)
library(parallel)
library(doParallel)

setwd("~/github/SingleCellPerturbations")

source("R/utils.R")

# Load PCA of molecular descriptions ----

train_features_source <- fread("data/pca170s0.1_dct_target_encoded_features_train.csv")
test_features_source <- fread("data/pca170s0.1_dct_target_encoded_features_test.csv")

adata_obs_meta <- fread("data/adata_obs_meta.csv")
sum_by_drug_and_cell <- adata_obs_meta[,
  .(n_cells = .N),
  by = c("sm_name", "cell_type")
]
sum_by_drug_and_cell <- sum_by_drug_and_cell[sm_name != "Dimethyl Sulfoxide"]

head(sum_by_drug_and_cell)
cat("\nShape: ", ncol(sum_by_drug_and_cell), " columns x ",
  nrow(sum_by_drug_and_cell), " rows",
  sep = ""
)

# Train Data ----
id_map <- fread("data/id_map.csv")
de_train <- read_parquet("data/de_train.parquet")
de_train <- setDT(de_train)

head(de_train)
cat("\nShape: ", ncol(de_train), " columns x ",
  nrow(de_train), " rows",
  sep = ""
)

# Train Target ----
target_cols <- c(
  "sm_name", # "cell_type",
  names(de_train)[unlist(lapply(de_train, is.numeric))]
)

train_targets <- de_train[, ..target_cols]

# Average targets by drugs
# train_targets <- train_targets[, lapply(.SD, mean),
#                                  .SDcols = names(train_targets)[-1],
#                                  by = "sm_name"]

cat("Train targets:\n")
head(train_targets[1:5, 1:10])
cat("\nShape: ", ncol(train_targets), " columns x ",
  nrow(train_targets), " rows",
  sep = ""
)

# Targets PCA ----

target_cols <- c(names(de_train)[unlist(lapply(de_train, is.numeric))])
mu <- colMeans(train_targets[, ..target_cols])
train_targets_pca <- prcomp(train_targets[, ..target_cols])

nComp <- 10
train_targets <- as.data.table(train_targets_pca$x[, 1:nComp])

head(train_targets)
cat("\nShape: ", ncol(train_targets), " columns x ",
  nrow(train_targets), " rows",
  sep = ""
)

options(repr.plot.width = 20, repr.plot.height = 10)

cumulative_var <- cumsum(train_targets_pca$sdev^2 / sum(train_targets_pca$sdev^2))

plot(cumulative_var,
  xlab = "Principal Component",
  ylab = "Cumulative Proportion of Variance Explained",
  ylim = c(0, 1), type = "b", cex.axis = 1.5, cex.lab = 1.5, xaxt = "n"
)
breaks <- seq(0, ncol(train_targets_pca$x), by = 10)
axis(1, at = breaks, labels = breaks)

# abline(v = c(which.max(cumulative_var >= 0.85)), col = "red", lty = 2)
# abline(h = 0.85, col = "red", lty = 2)
# abline(v = c(which.max(cumulative_var >= 0.90)), col = "red", lty = 2)
# abline(h = 0.90, col = "red", lty = 2)
abline(v = c(which.max(cumulative_var >= 0.95)), col = "red", lty = 2)
abline(h = 0.95, col = "red", lty = 2)

# Train features - target encoded ----


# select n features
n_pc <- 20
select_cols <- c(
  "sm_name", "cell_type",
  names(train_features_source)[grepl(
    "_drug",
    names(train_features_source)
  )][1:n_pc],
  names(train_features_source)[grepl(
    "_cell_type",
    names(train_features_source)
  )][1:n_pc]
)

train_features <- train_features_source[, ..select_cols]
train_features <- sum_by_drug_and_cell[train_features,
  on = c("sm_name", "cell_type")
]

cat("Train features:\n")
head(train_features)
cat("\nShape: ", ncol(train_features), " columns x ",
  nrow(train_features), " rows",
  sep = ""
)

# Test features - target encoded ----

select_cols <- c(
  "sm_name", "cell_type",
  names(test_features_source)[grepl(
    "_drug",
    names(test_features_source)
  )][1:n_pc],
  names(test_features_source)[grepl(
    "_cell_type",
    names(test_features_source)
  )][1:n_pc]
)

test_features <- test_features_source[, ..select_cols]
mean_test <- sum_by_drug_and_cell[cell_type %in% test_features$cell_type,
  .(n_cells = mean(n_cells)),
  by = "cell_type"
]
test_features <- mean_test[test_features, on = "cell_type"]
setcolorder(test_features, names(train_features))

cat("Test features:\n")
head(test_features)
cat("\nShape: ", ncol(test_features), " columns x ",
  nrow(test_features), " rows",
  sep = ""
)

# Remove ID columns ----

id_cols <- c("sm_name", "cell_type")
train_features <- train_features[, !..id_cols]
test_features <- test_features[, !..id_cols]

dim(train_targets)
dim(train_features)
dim(test_features)

# Ridge ----

n_folds <- 7
set.seed(1)
cv <- data.table(id = 1:train_features[, .N])
cv[, folds := sample(1:n_folds, .N, replace = TRUE)]
cv[order(folds), .N, by = folds]
head(cv)

cl <- makePSOCKcluster(detectCores())
registerDoParallel(cl)


submission <- fread("data/sample_submission.csv")
submission[1:5, 1:5]

preds_cv <- melt(submission,
  id.vars = "id",
  variable.name = "gene",
  value.name = "pred"
)
Y_oof <- NULL

for (fold in 1:n_folds) {
  tic()
  cat("\n", "fold:", fold, "\n")

  # train/test split for local validation
  train_ids <- cv[folds != fold, id]
  valid_ids <- cv[folds == fold, id]

  X_train <- as.matrix(train_features[train_ids])
  Y_train <- as.matrix(train_targets[train_ids])

  X_valid <- as.matrix(train_features[valid_ids])
  Y_valid <- as.matrix(train_targets[valid_ids])

  # Normalisation
  X_train <- scale(X_train)
  
  # Cross-validated lambda tuning
  cv_fit <- cv.glmnet(X_train, Y_train,
    family = "mgaussian", alpha = 1,
    parallel = TRUE
  )

  cat("Lambda min:", cv_fit$lambda.min, "\n")
  cat("Lambda 1se:", cv_fit$lambda.1se, "\n")
  cat("Minimum CV MSE:", min(cv_fit$cvm), "\n")

  # Predict train seen data
  Y_pred <- predict(cv_fit, newx = as.matrix(X_train), s = "lambda.min")
  Y_pred <- apply(Y_pred, 2, c)

  cat(
    "r2 train:",
    getR2(
      invTransPCA(Y_train),
      invTransPCA(Y_pred)
    ), "\n"
  )

  # Predict valid unseen data
  X_valid <- scale(X_valid,
                   center = attr(X_train, "scaled:center"),
                   scale = attr(X_train, "scaled:scale")
  )
  Y_pred <- predict(cv_fit, 
                    newx = as.matrix(X_valid), s = "lambda.min")
  Y_pred <- apply(Y_pred, 2, c)

  cat(
    "r2 valid:",
    getR2(
      invTransPCA(Y_valid),
      invTransPCA(Y_pred)
    ), "\n"
  )

  # Save the predictions for validation set of each fold
  Y_oof <- rbind(Y_oof, data.table(
    id = valid_ids,
    invTransPCA(Y_pred)
  ))
  # Predict test targets and average between folds
  X_submit <- as.matrix(test_features)
  X_submit <- scale(X_submit,
    center = attr(X_train, "scaled:center"),
    scale = attr(X_train, "scaled:scale")
  )
  preds <- predict(cv_fit, newx = X_submit, s = "lambda.min")
  preds <- apply(preds, 2, c)
  preds <- as.data.table(invTransPCA(preds))
  preds[, id := 0:(nrow(preds) - 1)]
  preds <- melt(preds, id.vars = "id")
  preds_cv[, pred := pred + preds[, value] / n_folds]
  toc()
}
stopCluster(cl)

# Evaluation on OOF train data ----

target_cols <- c(names(de_train)[unlist(lapply(de_train, is.numeric))])
true_targets <- as.matrix(de_train[, ..target_cols])

Y_oof2 <- as.matrix(Y_oof[order(id), !"id"])

rmse <- sqrt(mean((true_targets - Y_oof2)^2))
cat("\n", "RMSE: ", rmse, "\n", sep = "")

mrRMSE <- mean(sqrt(rowMeans((true_targets - Y_oof2)^2)))
cat("mrRMSE:", mrRMSE, "\n")

MAE <- mean(abs((true_targets - Y_oof2)))
cat("MAE:", MAE, "\n")

r_squared <- 1 - (sum((true_targets - Y_oof2)^2) / 
                    sum((true_targets - mean(Y_oof2))^2))
cat("R-squared (R2):", r_squared, "\n")

# predict for submit ---

submit <- dcast(preds_cv, 
                id ~ gene,
                value.var = "pred")
head(submit)
cat("\nShape: ", ncol(submit), " columns x ",
    nrow(submit), " rows", sep = "")

fwrite(submit, "submission.csv")
