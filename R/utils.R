dataSummary <- function(dt) {
  cat("\nShape: ", ncol(dt), " columns x ",
    nrow(dt), " rows",
    sep = ""
  )

  non_num_cols <- names(dt)[!unlist(lapply(dt, is.numeric))]
  unique_counts <- data.table(
    "column" = non_num_cols,
    "data_type" = unlist(lapply(dt[, ..non_num_cols], class),
      use.names = FALSE
    ),
    "unique_count" = lapply(dt[, ..non_num_cols], uniqueN),
    "is_NA?" = lapply(dt[, ..non_num_cols], function(x) any(is.na(x)))
  )
  cat("\nNon-numeric columns:\n\n")
  print(unique_counts)
}


getFig <- function(width, heigth) {
  options(repr.plot.width = width, repr.plot.height = heigth)
}


getR2 <- function(true, predicted) {
  # Calculate the residuals (differences between observed and predicted values)
  residuals <- true - predicted

  # the sum of squared residuals
  ss_residuals <- sum(residuals^2)

  # Calculate the total sum of squares
  mean_true <- mean(true)
  ss_total <- sum((true - mean_true)^2)

  # Calculate R-squared
  r2 <- 1 - (ss_residuals / ss_total)

  return(r2)
}


invTransPCA <- function(pred_matrix, n_Comp = nComp, col_means = mu) {
  pred <- pred_matrix %*% t(train_targets_pca$rotation[, 1:n_Comp])
  pred <- scale(pred, center = -col_means, scale = FALSE)
  return(pred)
}
