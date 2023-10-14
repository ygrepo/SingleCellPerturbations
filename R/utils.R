


dataSummary <- function(dt) {
  
  cat("\nShape: ", ncol(dt), " columns x ",
      nrow(dt), " rows", sep = "")  
  
  non_num_cols <- names(dt)[!unlist(lapply(dt, is.numeric))]
  unique_counts <- data.table(
    "column" = non_num_cols,
    "data_type" = unlist(lapply(dt[, ..non_num_cols], class),
                         use.names = FALSE),
    "unique_count" = lapply(dt[, ..non_num_cols], uniqueN),
    "is_NA?" = lapply(dt[, ..non_num_cols], function(x) any(is.na(x)))
  )
  cat("\nNon-numeric columns:\n\n")
  print(unique_counts)
}


getFig <- function(width, heigth){
  options(repr.plot.width = width, repr.plot.height = heigth)
}

