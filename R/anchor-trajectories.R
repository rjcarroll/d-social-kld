# anchor-trajectories.R
#
# Print full trajectories for candidate anchor firms from chains 1-3.

library(arrow)
library(data.table)

fy_tbl <- as.data.table(readRDS("data/intermediate/index-tables.rds")$fy_tbl)
pq_files <- sort(list.files("data/intermediate/chains-parquet",
                            pattern = "\\.parquet$", full.names = TRUE))

show_trajectory <- function(tick) {
  fy <- fy_tbl[ticker == tick][order(year_idx)]
  if (nrow(fy) == 0) { cat(sprintf("%s: not found\n", tick)); return(invisible()) }
  cols <- paste0("theta.", fy$fy_idx)
  means_mat <- matrix(0, nrow = length(cols), ncol = 3)
  for (i in 1:3) {
    pq <- read_parquet(pq_files[i], col_select = all_of(cols))
    means_mat[, i] <- colMeans(as.matrix(pq))
    rm(pq)
  }
  avg <- rowMeans(means_mat)
  fy[, year := 1990L + year_idx]
  fy[, theta := avg]
  cat(sprintf("\n=== %s (%d years) ===\n", tick, nrow(fy)))
  cat(sprintf("  Overall mean: %.4f  SD: %.4f  Range: [%.4f, %.4f]\n",
      mean(avg), sd(avg), min(avg), max(avg)))
  cat(sprintf("  All negative: %s  All positive: %s\n\n",
      all(avg < 0), all(avg > 0)))
  print(fy[, .(ticker, year, theta)])
}

# Candidate negative anchors
for (t in c("FAST", "CINF", "HP", "TSYS", "CTAS", "HAL")) {
  show_trajectory(t)
}
