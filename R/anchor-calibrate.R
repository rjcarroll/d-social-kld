# anchor-calibrate.R
#
# Check HAL and ABT trajectories across chains 1-3 (chain 4 reflected)
# to calibrate fixed anchor values.

library(arrow)
library(data.table)

fy_tbl <- as.data.table(readRDS("data/intermediate/index-tables.rds")$fy_tbl)

hal_fy <- fy_tbl[ticker == "HAL"][order(year_idx)]
abt_fy <- fy_tbl[ticker == "ABT"][order(year_idx)]

pq_files <- sort(list.files("data/intermediate/chains-parquet",
                            pattern = "\\.parquet$", full.names = TRUE))

# ── HAL trajectory ──────────────────────────────────────────────────────────

cat("=== HAL trajectory (chains 1-3 averaged) ===\n\n")
hal_cols <- paste0("theta.", hal_fy$fy_idx)
hal_means <- matrix(0, nrow = length(hal_cols), ncol = 3)
for (i in 1:3) {
  pq <- read_parquet(pq_files[i], col_select = all_of(hal_cols))
  hal_means[, i] <- colMeans(as.matrix(pq))
  rm(pq)
}
hal_avg <- rowMeans(hal_means)
hal_fy[, year := 1990L + year_idx]
hal_fy[, theta := hal_avg]
print(hal_fy[, .(ticker, year, theta)])
cat(sprintf("\nHAL overall mean: %.4f\n", mean(hal_avg)))
cat(sprintf("HAL first-year mean: %.4f\n", hal_avg[1]))
cat(sprintf("HAL range: [%.4f, %.4f]\n", min(hal_avg), max(hal_avg)))


# ── ABT trajectory ──────────────────────────────────────────────────────────

cat("\n=== ABT trajectory (chains 1-3 averaged) ===\n\n")
abt_cols <- paste0("theta.", abt_fy$fy_idx)
abt_means <- matrix(0, nrow = length(abt_cols), ncol = 3)
for (i in 1:3) {
  pq <- read_parquet(pq_files[i], col_select = all_of(abt_cols))
  abt_means[, i] <- colMeans(as.matrix(pq))
  rm(pq)
}
abt_avg <- rowMeans(abt_means)
abt_fy[, year := 1990L + year_idx]
abt_fy[, theta := abt_avg]
print(abt_fy[, .(ticker, year, theta)])
cat(sprintf("\nABT overall mean: %.4f\n", mean(abt_avg)))
cat(sprintf("ABT first-year mean: %.4f\n", abt_avg[1]))
cat(sprintf("ABT range: [%.4f, %.4f]\n", min(abt_avg), max(abt_avg)))


# ── Cross-chain consistency ─────────────────────────────────────────────────

cat("\n=== Per-chain first-year values ===\n\n")
for (i in 1:3) {
  cat(sprintf("Chain %d:  HAL init=%.4f  ABT init=%.4f\n",
      i, hal_means[1, i], abt_means[1, i]))
}

cat("\n=== Per-chain HAL/ABT overall means ===\n\n")
for (i in 1:3) {
  cat(sprintf("Chain %d:  HAL mean=%.4f  ABT mean=%.4f\n",
      i, mean(hal_means[, i]), mean(abt_means[, i])))
}
