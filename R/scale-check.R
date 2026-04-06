# scale-check.R
#
# Per-chain diagnostics to understand the scale non-identification.
# Reads from Parquet files (fast column-selective access).

library(arrow)
library(data.table)

fy_tbl   <- as.data.table(readRDS("data/intermediate/index-tables.rds")$fy_tbl)
item_tbl <- as.data.table(readRDS("data/intermediate/index-tables.rds")$item_tbl)

pq_files <- sort(list.files("data/intermediate/chains-parquet",
                            pattern = "\\.parquet$", full.names = TRUE))

# ── 1. Per-chain anchor values ───────────────────────────────────────────────

cat("=== Per-chain anchor values ===\n\n")
for (i in seq_along(pq_files)) {
  pq <- read_parquet(pq_files[i],
                     col_select = c("theta_hal_init", "theta_abt_init"))
  cat(sprintf("Chain %d:\n", i))
  cat(sprintf("  HAL init: mean=%.4f  sd=%.4f  [%.4f, %.4f]\n",
      mean(pq$theta_hal_init), sd(pq$theta_hal_init),
      quantile(pq$theta_hal_init, 0.05), quantile(pq$theta_hal_init, 0.95)))
  cat(sprintf("  ABT init: mean=%.4f  sd=%.4f  [%.4f, %.4f]\n",
      mean(pq$theta_abt_init), sd(pq$theta_abt_init),
      quantile(pq$theta_abt_init, 0.05), quantile(pq$theta_abt_init, 0.95)))
  rm(pq)
}


# ── 2. Per-chain sigma_init and sigma_global ─────────────────────────────────

cat("\n=== Per-chain sigma_init and sigma_global ===\n\n")
for (i in seq_along(pq_files)) {
  pq <- read_parquet(pq_files[i],
                     col_select = c("sigma_init", "sigma_global"))
  cat(sprintf("Chain %d:  sigma_init=%.4f (sd=%.4f)  sigma_global=%.5f (sd=%.5f)\n",
      i, mean(pq$sigma_init), sd(pq$sigma_init),
      mean(pq$sigma_global), sd(pq$sigma_global)))
  rm(pq)
}


# ── 3. Per-chain beta distribution ──────────────────────────────────────────

cat("\n=== Per-chain beta distribution (all items, posterior means) ===\n\n")
beta_cols <- paste0("beta.", seq_len(nrow(item_tbl)))
for (i in seq_along(pq_files)) {
  pq <- read_parquet(pq_files[i], col_select = beta_cols)
  beta_means <- colMeans(as.matrix(pq))
  cat(sprintf(paste0("Chain %d:  mean(beta)=%.3f  sd(beta)=%.3f  ",
                     "range=[%.3f, %.3f]  |beta|>5: %d\n"),
      i, mean(beta_means), sd(beta_means),
      min(beta_means), max(beta_means),
      sum(abs(beta_means) > 5)))
  rm(pq, beta_means)
  gc(verbose = FALSE)
}


# ── 4. Per-chain theta SD ───────────────────────────────────────────────────

cat("\n=== Per-chain theta SD (posterior means across all firm-years) ===\n\n")
theta_cols <- paste0("theta.", seq_len(nrow(fy_tbl)))
chunk_size <- 10000L
chunks <- split(theta_cols, ceiling(seq_along(theta_cols) / chunk_size))

for (i in seq_along(pq_files)) {
  all_means <- numeric(0)
  for (ch in chunks) {
    pq <- read_parquet(pq_files[i], col_select = ch)
    all_means <- c(all_means, colMeans(as.matrix(pq)))
    rm(pq)
  }
  gc(verbose = FALSE)
  cat(sprintf("Chain %d:  mean(theta)=%.4f  sd(theta)=%.4f  range=[%.3f, %.3f]\n",
      i, mean(all_means), sd(all_means),
      min(all_means), max(all_means)))
}


# ── 5. Scale ratio between chains ───────────────────────────────────────────
# If the issue is pure scale, the ratio of theta SDs across chains should
# explain the beta ratio too.

cat("\n=== Cross-chain scale ratios ===\n")
cat("(If pure scale issue, theta_sd ratio ≈ 1/beta_sd ratio)\n\n")

theta_sds <- numeric(length(pq_files))
beta_sds  <- numeric(length(pq_files))

for (i in seq_along(pq_files)) {
  # theta SD
  all_means <- numeric(0)
  for (ch in chunks) {
    pq <- read_parquet(pq_files[i], col_select = ch)
    all_means <- c(all_means, colMeans(as.matrix(pq)))
    rm(pq)
  }
  theta_sds[i] <- sd(all_means)

  # beta SD
  pq <- read_parquet(pq_files[i], col_select = beta_cols)
  beta_sds[i] <- sd(colMeans(as.matrix(pq)))
  rm(pq)
  gc(verbose = FALSE)
}

cat(sprintf("  Theta SDs: %s\n", paste(sprintf("%.4f", theta_sds), collapse = ", ")))
cat(sprintf("  Beta SDs:  %s\n", paste(sprintf("%.3f", beta_sds), collapse = ", ")))
cat(sprintf("  Theta ratios (vs chain 1): %s\n",
    paste(sprintf("%.3f", theta_sds / theta_sds[1]), collapse = ", ")))
cat(sprintf("  Beta ratios  (vs chain 1): %s\n",
    paste(sprintf("%.3f", beta_sds / beta_sds[1]), collapse = ", ")))
cat(sprintf("  Theta × Beta (should be constant if pure scale): %s\n",
    paste(sprintf("%.4f", theta_sds * beta_sds), collapse = ", ")))
