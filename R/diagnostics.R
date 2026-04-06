# diagnostics.R
#
# Post-run chain diagnostics for the dynamic IRT model.
#
# Strategy: reads from Parquet files (converted from chain CSVs by
# convert-to-parquet.R). Arrow's column-selective reads keep memory low
# and are much faster than the old awk-from-CSV approach.

rm(list = ls())

library(arrow)
library(posterior)

pq_files <- sort(list.files("data/intermediate/chains-parquet",
                            pattern = "\\.parquet$",
                            full.names = TRUE))
stopifnot("No parquet files found — run convert-to-parquet.R first" =
            length(pq_files) > 0)

index_tbls <- readRDS("data/intermediate/index-tables.rds")
fy_tbl     <- index_tbls$fy_tbl

n_chains <- length(pq_files)
message(sprintf("Found %d chain files.", n_chains))

# Get column names from first file
col_names <- ParquetFileReader$create(pq_files[1])$GetSchema()$names
message(sprintf("Total columns: %d", length(col_names)))


# ── helper: read selected columns from all chains into draws_array ───────────
#
# Uses Arrow's column-selective Parquet reads — only the requested columns
# are loaded into memory.

read_cols <- function(pq_files, col_names_wanted) {
  chain_list <- lapply(pq_files, function(f) {
    tbl <- read_parquet(f, col_select = all_of(col_names_wanted))
    as.matrix(tbl)
  })

  n_iter <- nrow(chain_list[[1]])
  arr <- array(
    NA_real_,
    dim = c(n_iter, length(pq_files), length(col_names_wanted)),
    dimnames = list(
      iteration = NULL,
      chain = seq_along(pq_files),
      variable = col_names_wanted
    )
  )
  for (c in seq_along(pq_files)) {
    arr[, c, ] <- chain_list[[c]]
  }
  rm(chain_list)
  gc(verbose = FALSE)
  as_draws_array(arr)
}


# ══════════════════════════════════════════════════════════════════════════════
# 1. Sampler diagnostics
# ══════════════════════════════════════════════════════════════════════════════

message("\n══ 1. Sampler diagnostics ═══════════════════════════════════════════")

sampler_cols <- c("divergent__", "treedepth__", "energy__")

div_counts   <- integer(n_chains)
maxtd_counts <- integer(n_chains)
ebfmi_vals   <- numeric(n_chains)

for (c in seq_along(pq_files)) {
  dt <- read_parquet(pq_files[c], col_select = all_of(sampler_cols))

  div_counts[c]   <- sum(dt$divergent__)
  maxtd_counts[c] <- sum(dt$treedepth__ >= 10L)

  e <- dt$energy__
  ebfmi_vals[c] <- var(diff(e)) / var(e)
  rm(dt, e)
}
n_draws_per_chain <- nrow(
  read_parquet(pq_files[1], col_select = 1L)
)
gc(verbose = FALSE)

message(sprintf("  Draws per chain: %d", n_draws_per_chain))
message(sprintf("  Divergences:     %s  (total: %d)",
                paste(div_counts, collapse = ", "), sum(div_counts)))
message(sprintf("  Max treedepth:   %s  (total: %d)",
                paste(maxtd_counts, collapse = ", "), sum(maxtd_counts)))
message(sprintf("  E-BFMI:          %s",
                paste(sprintf("%.4f", ebfmi_vals), collapse = ", ")))

n_div_total <- sum(div_counts)
if (n_div_total == 0) {
  message("  ✓ No divergences.")
} else {
  n_total <- n_draws_per_chain * n_chains
  message(sprintf("  ⚠ %d divergences out of %d draws (%.3f%%)",
                  n_div_total, n_total, 100 * n_div_total / n_total))
}

if (all(maxtd_counts == 0)) {
  message("  ✓ No max-treedepth hits.")
} else {
  message("  ⚠ Max-treedepth hits detected — consider increasing max_treedepth.")
}

if (all(ebfmi_vals > 0.3)) {
  message("  ✓ E-BFMI all > 0.3 (healthy).")
} else {
  bad <- which(ebfmi_vals <= 0.3)
  message(sprintf("  ⚠ Low E-BFMI on chain(s) %s — possible funnel/poor exploration.",
                  paste(bad, collapse = ", ")))
}


# ══════════════════════════════════════════════════════════════════════════════
# 2. R-hat and ESS: stratified parameter sample
# ══════════════════════════════════════════════════════════════════════════════

message("\n══ 2. Convergence spot checks (R-hat, ESS) ═════════════════════════")

# Helper to run a check on a named set of columns
check_params <- function(label, cols) {
  message(sprintf("\n── %s (%d params) ──", label, length(cols)))
  draws <- read_cols(pq_files, cols)
  summ  <- summarise_draws(draws, "mean", "sd", "rhat", "ess_bulk", "ess_tail")
  rm(draws); gc(verbose = FALSE)

  message(sprintf("  R-hat range:    [%.4f, %.4f]",
                  min(summ$rhat, na.rm = TRUE),
                  max(summ$rhat, na.rm = TRUE)))
  message(sprintf("  Bulk ESS range: [%.0f, %.0f]",
                  min(summ$ess_bulk, na.rm = TRUE),
                  max(summ$ess_bulk, na.rm = TRUE)))
  message(sprintf("  Tail ESS range: [%.0f, %.0f]",
                  min(summ$ess_tail, na.rm = TRUE),
                  max(summ$ess_tail, na.rm = TRUE)))

  bad <- summ$rhat > 1.01
  if (any(bad, na.rm = TRUE)) {
    message(sprintf("  ⚠ %d params with R-hat > 1.01:", sum(bad, na.rm = TRUE)))
    print(summ[which(bad), ], n = Inf)
  } else {
    message("  ✓ All R-hat ≤ 1.01")
  }

  low_bulk <- summ$ess_bulk < 400
  if (any(low_bulk, na.rm = TRUE)) {
    message(sprintf("  ⚠ %d params with bulk ESS < 400", sum(low_bulk, na.rm = TRUE)))
  }

  summ
}


# --- 2a. Hyperparameters and anchors ---

hyper_cols <- c("sigma_global", "sigma_init")
hyper_summ <- check_params("Hyperparameters & anchors", hyper_cols)
message("\n  Full table:")
print(hyper_summ, n = Inf)


# --- 2b. Alpha (item difficulty) — 50 evenly spaced ---

n_alpha   <- nrow(index_tbls$item_tbl)
alpha_idx <- unique(round(seq(1, n_alpha, length.out = 50)))
alpha_cols <- paste0("alpha.", alpha_idx)
alpha_summ <- check_params("Alpha (item difficulty)", alpha_cols)


# --- 2c. Beta (item discrimination) — 50 evenly spaced ---

beta_idx  <- unique(round(seq(1, n_alpha, length.out = 50)))
beta_cols <- paste0("beta.", beta_idx)
beta_summ <- check_params("Beta (item discrimination)", beta_cols)


# --- 2d. Theta (firm-year scores) — 200 stratified, in chunks of 50 ---

message("\n── Theta (firm-year scores) — 200 stratified ──")

n_theta   <- nrow(fy_tbl)
theta_idx <- unique(round(seq(1, n_theta, length.out = 200)))
theta_cols   <- paste0("theta.", theta_idx)
theta_chunks <- split(theta_cols, ceiling(seq_along(theta_cols) / 50))

theta_summ <- do.call(rbind, lapply(seq_along(theta_chunks), function(i) {
  message(sprintf("  Theta chunk %d/%d (%d params)...",
                  i, length(theta_chunks), length(theta_chunks[[i]])))
  draws <- read_cols(pq_files, theta_chunks[[i]])
  s <- summarise_draws(draws, "mean", "sd", "rhat", "ess_bulk", "ess_tail")
  rm(draws); gc(verbose = FALSE)
  s
}))

message(sprintf("\n  Theta overall (%d sampled):", nrow(theta_summ)))
message(sprintf("  R-hat range:    [%.4f, %.4f]",
                min(theta_summ$rhat, na.rm = TRUE),
                max(theta_summ$rhat, na.rm = TRUE)))
message(sprintf("  Bulk ESS range: [%.0f, %.0f]",
                min(theta_summ$ess_bulk, na.rm = TRUE),
                max(theta_summ$ess_bulk, na.rm = TRUE)))
message(sprintf("  Tail ESS range: [%.0f, %.0f]",
                min(theta_summ$ess_tail, na.rm = TRUE),
                max(theta_summ$ess_tail, na.rm = TRUE)))

bad_theta <- theta_summ$rhat > 1.01
if (any(bad_theta, na.rm = TRUE)) {
  message(sprintf("  ⚠ %d / %d sampled theta with R-hat > 1.01",
                  sum(bad_theta, na.rm = TRUE), nrow(theta_summ)))
  print(theta_summ[which(bad_theta), ], n = Inf)
} else {
  message("  ✓ All sampled theta R-hat ≤ 1.01")
}


# --- 2e. sigma_firm — 50 evenly spaced ---

n_firms <- length(unique(fy_tbl$ticker))
sf_idx  <- unique(round(seq(1, n_firms, length.out = 50)))
sf_cols <- paste0("sigma_firm.", sf_idx)
sf_summ <- check_params("sigma_firm", sf_cols)


# ══════════════════════════════════════════════════════════════════════════════
# 3. Anchor firm trajectories
# ══════════════════════════════════════════════════════════════════════════════

message("\n══ 3. Anchor firm posterior trajectories ════════════════════════════")

check_anchor <- function(ticker, label) {
  fy <- dplyr::filter(fy_tbl, ticker == !!ticker) |> dplyr::arrange(year_idx)
  cols <- paste0("theta.", fy$fy_idx)
  draws <- read_cols(pq_files, cols)
  summ <- summarise_draws(draws, "mean", "sd",
                          "rhat", "ess_bulk",
                          ~quantile(.x, 0.05),
                          ~quantile(.x, 0.95))
  summ$year <- 1990L + fy$year_idx
  rm(draws)
  gc(verbose = FALSE)

  message(sprintf("\n  %s (%s):", ticker, label))
  print(
    dplyr::select(summ, variable, year, mean, sd,
                  `5%`, `95%`, rhat, ess_bulk),
    n = Inf
  )

  if (label == "negative anchor") {
    if (all(summ$mean < 0)) {
      message(sprintf("  ✓ All %s posterior means negative.", ticker))
    } else {
      message(sprintf("  ⚠ %s has positive posterior mean(s).", ticker))
    }
  } else {
    if (all(summ$mean > 0)) {
      message(sprintf("  ✓ All %s posterior means positive.", ticker))
    } else {
      message(sprintf("  ⚠ %s has negative posterior mean(s).", ticker))
    }
  }
  summ
}

hal_summ <- check_anchor("HAL", "negative anchor")
abt_summ <- check_anchor("ABT", "positive anchor")


# ══════════════════════════════════════════════════════════════════════════════
# 4. Overall summary
# ══════════════════════════════════════════════════════════════════════════════

message("\n══ 4. Overall summary ═══════════════════════════════════════════════")

all_rhats <- c(hyper_summ$rhat, alpha_summ$rhat, beta_summ$rhat,
               theta_summ$rhat, sf_summ$rhat)
all_bulk  <- c(hyper_summ$ess_bulk, alpha_summ$ess_bulk, beta_summ$ess_bulk,
               theta_summ$ess_bulk, sf_summ$ess_bulk)
all_tail  <- c(hyper_summ$ess_tail, alpha_summ$ess_tail, beta_summ$ess_tail,
               theta_summ$ess_tail, sf_summ$ess_tail)

n_checked  <- length(all_rhats)
n_bad_rhat <- sum(all_rhats > 1.01, na.rm = TRUE)
n_low_bulk <- sum(all_bulk < 400, na.rm = TRUE)
n_low_tail <- sum(all_tail < 400, na.rm = TRUE)

message(sprintf("  Parameters spot-checked: %d", n_checked))
message(sprintf("  R-hat > 1.01:            %d (%.1f%%)",
                n_bad_rhat, 100 * n_bad_rhat / n_checked))
message(sprintf("  Bulk ESS < 400:          %d (%.1f%%)",
                n_low_bulk, 100 * n_low_bulk / n_checked))
message(sprintf("  Tail ESS < 400:          %d (%.1f%%)",
                n_low_tail, 100 * n_low_tail / n_checked))
message(sprintf("  Overall R-hat range:     [%.4f, %.4f]",
                min(all_rhats, na.rm = TRUE), max(all_rhats, na.rm = TRUE)))
message(sprintf("  Overall min bulk ESS:    %.0f", min(all_bulk, na.rm = TRUE)))
message(sprintf("  Overall min tail ESS:    %.0f", min(all_tail, na.rm = TRUE)))

if (n_bad_rhat == 0 && n_div_total == 0 && all(ebfmi_vals > 0.3)) {
  message("\n  ✓ All diagnostics pass. Chains look healthy.")
} else {
  message("\n  ⚠ Some diagnostics flagged — review details above.")
}

message("\ndiagnostics.R complete.")
