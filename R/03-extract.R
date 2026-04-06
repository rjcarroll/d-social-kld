# 03-extract.R
#
# Extracts posterior summaries for all three sets of estimated quantities:
#
#   theta  — firm-year latent CSR scores (the D-SOCIAL-KLD scores)
#   alpha  — item difficulty parameters
#   beta   — item discrimination parameters
#
# Outputs (all committed to the repo):
#
#   data/output/d-social-kld_scores.csv        — firm-year theta summaries
#   data/output/d-social-kld_item-params.csv   — item alpha and beta summaries
#
# Both files report the posterior mean, SD, and a set of percentiles. The
# scores file is the primary output for downstream research use.
#
# ── Memory strategy ──────────────────────────────────────────────────────────
#
# Reads from Parquet files (converted from chain CSVs by convert-to-parquet.R).
# Arrow's column-selective reads keep memory low and are much faster than the
# old awk-from-CSV approach. Theta is still processed in chunks.


# ── setup ─────────────────────────────────────────────────────────────────────

rm(list = ls())

library(arrow)
library(tidyverse)

# ── load ──────────────────────────────────────────────────────────────────────

pq_files <- sort(list.files("data/intermediate/chains-parquet",
                            pattern = "\\.parquet$", full.names = TRUE))
stopifnot("No parquet files found — run convert-to-parquet.R first" =
            length(pq_files) > 0)

index_tbls  <- readRDS("data/intermediate/index-tables.rds")
kld_meta    <- read_csv("data/intermediate/firm-year-level.csv",
                        show_col_types = FALSE)

firm_tbl <- index_tbls$firm_tbl   # ticker → firm_idx
item_tbl <- index_tbls$item_tbl   # item_col × year_idx → item_idx
fy_tbl   <- index_tbls$fy_tbl     # ticker × year_idx → fy_idx

n_chains <- length(pq_files)
message(sprintf("Found %d chain files.", n_chains))


# ── helper: read selected columns from all chains ───────────────────────────
#
# Uses Arrow's column-selective Parquet reads, then row-binds across chains.
# Returns a matrix with (n_draws_total × n_cols) — all chains concatenated.

read_cols_flat <- function(pq_files, col_names_wanted) {
  chain_mats <- lapply(pq_files, function(f) {
    as.matrix(read_parquet(f, col_select = all_of(col_names_wanted)))
  })
  do.call(rbind, chain_mats)
}


# ── helper: posterior percentiles ────────────────────────────────────────────

#' Summarise a draws matrix (iterations × parameters) into a tidy data frame.
#' Operates column-wise to avoid the memory cost of pivot_longer on large draws.
summarise_posterior <- function(draws_mat, var_names) {
  tibble(
    param  = var_names,
    mean   = colMeans(draws_mat),
    sd     = apply(draws_mat, 2, sd),
    min    = apply(draws_mat, 2, min),
    pct05  = apply(draws_mat, 2, quantile, 0.05),
    pct10  = apply(draws_mat, 2, quantile, 0.10),
    pct25  = apply(draws_mat, 2, quantile, 0.25),
    pct50  = apply(draws_mat, 2, quantile, 0.50),
    pct75  = apply(draws_mat, 2, quantile, 0.75),
    pct90  = apply(draws_mat, 2, quantile, 0.90),
    pct95  = apply(draws_mat, 2, quantile, 0.95),
    max    = apply(draws_mat, 2, max)
  )
}


# ── extract theta (firm-year scores) ─────────────────────────────────────────

message("Extracting theta (firm-year scores)...")

# CmdStan CSV uses dot notation: theta.1, theta.2, ...
theta_csv_names <- paste0("theta.", seq_len(nrow(fy_tbl)))

chunk_size <- 2000L
chunks     <- split(theta_csv_names, ceiling(seq_along(theta_csv_names) / chunk_size))

theta_summary <- bind_rows(lapply(seq_along(chunks), function(i) {
  message(sprintf("  theta chunk %d/%d (%d params)...",
                  i, length(chunks), length(chunks[[i]])))
  draws <- read_cols_flat(pq_files, chunks[[i]])
  # Use original bracket names for downstream parsing
  bracket_names <- sub("^theta\\.", "theta[", chunks[[i]]) |>
    paste0("]")
  result <- summarise_posterior(draws, bracket_names)
  rm(draws); gc(verbose = FALSE)
  result
}))

theta_summary <-
  theta_summary |>
  # Parse the firm-year index from the parameter name, e.g. "theta[42]" → 42.
  mutate(fy_idx = as.integer(str_extract(param, "(?<=\\[)\\d+(?=\\])"))) |>
  # Join back to tickers and year indices.
  left_join(fy_tbl, by = "fy_idx") |>
  # Convert year index to calendar year.
  mutate(year = 1990L + year_idx) |>
  # Join company names.
  left_join(
    distinct(kld_meta, name, ticker, year),
    by = c("ticker", "year")
  ) |>
  select(name, ticker, year, mean, sd, min,
         pct05, pct10, pct25, pct50, pct75, pct90, pct95, max) |>
  arrange(ticker, year)

dir.create("data/output", recursive = TRUE, showWarnings = FALSE)
write_csv(theta_summary, "data/output/d-social-kld_scores.csv")
message(sprintf("Wrote %d firm-year score rows.", nrow(theta_summary)))


# ── extract alpha and beta (item parameters) ──────────────────────────────────
#
# alpha[j] is the difficulty of item j: how much baseline CSR a firm needs to
# have a 50% probability of adopting this indicator (holding beta constant).
# Higher alpha → fewer firms adopt the indicator (it's "harder").
#
# beta[j] is the discrimination of item j: how steeply the probability of
# adoption rises with CSR level. High positive beta → strength indicators that
# strongly distinguish responsible firms. Negative beta → concern indicators
# (more-responsible firms are *less* likely to have the concern coded 1).

message("Extracting alpha and beta (item parameters)...")

n_items <- nrow(item_tbl)

# Alpha — small enough to load in one shot
alpha_csv_names <- paste0("alpha.", seq_len(n_items))
alpha_draws     <- read_cols_flat(pq_files, alpha_csv_names)
alpha_bracket   <- paste0("alpha[", seq_len(n_items), "]")
alpha_summary   <- summarise_posterior(alpha_draws, alpha_bracket) |>
  mutate(item_idx = as.integer(str_extract(param, "(?<=\\[)\\d+(?=\\])")))
rm(alpha_draws); gc(verbose = FALSE)

# Beta
beta_csv_names <- paste0("beta.", seq_len(n_items))
beta_draws     <- read_cols_flat(pq_files, beta_csv_names)
beta_bracket   <- paste0("beta[", seq_len(n_items), "]")
beta_summary   <- summarise_posterior(beta_draws, beta_bracket) |>
  mutate(item_idx = as.integer(str_extract(param, "(?<=\\[)\\d+(?=\\])")))
rm(beta_draws); gc(verbose = FALSE)

item_summary <-
  item_tbl |>
  # Join alpha.
  left_join(
    rename_with(alpha_summary, ~ paste0("alpha_", .), .cols = -c(param, item_idx)),
    by = "item_idx"
  ) |>
  # Join beta.
  left_join(
    rename_with(beta_summary, ~ paste0("beta_", .), .cols = -c(param, item_idx)),
    by = "item_idx"
  ) |>
  mutate(year = 1990L + year_idx) |>
  select(item_col, year, item_idx,
         alpha_mean, alpha_sd, alpha_pct05, alpha_pct50, alpha_pct95,
         beta_mean,  beta_sd,  beta_pct05,  beta_pct50,  beta_pct95) |>
  arrange(year, item_col)

write_csv(item_summary, "data/output/d-social-kld_item-params.csv")
message(sprintf("Wrote %d item-parameter rows.", nrow(item_summary)))


# ── quick sanity checks ───────────────────────────────────────────────────────

message("\n── Sanity checks ────────────────────────────────────────────────────")

# 1. Anchors should respect their sign constraints in the posterior means.
hal_check <- filter(theta_summary, ticker == "HAL")
abt_check <- filter(theta_summary, ticker == "ABT")

hal_neg <- sum(hal_check$mean < 0)
message(sprintf(
  "HAL (Halliburton): %d/%d years with negative posterior mean (pinned at -1 in first year only).",
  hal_neg, nrow(hal_check)
))

if (all(abt_check$mean > 0)) {
  message("ABT (Abbott Labs): all posterior means positive. \u2713")
} else {
  warning("ABT has year(s) with negative posterior mean — check identification.")
}

# 2. Distribution of beta signs: most strengths should be positive, most
#    concerns negative. (Not a hard constraint, but a useful check.)
beta_sign_check <-
  item_summary |>
  mutate(
    is_concern  = str_detect(item_col, "_con_"),
    beta_neg    = beta_mean < 0
  ) |>
  group_by(is_concern, beta_neg) |>
  summarise(n = n(), .groups = "drop")

message("\nBeta sign distribution (concern vs strength indicators):")
print(beta_sign_check)

message("\n03-extract.R complete.")
