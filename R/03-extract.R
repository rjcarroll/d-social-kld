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


# ── setup ─────────────────────────────────────────────────────────────────────

rm(list = ls())

library(tidyverse)
library(cmdstanr)
library(posterior)   # for as_draws_df(), summarise_draws()

# ── load ──────────────────────────────────────────────────────────────────────

# Reconstruct the fit from chain CSVs (avoids loading a monolithic RDS that
# exceeds available RAM). Only the CSV paths are needed; draws are read lazily
# when we call fit$draws() below.
csv_files <- list.files("data/intermediate/chains", pattern = "\\.csv$",
                        full.names = TRUE)
stopifnot("No chain CSVs found in data/intermediate/chains/" = length(csv_files) > 0)
fit <- as_cmdstan_fit(csv_files)

index_tbls  <- readRDS("data/intermediate/index-tables.rds")
kld_meta    <- read_csv("data/intermediate/firm-year-level.csv",
                        show_col_types = FALSE)

firm_tbl <- index_tbls$firm_tbl   # ticker → firm_idx
item_tbl <- index_tbls$item_tbl   # item_col × year_idx → item_idx
fy_tbl   <- index_tbls$fy_tbl     # ticker × year_idx → fy_idx


# ── helper: posterior percentiles ────────────────────────────────────────────

#' Summarise a draws matrix (iterations × parameters) into a tidy data frame.
#' Operates column-wise to avoid the memory cost of pivot_longer on large draws.
summarise_posterior <- function(draws_df, vars) {
  draws_mat <- as.matrix(draws_df[, vars, drop = FALSE])
  tibble(
    param  = vars,
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

# theta has ~51K parameters × 10K draws — too large to load at once on 16 GB.
# Process in chunks, summarise each, then combine.
all_theta_vars <- fit$metadata()$stan_variables
all_theta_vars <- grep("^theta$", all_theta_vars, value = TRUE)

# Get the full list of theta[i] column names from one chain's metadata.
theta_names <- paste0("theta[", seq_len(nrow(fy_tbl)), "]")

chunk_size <- 5000L
chunks     <- split(theta_names, ceiling(seq_along(theta_names) / chunk_size))

theta_summary <- bind_rows(lapply(seq_along(chunks), function(i) {
  message(sprintf("  theta chunk %d/%d (%d params)...",
                  i, length(chunks), length(chunks[[i]])))
  draws <- as_draws_df(fit$draws(variables = chunks[[i]]))
  result <- summarise_posterior(draws, chunks[[i]])
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

# alpha and beta are much smaller (~2K each) — safe to load in one shot.
alpha_draws <- as_draws_df(fit$draws("alpha"))
alpha_vars  <- grep("^alpha\\[", names(alpha_draws), value = TRUE)
alpha_summary <-
  summarise_posterior(alpha_draws, alpha_vars) |>
  mutate(item_idx = as.integer(str_extract(param, "(?<=\\[)\\d+(?=\\])")))
rm(alpha_draws); gc(verbose = FALSE)

beta_draws <- as_draws_df(fit$draws("beta"))
beta_vars  <- grep("^beta\\[",  names(beta_draws),  value = TRUE)
beta_summary <-
  summarise_posterior(beta_draws, beta_vars) |>
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

if (all(hal_check$mean < 0)) {
  message("HAL (Halliburton): all posterior means negative. \u2713")
} else {
  warning("HAL has year(s) with positive posterior mean — check identification.")
}

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
