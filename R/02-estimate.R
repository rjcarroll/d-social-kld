# 02-estimate.R
#
# Prepares the Stan data list from the reshaped KLD matrix and runs the dynamic
# IRT model. Outputs:
#
#   data/intermediate/stan-data.rds          — the list passed to Stan
#   data/intermediate/stan-fit.rds           — the CmdStanR fit object
#
# ── Before running ────────────────────────────────────────────────────────────
#
# This script is computationally intensive. On the full 1991–2018 dataset
# expect a runtime of many hours (exact time depends on hardware and chain
# configuration). A `demo_mode` flag is provided to run a fast subset suitable
# for testing the pipeline end-to-end.
#
# You will need cmdstanr and a local Stan installation:
#   install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/",
#                                          getOption("repos")))
#   cmdstanr::install_cmdstan()


# ── setup ─────────────────────────────────────────────────────────────────────

rm(list = ls())

library(tidyverse)
library(cmdstanr)

source("R/utils.R")

# ── options ───────────────────────────────────────────────────────────────────

# Set demo_mode = TRUE for a fast run with reduced data and iterations.
# Useful for verifying the pipeline and getting a sense of sampler speed
# before committing to a full overnight run.
#
# demo_years controls which calendar years to include. The smaller the range,
# the faster the run — the early KLD years (1991–2000) cover only ~650 firms,
# making them a good smoke-test. Set to NULL to use all available years.
demo_mode  <- FALSE
demo_years <- 1991:2000

# MCMC settings. Defaults match the original Carroll et al. (2016) paper:
# 5 000 iterations with 1 000 burn-in and thinning by 2 → 2 000 retained draws
# per chain. With multiple chains you get more total draws; the paper used one
# chain (MCMCpack does not easily parallelize). Here we default to 4 chains.
n_chains  <- 4L
n_warmup  <- 1000L
n_iter    <- 2500L   # per chain; 2 500 × 4 chains ≈ 10 000 draws total
n_thin    <- 1L      # Stan's HMC is already efficient; thinning rarely helps

# Random seed for reproducibility.
mcmc_seed <- 61802L


# ── load data ─────────────────────────────────────────────────────────────────

load("data/intermediate/data-reshaped.RData")
# Objects now in environment: dat (tibble), yrs (int vec), all_tickers (chr vec)

# ── optional demo subset ──────────────────────────────────────────────────────

if (demo_mode) {

  yr_min <- min(1990L + yrs)   # earliest calendar year in the reshaped data
  yr_max <- max(1990L + yrs)

  if (is.null(demo_years)) {
    demo_years <- seq(yr_min, yr_max)
  }

  demo_years <- as.integer(demo_years)
  demo_idx   <- demo_years - 1990L   # convert to year-index space

  message(sprintf(
    "Running in demo mode: %d–%d (%d years), reduced iterations.",
    min(demo_years), max(demo_years), length(demo_years)
  ))

  # Keep only columns whose year index falls in the requested range.
  keep_cols <- which(yrs %in% demo_idx)
  dat <- dat[, keep_cols]
  yrs <- yrs[keep_cols]

  # Reduce iterations — enough to see sampler behaviour, not for inference.
  n_warmup <- 200L
  n_iter   <- 300L

}

# ── convert to matrix ─────────────────────────────────────────────────────────

wide_mat           <- as.matrix(dat)
rownames(wide_mat) <- all_tickers

# ── build long-format observation table ───────────────────────────────────────
#
# The Stan model works on a flat list of non-missing (firm-year, item, y)
# triples rather than the full sparse matrix. This avoids passing millions of
# NA placeholders to Stan.

message("Building long observation table...")
obs <- wide_to_long(wide_mat, yrs)
# obs columns: ticker, item_col, year_idx, y

# ── build index tables ────────────────────────────────────────────────────────

# Firm index (alphabetical, matching all_tickers).
firm_tbl <- tibble(ticker = all_tickers, firm_idx = seq_along(all_tickers))

# Item index: each unique (indicator column name, year index) pair.
item_tbl <-
  distinct(obs, item_col, year_idx) |>
  arrange(year_idx, item_col) |>
  mutate(item_idx = row_number())

# Firm-year index: each unique (ticker, year_idx) pair that has at least one
# non-missing observation.
fy_tbl <-
  distinct(obs, ticker, year_idx) |>
  arrange(ticker, year_idx) |>
  mutate(fy_idx = row_number()) |>
  left_join(firm_tbl, by = "ticker")

# Join indices onto the observation table.
obs <-
  obs |>
  left_join(fy_tbl,   by = c("ticker", "year_idx")) |>
  left_join(item_tbl, by = c("item_col", "year_idx"))

# ── identify anchor firms ─────────────────────────────────────────────────────
#
# HAL (Halliburton) and ABT (Abbott Laboratories) serve as identification anchors.
# HAL is constrained to be negative (low-CSR anchor) and ABT positive
# (high-CSR anchor). Both firms have non-missing KLD indicator data for all 28
# years (1991–2018), making them robust anchors across the full panel.
#
# Anchor selection rationale:
#   HAL — Halliburton is consistently associated with environmental and
#         governance controversies, making it a natural low-CSR reference point.
#         It was also used as the negative anchor in Carroll, Primo & Richter
#         (2016), giving continuity with the original measure.
#   ABT — Abbott Laboratories (spun off AbbVie in 2013; the Abbott entity
#         retained the ABT ticker) ranks among the top corporate citizens in
#         virtually every major CSR index (DJSI, MSCI, Fortune). Unlike MSFT or
#         JNJ, ABT has full KLD indicator coverage back to 1991.
#
# If either anchor is absent from the data (after the demo subsetting above),
# the script stops with an informative error.

anchor_neg <- "HAL"
anchor_pos <- "ABT"

stopifnot(
  "HAL (Halliburton) not found in data — cannot identify model direction" =
    anchor_neg %in% fy_tbl$ticker,
  "ABT (Abbott Labs) not found in data — cannot identify model direction" =
    anchor_pos %in% fy_tbl$ticker
)

# ── dynamic model arrays ──────────────────────────────────────────────────────

message("Building dynamic model arrays...")

# "Fresh start" firm-years: first appearances AND years immediately following
# a gap in a firm's panel history. In the NCP model, theta is built
# deterministically in transformed parameters, so every fy_idx must be covered
# by either an init assignment or a transition. Gap years are not connected by
# transitions (build_transitions only links consecutive years), so they must be
# treated as additional fresh starts.
fresh_fy_tbl <-
  fy_tbl |>
  arrange(ticker, year_idx) |>
  group_by(ticker) |>
  mutate(prev_year = lag(year_idx),
         is_fresh  = is.na(prev_year) | year_idx > prev_year + 1L) |>
  filter(is_fresh) |>
  ungroup()

# The identification constraints apply only to the anchors' first appearances.
# Anchor gap years (if any) go into other_init_fy and are treated as free.
hal_first_year <- filter(fresh_fy_tbl, ticker == anchor_neg) |>
  pull(year_idx) |> min()
abt_first_year <- filter(fresh_fy_tbl, ticker == anchor_pos) |>
  pull(year_idx) |> min()

hal_init_fy_idx   <- filter(fresh_fy_tbl,
                              ticker == anchor_neg, year_idx == hal_first_year) |>
  pull(fy_idx)
abt_init_fy_idx   <- filter(fresh_fy_tbl,
                              ticker == anchor_pos, year_idx == abt_first_year) |>
  pull(fy_idx)
other_init_fy_idx <- filter(fresh_fy_tbl,
                              !(ticker == anchor_neg & year_idx == hal_first_year),
                              !(ticker == anchor_pos & year_idx == abt_first_year)) |>
  pull(fy_idx)

# Consecutive transitions.
trans_tbl <- build_transitions(fy_tbl)

# ── compile the Stan model ────────────────────────────────────────────────────

message(sprintf("[%s] Compiling Stan model...", Sys.time()))
model <- cmdstan_model(stan_file = "stan/dynamic_irt.stan")

# ── assemble the Stan data list ───────────────────────────────────────────────

stan_data <- list(

  N        = nrow(obs),
  N_fy     = nrow(fy_tbl),
  N_items  = nrow(item_tbl),
  N_firms  = nrow(firm_tbl),

  y      = obs$y,
  fy_obs = obs$fy_idx,
  it_obs = obs$item_idx,

  fy_firm = fy_tbl$firm_idx,

  hal_init_fy   = hal_init_fy_idx,
  abt_init_fy   = abt_init_fy_idx,
  N_init_other  = length(other_init_fy_idx),
  other_init_fy = other_init_fy_idx,

  N_trans    = nrow(trans_tbl),
  trans_prev = trans_tbl$trans_prev,
  trans_curr = trans_tbl$trans_curr,
  trans_firm = trans_tbl$trans_firm

)

message(sprintf(
  "Stan data: %d observations | %d firm-years | %d items | %d firms",
  stan_data$N, stan_data$N_fy, stan_data$N_items, stan_data$N_firms
))

saveRDS(stan_data, file = "data/intermediate/stan-data.rds")

# Also save the index tables — 03-extract.R needs them to map parameter indices
# back to tickers, years, and indicator names.
saveRDS(
  list(firm_tbl = firm_tbl, item_tbl = item_tbl, fy_tbl = fy_tbl),
  file = "data/intermediate/index-tables.rds"
)

# ── run the sampler ───────────────────────────────────────────────────────────
#
# Note on `adapt_delta`: the default (0.8) is usually fine for well-specified
# models. If you see divergences in the diagnostics, try 0.9 or 0.95.
#
# Note on `max_treedepth`: the default is 10. If you see many "max treedepth"
# warnings, increase to 12. This adds compute time but improves mixing.

message(sprintf("[%s] Starting MCMC... (this may take many hours on the full dataset)", Sys.time()))

t_mcmc_start <- proc.time()

dir.create("data/intermediate/chains", recursive = TRUE, showWarnings = FALSE)
fit <- model$sample(
  data            = stan_data,
  seed            = mcmc_seed,
  chains          = n_chains,
  parallel_chains = min(n_chains, parallel::detectCores()),
  iter_warmup     = n_warmup,
  iter_sampling   = n_iter,
  thin            = n_thin,
  adapt_delta     = 0.9,
  max_treedepth   = 10,
  init            = 0,
  refresh         = 10,
  output_dir      = "data/intermediate/chains"
)

t_mcmc_elapsed <- (proc.time() - t_mcmc_start)[["elapsed"]]

# ── timing summary ──────────────────────────────────────────────────────────
#
# Per-chain wall-clock times from Stan, plus projections for larger runs.

chain_times <- fit$time()$chains
total_iter  <- n_warmup + n_iter

message("\n── Timing ──────────────────────────────────────────────────")
message(sprintf("  Wall clock (R):       %.1f s", t_mcmc_elapsed))
message(sprintf("  Observations:         %s", format(stan_data$N, big.mark = ",")))
message(sprintf("  Firm-years:           %s", format(stan_data$N_fy, big.mark = ",")))
message(sprintf("  Items:                %d", stan_data$N_items))
message(sprintf("  Chains:               %d", n_chains))
message(sprintf("  Iterations per chain: %d (%d warmup + %d sampling)",
                total_iter, n_warmup, n_iter))

# Per-chain breakdown.
for (i in seq_len(nrow(chain_times))) {
  message(sprintf(
    "  Chain %d: %.1f s warmup, %.1f s sampling, %.1f s total",
    i, chain_times$warmup[i], chain_times$sampling[i],
    chain_times$total[i]
  ))
}

# Cost per iteration (using the slowest chain, since chains run in parallel).
slowest    <- max(chain_times$total)
sec_per_it <- slowest / total_iter
message(sprintf(
  "  Slowest chain: %.2f s/iteration (%.1f ms per obs per iteration)",
  sec_per_it, sec_per_it / stan_data$N * 1000
))

# Projection for full run (1991-2018).
full_obs <- 1.5e6   # approximate non-missing obs in full panel
full_warmup  <- 1000L
full_sampling <- 1500L
full_iter <- full_warmup + full_sampling
proj_sec  <- sec_per_it * (full_obs / stan_data$N) * full_iter
message(sprintf(
  "\n  ── Projection for full run (%.0fM obs, %d iterations) ──",
  full_obs / 1e6, full_iter
))
message(sprintf("  Estimated: %.1f hours (%.1f days)",
                proj_sec / 3600, proj_sec / 86400))
message("  (Rough: assumes linear scaling in obs and iterations.)")

# ── save fit ──────────────────────────────────────────────────────────────────
#
# Save before diagnostics so the fit is never lost to a downstream memory error.

message(sprintf("[%s] Saving fit object...", Sys.time()))
fit$save_object(file = "data/intermediate/stan-fit.rds")

# ── basic diagnostics ─────────────────────────────────────────────────────────
#
# Always check these before trusting the output:
#   - No divergences (or very few: < 0.1% of post-warmup draws)
#   - R-hat < 1.01 for all parameters
#   - Bulk and tail ESS > 100 per chain (ideally > 400 total)

message(sprintf(
  "\n[%s] ── Sampler diagnostics ──────────────────────────────",
  Sys.time()
))
fit$diagnostic_summary()

diag_tbl <- fit$summary(
  variables = NULL,
  "mean", "sd", "rhat", "ess_bulk", "ess_tail"
)
message(sprintf(
  "R-hat range: [%.4f, %.4f]",
  min(diag_tbl$rhat, na.rm = TRUE),
  max(diag_tbl$rhat, na.rm = TRUE)
))
message(sprintf(
  "Min bulk ESS: %.0f | Min tail ESS: %.0f",
  min(diag_tbl$ess_bulk, na.rm = TRUE),
  min(diag_tbl$ess_tail, na.rm = TRUE)
))

message(sprintf("[%s] 02-estimate.R complete.", Sys.time()))
