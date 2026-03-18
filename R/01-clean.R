# 01-clean.R
#
# Reads the raw KLD STATS data, resolves data-quality issues, and reshapes it
# into the wide firm × item matrix that the Stan model expects.
#
# Ticker disambiguation strategy
# ────────────────────────────────
# Tickers are unstable: companies change them, and distinct KLD "domiciles"
# (e.g. US vs. international coverage) reuse the same symbol for different
# firms. The strategy here is:
#
#   • Where `issuerid` is populated (2017–2018 in this dataset), we use it to
#     identify collisions between truly distinct firms sharing a ticker.
#
#   • For all years, when a ticker is claimed by more than one firm name, the
#     firm that *first* appeared under that ticker keeps the original symbol;
#     later entrants receive a deterministic numeric suffix: TICK.2, TICK.3, …
#     A log carries these suffixes forward across years.
#
# This replaces the earlier random-tag approach, which was non-reproducible
# without a fixed seed and susceptible to hash collisions.
#
# Outputs
# ────────
#   data/intermediate/firm-year-level.csv   — long-format firm-year metadata
#   data/intermediate/data-reshaped.RData   — wide matrix + year index for Stan


# ── setup ─────────────────────────────────────────────────────────────────────

rm(list = ls())

library(tidyverse)
library(foreach)

source("R/utils.R")

# ── load ──────────────────────────────────────────────────────────────────────

kld <-
  read_csv(
    file      = "data/kld/KLD_stats_all_through_2018.csv",
    na        = c("", "R", "NA"),
    col_types = cols(
      issuerid         = col_character(),   # kept: used for disambiguation
      CUSIP            = col_skip(),
      legacy_companyID = col_skip(),
      domicile         = col_skip()
    )
  ) |>
  # HUM_str_num has a parsing issue when its value is 2; clamp to 1.
  mutate(HUM_str_num = if_else(HUM_str_num >= 1, 1, 0)) |>
  # Drop cumulative count columns — we only want the binary indicators.
  select(-contains("num")) |>
  select(CompanyName, Ticker, issuerid, year, everything()) |>
  rename(name = CompanyName, ticker = Ticker)

# ── standardize to 0/1/NA integers ────────────────────────────────────────────

kld <-
  kld |>
  mutate(across(where(is.logical), logical_to_int)) |>
  mutate(across(where(is.double) & !year, as.integer)) |>
  mutate(across(where(is.integer) & !year, clamp_binary))

# ── restrict to firms present 1991–2012 ───────────────────────────────────────
#
# Keeps the extended (through 2018) scores comparable to the published measure.
# Firms that enter the KLD database only after 2012 are excluded.

universe <-
  kld |>
  filter(year <= 2012) |>
  distinct(ticker) |>
  pull(ticker)

kld <- filter(kld, ticker %in% universe)
rm(universe)

# ── fix a small set of known upfront ticker problems ──────────────────────────
#
# These are cases where the raw data has an outright error (wrong ticker) rather
# than a collision between two legitimately different firms. Everything else is
# handled by the systematic disambiguation loop below.

# "X" (U.S. Steel NYSE ticker) is read as NA by readr; restore it.
kld <- mutate(kld, ticker = if_else(ticker == "X", "USX", ticker))

# Firms whose ticker was missing due to the NA read-in rule above.
kld$ticker[kld$name %in% c("Ryder System, Inc.", "RYDER SYSTEM, INC.")] <- "R"

# HomeAway was incorrectly coded under Expedia's ticker.
kld$ticker[kld$name == "HOMEAWAY, INC."] <- "AWAY"

# ── identify which metrics appear in each year ─────────────────────────────────
#
# Indicators are added to KLD over time; a column that is all-NA in a given
# year is effectively absent and should be excluded from that year's items.

has_data <- function(x) any(!is.na(x))

mets_by_year <-
  foreach(yr = seq(min(kld$year), max(kld$year))) %do% {
    kld |>
      filter(year == yr) |>
      select(-name, -ticker, -issuerid, -year) |>
      select(where(has_data)) |>
      names()
  }

# ── build and clean the per-year firm list ─────────────────────────────────────

yr_min    <- min(kld$year)
yr_max    <- max(kld$year)
yr_offset <- yr_min - 1L
n_years   <- yr_max - yr_offset

firms_in_year <-
  foreach(yr = seq(yr_min, yr_max)) %do% {
    kld |> filter(year == yr) |> pull(ticker) |> sort()
  }

# ── ticker disambiguation ─────────────────────────────────────────────────────

ticker_log <- dplyr::tibble(
  name           = character(),
  old_ticker     = character(),
  new_ticker     = character(),
  year_created   = integer(),
  year_last_used = integer()
)

for (k in seq_len(n_years)) {

  yr <- yr_offset + k

  # Propagate previously assigned suffixed tickers into this year.
  if (nrow(ticker_log) > 0) {
    out        <- propagate_ticker_log(kld, ticker_log, yr)
    kld        <- out$kld
    ticker_log <- out$log
  }

  # Where issuerid is available, use it to identify true firm-level collisions.
  out        <- resolve_with_issuerid(kld, yr, ticker_log)
  kld        <- out$kld
  ticker_log <- out$log

  # Resolve remaining name-based collisions with deterministic suffixes.
  out        <- resolve_duplicate_tickers(kld, yr, ticker_log)
  kld        <- out$kld
  ticker_log <- out$log

  # Rebuild the firm list for this year from the (now updated) kld.
  firms_in_year[[k]] <-
    kld |>
    filter(year == yr) |>
    distinct(ticker) |>   # distinct() handles any exact-duplicate rows
    pull(ticker) |>
    sort()
}

# Verify: all tickers are now distinct within every year.
bad_years <- which(!vapply(firms_in_year, all_distinct, logical(1L)))
if (length(bad_years) > 0) {
  stop("Ticker disambiguation incomplete for year(s): ",
       paste(bad_years + yr_offset, collapse = ", "),
       "\nInspect kld for those years and extend the resolution logic.")
}
message("Ticker disambiguation complete. All years have distinct tickers.")

# Drop exact-duplicate firm-year rows (same ticker + year, identical data).
kld <-
  kld |>
  unite(id, ticker, year, remove = FALSE) |>
  distinct(id, .keep_all = TRUE) |>
  select(-id)

# ── save long-format metadata for 03-extract.R ────────────────────────────────

write_csv(
  x    = select(kld, name, ticker, issuerid, year),
  file = "data/intermediate/firm-year-level.csv"
)

# ── reshape to wide: rows = firms, cols = indicator_year ──────────────────────
#
# Each column corresponds to one (indicator, year) pair. Firms absent in a
# given year receive NA rows for that year's columns.

all_tickers <- sort(unique(kld$ticker))

for (k in seq_len(n_years)) {

  yr <- yr_offset + k

  in_year <-
    kld |>
    filter(year == yr, ticker %in% firms_in_year[[k]]) |>
    select(ticker, all_of(mets_by_year[[k]])) |>
    rename_with(~ paste(., yr, sep = "_"), .cols = -ticker) |>
    arrange(ticker)

  out_of_year <- tibble(ticker = setdiff(all_tickers, firms_in_year[[k]]))

  block      <- bind_rows(in_year, out_of_year) |> arrange(ticker) |> select(-ticker)
  year_index <- rep(k, ncol(block))

  if (k == 1L) { dat <- block; yrs <- year_index }
  else         { dat <- bind_cols(dat, block); yrs <- c(yrs, year_index) }
}

# ── final check: all values are 0, 1, or NA ───────────────────────────────────

dat <- mutate(dat, across(everything(), as.integer))

has_bad <- function(x) any(!(is.na(x) | x %in% c(0L, 1L)))

bad_cols <- dat |>
  summarise(across(everything(), has_bad)) |>
  pivot_longer(everything(), names_to = "col", values_to = "flag") |>
  filter(flag) |>
  pull(col)

if (length(bad_cols) > 0) {
  message("Clamping residual non-binary values in: ", paste(bad_cols, collapse = ", "))
  dat <- mutate(dat, across(all_of(bad_cols), clamp_binary))
}

# ── save ──────────────────────────────────────────────────────────────────────

save(dat, yrs, all_tickers,
     file = "data/intermediate/data-reshaped.RData")

message(sprintf("01-clean.R complete: %d firms × %d item-year columns.",
                nrow(dat), ncol(dat)))
