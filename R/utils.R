# utils.R
#
# Shared helper functions for the D-SOCIAL-KLD pipeline.
# Sourced by 01-clean.R and 02-estimate.R as needed.


# ── Data-cleaning helpers ─────────────────────────────────────────────────────

#' Coerce a logical column to 0/1 integer, preserving NA.
logical_to_int <- function(x) {
  y <- integer(length(x))
  y[is.na(x)]    <- NA_integer_
  y[x == FALSE]  <- 0L
  y[x == TRUE]   <- 1L
  y
}

#' Clamp values ≥ 2 down to 1, coercing to integer. Preserves NA.
clamp_binary <- function(x) {
  y <- as.integer(x)
  y[!is.na(y) & y >= 2L] <- 1L
  y
}

#' Return TRUE if a column has any non-NA values.
has_data <- function(x) any(!is.na(x))

#' Return TRUE if all elements of x are distinct (no duplicates).
all_distinct <- function(x) identical(sort(x), sort(unique(x)))


# ── Duplicate-ticker resolution ───────────────────────────────────────────────
#
# Several firms from distinct KLD coverage universes share the same ticker. The
# strategy is:
#
#   1. For years with an `issuerid` column populated (2017–2018 in practice):
#      use issuerid to assign deterministic, stable firm keys. Firms that share
#      a ticker but have different issuerids are kept separate.
#
#   2. For earlier years (issuerid absent): the firm that *first* appears under
#      a ticker keeps it; later entrants receive a numeric suffix: TICK.2,
#      TICK.3, etc., ordered by first appearance year. This is deterministic
#      and collision-free by construction.
#
# A log is maintained so that once a firm receives a suffixed ticker it carries
# that ticker forward into all subsequent years.

#' Resolve duplicate tickers within a single year using deterministic suffixes.
#'
#' @param kld       Full KLD tibble (columns: name, ticker, year, ...).
#' @param yr        Integer: the year being processed.
#' @param log       Tibble: running disambiguation log with columns
#'                  `name`, `old_ticker`, `new_ticker`, `year_created`,
#'                  `year_last_used`.
#'
#' @return A list with `kld` (updated) and `log` (updated).
resolve_duplicate_tickers <- function(kld, yr, log) {

  # Find tickers that appear for more than one distinct firm name in this year.
  collisions <-
    dplyr::filter(kld, year == yr) |>
    dplyr::distinct(ticker, name) |>
    dplyr::group_by(ticker) |>
    dplyr::filter(dplyr::n() > 1) |>
    dplyr::ungroup()

  if (nrow(collisions) == 0L) return(list(kld = kld, log = log))

  for (tick_i in unique(collisions$ticker)) {

    # Names already seen under this ticker in earlier years.
    historical <-
      dplyr::filter(kld, year < yr, ticker == tick_i) |>
      dplyr::distinct(name) |>
      dplyr::pull(name)

    # All distinct names claiming this ticker in the current year.
    current <-
      dplyr::filter(kld, year == yr, ticker == tick_i) |>
      dplyr::distinct(name) |>
      dplyr::pull(name)

    # New firms: not seen before under this ticker, not already in the log.
    already_logged <- dplyr::filter(log, old_ticker == tick_i) |> dplyr::pull(name)
    new_names <- setdiff(current, c(historical, already_logged))

    if (length(new_names) == 0L) next

    # Assign the next available numeric suffix for each new firm.
    # Suffixes already used for this base ticker (from the log).
    used_suffixes <-
      dplyr::filter(log, old_ticker == tick_i) |>
      dplyr::pull(new_ticker) |>
      stringr::str_extract("(?<=\\.)[0-9]+$") |>
      as.integer() |>
      stats::na.omit()

    next_suffix <- if (length(used_suffixes) == 0L) 2L else max(used_suffixes) + 1L

    new_tickers <- paste0(tick_i, ".", seq(next_suffix, next_suffix + length(new_names) - 1L))

    new_rows <- dplyr::tibble(
      name           = new_names,
      old_ticker     = tick_i,
      new_ticker     = new_tickers,
      year_created   = yr,
      year_last_used = yr
    )
    log <- dplyr::bind_rows(log, new_rows)

    for (j in seq_along(new_names)) {
      kld$ticker[kld$year == yr & kld$name == new_names[j]] <- new_tickers[j]
    }
  }

  list(kld = kld, log = log)
}

#' Propagate previously assigned suffixed tickers forward into a new year.
#'
#' @param kld Full KLD tibble.
#' @param log Disambiguation log from `resolve_duplicate_tickers()`.
#' @param yr  Year to update.
#'
#' @return List with updated `kld` and `log` (`year_last_used` refreshed).
propagate_ticker_log <- function(kld, log, yr) {

  for (i in seq_len(nrow(log))) {
    rows_i <- kld$year == yr & kld$name == log$name[i]
    if (any(rows_i)) {
      kld$ticker[rows_i]    <- log$new_ticker[i]
      log$year_last_used[i] <- yr
    }
  }

  list(kld = kld, log = log)
}

#' Resolve duplicate tickers using issuerid where available.
#'
#' When `issuerid` is populated for firms in a given year, two firms sharing
#' a ticker but with different issuerids are unambiguously distinct and receive
#' deterministic suffixes derived from their issuerid.
#'
#' @param kld Full KLD tibble (must have an `issuerid` column).
#' @param yr  Integer: year to process.
#' @param log Disambiguation log.
#'
#' @return List with updated `kld` and `log`.
resolve_with_issuerid <- function(kld, yr, log) {

  yr_data <- dplyr::filter(kld, year == yr, !is.na(issuerid), issuerid != "")

  if (nrow(yr_data) == 0L) return(list(kld = kld, log = log))

  # Find tickers with multiple distinct issuerids in this year.
  collisions <-
    yr_data |>
    dplyr::distinct(ticker, issuerid, name) |>
    dplyr::group_by(ticker) |>
    dplyr::filter(dplyr::n() > 1) |>
    dplyr::ungroup()

  if (nrow(collisions) == 0L) return(list(kld = kld, log = log))

  for (tick_i in unique(collisions$ticker)) {

    entries <- dplyr::filter(collisions, ticker == tick_i) |>
      dplyr::arrange(issuerid)   # deterministic order

    # The firm with the lowest issuerid (alphabetically) keeps the original
    # ticker; the rest receive .2, .3, … suffixes.
    historical_issuerids <-
      dplyr::filter(kld, year < yr, ticker == tick_i) |>
      dplyr::pull(issuerid) |>
      stats::na.omit() |>
      unique()

    keeper <- if (length(historical_issuerids) > 0) {
      dplyr::filter(entries, issuerid %in% historical_issuerids)$issuerid[1]
    } else {
      entries$issuerid[1]
    }

    new_entries <- dplyr::filter(entries, issuerid != keeper)
    if (nrow(new_entries) == 0L) next

    used_suffixes <-
      dplyr::filter(log, old_ticker == tick_i) |>
      dplyr::pull(new_ticker) |>
      stringr::str_extract("(?<=\\.)[0-9]+$") |>
      as.integer() |>
      stats::na.omit()

    next_suffix <- if (length(used_suffixes) == 0L) 2L else max(used_suffixes) + 1L

    for (j in seq_len(nrow(new_entries))) {
      new_tick <- paste0(tick_i, ".", next_suffix + j - 1L)
      iid_j    <- new_entries$issuerid[j]
      name_j   <- new_entries$name[j]

      kld$ticker[kld$year == yr & kld$issuerid == iid_j] <- new_tick

      log <- dplyr::bind_rows(log, dplyr::tibble(
        name           = name_j,
        old_ticker     = tick_i,
        new_ticker     = new_tick,
        year_created   = yr,
        year_last_used = yr
      ))
    }
  }

  list(kld = kld, log = log)
}


# ── Stan data preparation helpers ─────────────────────────────────────────────

#' Convert the wide firm × item matrix to a long data frame of non-missing obs.
#'
#' @param wide_mat  A numeric matrix (rows = firms, cols = items). Row names
#'                  are firm tickers; column names are "INDICATOR_YEAR".
#' @param yrs       Integer vector of length ncol(wide_mat): the time-period
#'                  index (1-based) for each column.
#'
#' @return A tibble with columns: ticker, item_col, year_idx, y (0/1).
wide_to_long <- function(wide_mat, yrs) {

  tickers   <- rownames(wide_mat)
  item_cols <- colnames(wide_mat)

  rows <- vector("list", ncol(wide_mat))
  for (j in seq_len(ncol(wide_mat))) {
    col_j <- wide_mat[, j]
    obs_i <- which(!is.na(col_j))
    if (length(obs_i) == 0L) next
    rows[[j]] <- dplyr::tibble(
      ticker   = tickers[obs_i],
      item_col = item_cols[j],
      year_idx = yrs[j],
      y        = as.integer(col_j[obs_i])
    )
  }

  dplyr::bind_rows(rows)
}

#' Build the consecutive-year transition arrays needed by the Stan model.
#'
#' For each firm, identify pairs of firm-year entries (prev, curr) where the
#' curr year immediately follows prev in the data for that firm.
#'
#' @param fy_tbl  A tibble with columns `fy_idx` (integer), `ticker` (chr),
#'                `year_idx` (integer), `firm_idx` (integer).
#'
#' @return A tibble with columns: trans_prev (fy_idx), trans_curr (fy_idx),
#'         trans_firm (firm_idx).
build_transitions <- function(fy_tbl) {

  fy_tbl |>
    dplyr::arrange(ticker, year_idx) |>
    dplyr::group_by(ticker) |>
    dplyr::mutate(
      prev_fy_idx  = dplyr::lag(fy_idx),
      prev_year    = dplyr::lag(year_idx),
      is_consec    = !is.na(prev_year) & (year_idx == prev_year + 1L)
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(is_consec) |>
    dplyr::transmute(
      trans_prev = prev_fy_idx,
      trans_curr = fy_idx,
      trans_firm = firm_idx
    )
}
