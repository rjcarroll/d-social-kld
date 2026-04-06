# convert-to-parquet.R
#
# Converts chain CSVs to Parquet format for fast column-selective reads.
# Drops NCP internals (z_init, z_trans, sigma_firm_raw) that are never used
# downstream, cutting the column count from ~120K to ~56K and file size roughly
# in half. Processes one chain at a time (~1 GB peak memory).

library(arrow)

csv_files <- sort(list.files("data/intermediate/chains", pattern = "\\.csv$",
                             full.names = TRUE))
csv_files <- csv_files[!grepl("profile", csv_files)]
stopifnot("No chain CSVs found" = length(csv_files) > 0)

# Parse header to identify columns to keep
hdr <- strsplit(
  system(sprintf("grep -v '^#' '%s' | head -1", csv_files[1]), intern = TRUE),
  ","
)[[1]]

# Keep: sampler internals, anchors, alpha, beta, sigma_firm, sigma_global,
#        sigma_init, theta. Drop: z_init, z_trans, sigma_firm_raw.
keep <- !grepl("^(z_init|z_trans|sigma_firm_raw)\\.", hdr)
keep_names <- hdr[keep]

cat(sprintf("Total columns: %d\n", length(hdr)))
cat(sprintf("Keeping: %d  Dropping: %d\n", sum(keep), sum(!keep)))

dir.create("data/intermediate/chains-parquet", showWarnings = FALSE)

for (i in seq_along(csv_files)) {
  cat(sprintf("\nChain %d/%d: %s\n", i, length(csv_files), basename(csv_files[i])))

  # Strip comment lines to a temp file
  tmp <- tempfile(fileext = ".csv")
  system2("grep", args = c("-v", "^#", csv_files[i]), stdout = tmp)

  # Read with arrow (128 MB block size to accommodate the 1.3 MB header)
  tbl <- read_csv_arrow(tmp,
    read_options = CsvReadOptions$create(block_size = 128L * 1024L * 1024L))
  cat(sprintf("  Read: %d rows x %d cols\n", nrow(tbl), ncol(tbl)))

  # Select only the columns we need
  tbl <- tbl[, keep_names]
  cat(sprintf("  After drop: %d cols\n", ncol(tbl)))

  # Write as parquet
  out_name <- sub("\\.csv$", ".parquet", basename(csv_files[i]))
  out_path <- file.path("data/intermediate/chains-parquet", out_name)
  write_parquet(tbl, out_path)
  cat(sprintf("  Wrote: %s (%.2f GB)\n", out_name, file.size(out_path) / 1e9))

  unlink(tmp)
  rm(tbl); gc(verbose = FALSE)
}

cat("\nDone. Parquet files:\n")
for (f in list.files("data/intermediate/chains-parquet", full.names = TRUE)) {
  cat(sprintf("  %s  %.2f GB\n", basename(f), file.size(f) / 1e9))
}
