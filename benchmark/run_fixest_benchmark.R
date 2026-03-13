args <- commandArgs(trailingOnly = TRUE)
options <- list(
  input = file.path("benchmark", "artifacts", "benchmark_panel.csv"),
  output = file.path("benchmark", "artifacts", "fixest_results.csv"),
  threads = "0"
)
for (arg in args) {
  if (startsWith(arg, "--") && grepl("=", arg, fixed = TRUE)) {
    pieces <- strsplit(substring(arg, 3), "=", fixed = TRUE)[[1]]
    options[[pieces[1]]] <- pieces[2]
  }
}

if (!requireNamespace("fixest", quietly = TRUE)) {
  stop("fixest is not installed. Run benchmark/setup_fixest.R or benchmark/run_all.ps1 first.")
}

library(fixest)
requested_threads <- as.integer(options$threads)
if (is.na(requested_threads) || requested_threads <= 0) {
  requested_threads <- parallel::detectCores(logical = TRUE)
}
setFixest_nthreads(requested_threads)

input_path <- options$input
output_path <- options$output

df <- read.csv(input_path, na.strings = c("", "NA"))

time_spec <- function(formula_obj, data) {
  invisible(feols(formula_obj, data = data))
  gc()
  elapsed <- system.time(
    invisible(feols(formula_obj, data = data))
  )[["elapsed"]]
  elapsed
}

specs <- list(
  # Match Julia/Stata's omitted periods: -1 as the reference period and -17
  # as the trimmed lead just outside the checked event window.
  list(name = "id + t", formula = Y ~ sunab(first_treat_fixest, t, ref.p = c(-1, -17)) | id + t),
  list(name = "id + id1 + id2 + t", formula = Y ~ sunab(first_treat_fixest, t, ref.p = c(-1, -17)) | id + id1 + id2 + t)
)

results <- data.frame(engine = character(), spec = character(), seconds = double(), nobs = integer())
for (spec in specs) {
  elapsed <- time_spec(spec$formula, df)
  results <- rbind(results, data.frame(
    engine = "R fixest",
    spec = spec$name,
    seconds = elapsed,
    nobs = nrow(df)
  ))
  cat(sprintf("fixest | %s | %.3f seconds\n", spec$name, elapsed))
}

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write.csv(results, output_path, row.names = FALSE)
