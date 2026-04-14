#!/usr/bin/env Rscript
################################################################################
# Run PhyloWGS on all DREAM samples and record algorithm runtime
#
# 1. Read PhyloWGS_runtime_analysis/mutation_table (*_ssm_data.txt) and
#    copy_number_file (*_cnv_data.txt). Samples with both files are run.
# 2. Before each run: if PhyloWGS_timecount.csv exists, skip samples that
#    already have an (algorithm="PhyloWGS", sample=<id>) row.
# 3. For each remaining sample:
#    - Time only the PhyloWGS run (multievolve.py), excluding data prep,
#      result writing, and any plotting.
# 4. Append one row per sample to PhyloWGS_timecount.csv with columns:
#    algorithm, sample, time, diploid_mutation_count, mutation_count, N_clusters
#    (N_clusters set to NA; diploid_mutation_count = mutation_count for PhyloWGS
#     since we do not filter by diploid here).
#
# Prerequisites:
#   - Run create_phylowgs_inputs_from_vcf_battenberg.r first to generate
#     mutation_table/*_ssm_data.txt and copy_number_file/*_cnv_data.txt.
#   - Python 2 with PhyloWGS dependencies (numpy, scipy, ete2, GSL), and
#     multievolve.py/evolve.py run from PhyloWGS repo root.
################################################################################

##========================= Paths =========================##

project_root <- "/Users/khanhngocdinh/Documents/YanjieChen/GitHub/PhyloWGS"
setwd(project_root)

base_dir          <- file.path(project_root, "PhyloWGS_runtime_analysis")
mutation_dir     <- file.path(base_dir, "mutation_table")
copy_number_dir   <- file.path(base_dir, "copy_number_file")
results_dir      <- file.path(base_dir, "results_phylowgs")
timecount_csv    <- file.path(base_dir, "PhyloWGS_timecount.csv")

# PhyloWGS scripts (run from project_root)
# Prefer python2, fall back to python2.7 (e.g. pyenv or system)
python2          <- Sys.which("python2")
if (!nzchar(python2)) python2 <- Sys.which("python2.7")
# If not in PATH, try pyenv shims (e.g. when R is started from IDE)
if (!nzchar(python2)) {
  pyenv_root <- Sys.getenv("PYENV_ROOT", NA_character_)
  if (is.na(pyenv_root) || !nzchar(pyenv_root)) pyenv_root <- path.expand("~/.pyenv")
  shims <- file.path(pyenv_root, "shims")
  if (dir.exists(shims)) {
    old_path <- Sys.getenv("PATH", "")
    Sys.setenv(PATH = paste(shims, old_path, sep = .Platform$path.sep))
    python2 <- Sys.which("python2.7")
    if (!nzchar(python2)) python2 <- Sys.which("python2")
    Sys.setenv(PATH = old_path)
  }
}
if (!nzchar(python2)) {
  stop("Neither 'python2' nor 'python2.7' found in PATH. ",
       "Install Python 2.7 (e.g. pyenv install 2.7.18) or see PhyloWGS_runtime_analysis/INSTALL_PYTHON2.md")
}
python2          <- as.character(python2)
multievolve_script <- file.path(project_root, "multievolve.py")

# MCMC settings (use smaller values for faster test runs)
NUM_CHAINS       <- 4L
BURNIN_SAMPLES   <- 1000L   # default 1000
MCMC_SAMPLES     <- 2500L   # default 2500

if (!dir.exists(base_dir))        dir.create(base_dir, recursive = TRUE)
if (!dir.exists(mutation_dir)) {
  stop("mutation_table directory not found: ", mutation_dir,
       ". Run create_phylowgs_inputs_from_vcf_battenberg.r first.")
}
if (!dir.exists(copy_number_dir)) {
  stop("copy_number_file directory not found: ", copy_number_dir,
       ". Run create_phylowgs_inputs_from_vcf_battenberg.r first.")
}
if (!dir.exists(results_dir))     dir.create(results_dir, recursive = TRUE)
if (!file.exists(multievolve_script)) {
  stop("multievolve.py not found: ", multievolve_script)
}

##========================= Helpers =========================##

#' Discover DREAM samples that have both *_ssm_data.txt and *_cnv_data.txt
find_dream_samples <- function(mut_dir, cnv_dir) {
  ssm_files <- list.files(mut_dir, pattern = "_ssm_data\\.txt$", full.names = FALSE)
  sample_from_ssm <- gsub("_ssm_data\\.txt$", "", ssm_files)
  out <- character(0)
  for (s in sample_from_ssm) {
    cnv_file <- file.path(cnv_dir, paste0(s, "_cnv_data.txt"))
    if (file.exists(cnv_file)) out <- c(out, s)
  }
  sort(unique(out))
}

#' Read existing timecount CSV and return character vector of sample ids
#' that already have a row with algorithm == algorithm_name.
samples_already_run <- function(csv_path, algorithm_name = "PhyloWGS") {
  if (!file.exists(csv_path)) return(character(0))
  tbl <- read.csv(csv_path, stringsAsFactors = FALSE)
  if (!"algorithm" %in% names(tbl) || !"sample" %in% names(tbl)) return(character(0))
  done <- tbl[tbl$algorithm == algorithm_name, "sample"]
  as.character(unique(done))
}

#' Count number of SSM rows (data rows, excluding header) in ssm_data file
count_ssm_rows <- function(ssm_path) {
  if (!file.exists(ssm_path)) return(0L)
  n <- length(readLines(ssm_path))
  max(0L, n - 1L)
}

#' Append one row to timecount CSV (same structure as DECODE/SciClone)
append_timecount <- function(algorithm,
                             sample,
                             time_sec,
                             diploid_mutation_count,
                             mutation_count,
                             N_clusters,
                             csv_path = timecount_csv) {
  new_row <- data.frame(
    algorithm = algorithm,
    sample = sample,
    time = time_sec,
    diploid_mutation_count = diploid_mutation_count,
    mutation_count = mutation_count,
    N_clusters = N_clusters,
    stringsAsFactors = FALSE
  )
  if (file.exists(csv_path)) {
    prev <- read.csv(csv_path, stringsAsFactors = FALSE)
    out  <- rbind(prev, new_row)
  } else {
    out  <- new_row
  }
  write.csv(out, csv_path, row.names = FALSE)
}

##========================= Single-sample PhyloWGS run (timed) =========================##

run_phylowgs_single_dream <- function(sample_id) {
  cat("\n========================================\n")
  cat("PhyloWGS DREAM runtime — sample: ", sample_id, "\n", sep = "")

  ssm_path <- file.path(mutation_dir, paste0(sample_id, "_ssm_data.txt"))
  cnv_path <- file.path(copy_number_dir, paste0(sample_id, "_cnv_data.txt"))

  if (!file.exists(ssm_path)) stop("SSM file not found: ", ssm_path)
  if (!file.exists(cnv_path)) stop("CNV file not found: ", cnv_path)

  mutation_count <- count_ssm_rows(ssm_path)
  if (mutation_count == 0) {
    stop("Sample ", sample_id, " has no SSM rows in ", ssm_path)
  }

  # Output directory for this sample (multievolve creates chain_0, chain_1, ... inside)
  sample_out_dir <- file.path(results_dir, sample_id)
  if (!dir.exists(sample_out_dir)) dir.create(sample_out_dir, recursive = TRUE)

  # Build command: run from project_root so PhyloWGS imports work
  args <- c(
    multievolve_script,
    "--num-chains", as.character(NUM_CHAINS),
    "--ssms", normalizePath(ssm_path, mustWork = TRUE),
    "--cnvs", normalizePath(cnv_path, mustWork = TRUE),
    "-O", normalizePath(sample_out_dir, mustWork = FALSE),
    "--burnin-samples", as.character(BURNIN_SAMPLES),
    "--mcmc-samples", as.character(MCMC_SAMPLES)
  )

  cat("  Running PhyloWGS (timing started)...\n")
  t_start <- Sys.time()
  owd <- getwd()
  on.exit(setwd(owd), add = TRUE)
  setwd(project_root)
  ret <- system2(
    python2,
    args,
    stdout = file.path(sample_out_dir, "phylowgs_stdout.txt"),
    stderr = file.path(sample_out_dir, "phylowgs_stderr.txt")
  )
  setwd(owd)
  on.exit(NULL)
  t_end <- Sys.time()
  elapsed_sec <- as.numeric(difftime(t_end, t_start, units = "secs"))

  if (ret != 0) {
    stop("PhyloWGS exited with code ", ret, ". Check ", sample_out_dir, " for logs.")
  }

  cat("  PhyloWGS completed in ", round(elapsed_sec, 2), " seconds\n", sep = "")

  # PhyloWGS outputs trees, not a single N_clusters; use NA to match CSV structure
  append_timecount(
    algorithm = "PhyloWGS",
    sample = sample_id,
    time_sec = elapsed_sec,
    diploid_mutation_count = mutation_count,
    mutation_count = mutation_count,
    N_clusters = NA_integer_,
    csv_path = timecount_csv
  )
  cat("  PhyloWGS_timecount.csv updated.\n")

  invisible(list(time_sec = elapsed_sec, mutation_count = mutation_count))
}

##========================= Main =========================##

cat("=== PhyloWGS Runtime Analysis for DREAM ===\n")
cat("mutation_table :", mutation_dir, "\n")
cat("copy_number    :", copy_number_dir, "\n")
cat("results        :", results_dir, "\n")
cat("timecount      :", timecount_csv, "\n\n")

all_samples <- find_dream_samples(mutation_dir, copy_number_dir)
if (length(all_samples) == 0) {
  stop("No samples found (need both *_ssm_data.txt and *_cnv_data.txt).")
}

already_done <- samples_already_run(timecount_csv, "PhyloWGS")
to_run <- setdiff(all_samples, already_done)

if (length(already_done) > 0) {
  cat("Samples already in timecount (skipped): ", length(already_done), "\n  ",
      paste(already_done, collapse = ", "), "\n\n", sep = "")
}
cat("Samples to run: ", length(to_run), "\n  ", paste(to_run, collapse = ", "), "\n\n", sep = "")

if (length(to_run) == 0) {
  cat("Nothing to run. Exiting.\n")
  quit(save = "no", status = 0)
}

successful <- character(0)
failed     <- character(0)

for (S in to_run) {
  cat("\n----------------------------------------\n")
  cat("Sample:", S, "\n")
  ssm_path <- file.path(mutation_dir, paste0(S, "_ssm_data.txt"))
  mut_count <- count_ssm_rows(ssm_path)
  if (mut_count > 5000L) {
    cat("  Skipping sample ", S,
        " because mutation_count=", mut_count,
        " (> 5000); expected runtime too long.\n", sep = "")
    next
  }
  res <- tryCatch({
    run_phylowgs_single_dream(S)
    successful <<- c(successful, S)
    TRUE
  }, error = function(e) {
    cat("  ERROR: ", conditionMessage(e), "\n", sep = "")
    failed <<- c(failed, S)
    FALSE
  })
}

successful <- unique(successful)
failed     <- unique(failed)

cat("\n=== Summary (PhyloWGS DREAM runtime) ===\n")
cat("Successful:", length(successful), "\n")
if (length(successful) > 0) cat("  ", paste(successful, collapse = ", "), "\n")
cat("Failed    :", length(failed), "\n")
if (length(failed) > 0) cat("  ", paste(failed, collapse = ", "), "\n")
cat("Time counts written to:", timecount_csv, "\n")
