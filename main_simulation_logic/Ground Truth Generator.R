
library(data.table)

### load ground truth ######
generate_pbmc_tirtl_long <- function(file_path, verbose = TRUE) {
  if (verbose) {
    cat("=== generate_pbmc_tirtl_long() ===\n")
    cat("Reading RDS file from:\n  ", file_path, "\n\n", sep = "")
  }
  
  PBMC_tirtl_cleaned_df <- readRDS(file_path)
  setDT(PBMC_tirtl_cleaned_df)
  
  if (verbose) {
    cat("Loaded object has", nrow(PBMC_tirtl_cleaned_df), "rows.\n")
    cat("Columns:\n")
    print(names(PBMC_tirtl_cleaned_df))
    cat("\n")
  }
  
  required_cols <- c("clonotype_id", "alpha_nuc", "beta_nuc",
                     "cdr3a", "cdr3b", "madhype_score", "tshell_padj",
                     "wells", "intensities", "rel_intensities")
  missing_cols <- setdiff(required_cols, names(PBMC_tirtl_cleaned_df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in PBMC_tirtl_cleaned_df: ",
         paste(missing_cols, collapse = ", "))
  }
  
  if (verbose) cat("Step 1: Expanding clonotypes into per-well rows...\n")
  
  PBMC_tirtl_long <- PBMC_tirtl_cleaned_df[
    , {
      ws   <- trimws(strsplit(wells, ",")[[1]])
      ints <- as.numeric(trimws(strsplit(intensities, ",")[[1]]))
      rel  <- as.numeric(trimws(strsplit(rel_intensities, ",")[[1]]))
      
      # handle length mismatches
      len <- max(length(ws), length(ints), length(rel))
      length(ws)   <- len
      length(ints) <- len
      length(rel)  <- len
      
      data.table(
        clonotype_id   = clonotype_id,
        alpha_nuc      = alpha_nuc,
        beta_nuc       = beta_nuc,
        cdr3a          = cdr3a,
        cdr3b          = cdr3b,
        madhype_score  = madhype_score,
        tshell_padj    = tshell_padj,
        well           = ws,
        intensity      = ints,   # NEW: raw TIRTL read counts per clonotype–well
        rel_intensity  = rel
      )
    },
    by = clonotype_id
  ][
    , clonotype_well := paste0("clone", clonotype_id, "_", well)
  ][
    , .(clonotype_well, clonotype_id,
        alpha_nuc, beta_nuc, cdr3a, cdr3b,
        madhype_score, tshell_padj,
        well, intensity, rel_intensity)
  ]
  
  if (verbose) {
    cat("Step 2: Done generating PBMC_tirtl_long.\n")
    cat("  Rows:", nrow(PBMC_tirtl_long), "\n")
    cat("  Unique wells:", length(unique(PBMC_tirtl_long$well)), "\n")
    cat("  Unique clonotypes:", length(unique(PBMC_tirtl_long$clonotype_id)), "\n\n")
    cat("Preview (head):\n")
    print(head(PBMC_tirtl_long))
    cat("\n=== generate_pbmc_tirtl_long() finished ===\n\n")
  }
  
  return(PBMC_tirtl_long)
}



############# simulate on wells

simulate_cxcl13_on_wells <- function(PBMC_tirtl_long,
                                     subset_rows        = NULL,   # e.g. c("A","B","C","D")
                                     subset_wells       = NULL,   # e.g. c("A11","B5")
                                     # log2-scale parameters for clonotype-level CXCL13
                                     nonreact_mean_log2 = 0.3,    # ~2 units, near-background
                                     nonreact_sd_log2   = 0.4,
                                     react_mean_log2    = 3.5,    # ~45 units, ~25x higher
                                     react_sd_log2      = 0.7,
                                     reactive_fraction  = 0.01,   # fraction of clonotypes set as reactive
                                     # per-row noise on log2 scale (multiplicative noise)
                                     per_row_noise_sd_log2 = 0.3,
                                     seed_clono_noise   = 777,
                                     seed_row_noise     = 12,
                                     verbose            = TRUE) {
  dt <- as.data.table(PBMC_tirtl_long)
  
  if (verbose) {
    cat("=== simulate_cxcl13_on_wells() ===\n")
    cat("Input rows:", nrow(dt), "\n")
    cat("Unique wells before subsetting:", length(unique(dt$well)), "\n\n")
  }
  
  # -------- 1) Subsetting wells (optional) --------
  if (!is.null(subset_rows)) {
    if (verbose) {
      cat("Subsetting to rows (first letter of well) in: ",
          paste(subset_rows, collapse = ", "), "\n")
    }
    dt <- dt[substr(well, 1, 1) %in% subset_rows]
  }
  
  if (!is.null(subset_wells)) {
    if (verbose) {
      cat("Further subsetting to specific wells in: ",
          paste(subset_wells, collapse = ", "), "\n")
    }
    dt <- dt[well %in% subset_wells]
  }
  
  if (verbose) {
    cat("After subsetting:\n")
    cat("  Rows:", nrow(dt), "\n")
    cat("  Unique wells:", length(unique(dt$well)), "\n")
    cat("  Unique clonotypes:", length(unique(dt$clonotype_id)), "\n\n")
  }
  
  if (nrow(dt) == 0) {
    stop("No rows left after subsetting. Check subset_rows/subset_wells.")
  }
  
  # -------- 2) Determine which clonotypes are reactive --------
  if (verbose) {
    cat("Step 2: Selecting reactive clonotypes by highest total rel_intensity...\n")
    cat("  Reactive fraction requested:", reactive_fraction, "\n")
    cat("  Non-reactive log2 mean, sd: ", nonreact_mean_log2, ", ", nonreact_sd_log2, "\n", sep = "")
    cat("  Reactive log2 mean, sd:     ", react_mean_log2,    ", ", react_sd_log2,    "\n", sep = "")
  }
  
  reactive_fraction <- max(min(reactive_fraction, 1), 0)
  
  clonos  <- sort(unique(dt$clonotype_id))
  n_clono <- length(clonos)
  
  if (n_clono == 0) {
    stop("No clonotypes present after subsetting.")
  }
  
  # Rank clonotypes by total rel_intensity (descending)
  rank_dt <- dt[, .(total_rel_intensity = sum(rel_intensity, na.rm = TRUE)),
                by = clonotype_id]
  setorder(rank_dt, -total_rel_intensity, clonotype_id)
  
  n_reactive <- floor(reactive_fraction * n_clono)
  if (n_reactive > 0) {
    reactive_ids <- rank_dt$clonotype_id[seq_len(n_reactive)]
  } else {
    reactive_ids <- integer(0)
  }
  
  is_reactive_vec <- clonos %in% reactive_ids
  
  if (verbose) {
    cat("  Total clonotypes:", n_clono, "\n")
    cat("  Number set as reactive:", n_reactive, "\n")
    if (n_reactive > 0) {
      cat("  Example reactive clonotype_ids (up to 10):\n")
      print(head(reactive_ids, 10))
    }
    cat("\n")
  }
  
  # -------- 3) Simulate per-clonotype CXCL13 on log2 scale --------
  set.seed(seed_clono_noise)
  
  # OLD normal-based simulation (kept for reference, now commented out):
  # nonreact_mean <- 6.9
  # nonreact_sd   <- 26.4
  # react_mean    <- 21.6
  # react_sd      <- 32.4
  #
  # cx_vals <- numeric(n_clono)
  # cx_vals[!is_reactive_vec] <- rnorm(sum(!is_reactive_vec),
  #                                    mean = nonreact_mean, sd = nonreact_sd)
  # if (n_reactive > 0) {
  #   cx_vals[is_reactive_vec] <- rnorm(sum(is_reactive_vec),
  #                                     mean = react_mean, sd = react_sd)
  # }
  # cx_vals <- pmax(0, cx_vals)  # clamp at 0
  
  # NEW log-normal style simulation:
  log2_vals <- numeric(n_clono)
  
  # Non-reactive (resting / low CXCL13 TCRs)
  log2_vals[!is_reactive_vec] <- rnorm(
    n    = sum(!is_reactive_vec),
    mean = nonreact_mean_log2,
    sd   = nonreact_sd_log2
  )
  
  # Reactive (CXCL13-high, activated / tumor-reactive TCRs)
  if (n_reactive > 0) {
    log2_vals[is_reactive_vec] <- rnorm(
      n    = sum(is_reactive_vec),
      mean = react_mean_log2,
      sd   = react_sd_log2
    )
  }
  
  # Convert log2 expression to linear scale (always positive)
  cx_vals <- 2^log2_vals
  
  CXCL13_simulated <- data.table(
    clonotype_id             = clonos,
    CXCL13_simulated_exp_raw = cx_vals,
    is_reactive              = is_reactive_vec
  )
  
  if (verbose) {
    cat("  Summary of CXCL13_simulated_exp_raw (linear scale, all clonotypes):\n")
    print(summary(CXCL13_simulated$CXCL13_simulated_exp_raw))
    cat("  Reactive vs non-reactive counts:\n")
    print(table(CXCL13_simulated$is_reactive))
    cat("\n")
  }
  
  # Attach per-clonotype CXCL13 + reactivity flag to each clonotype–well row
  dt[
    CXCL13_simulated,
    `:=`(
      CXCL13_simulated_exp_raw = i.CXCL13_simulated_exp_raw,
      is_reactive              = i.is_reactive
    ),
    on = "clonotype_id"
  ]
  
  # -------- 4) Add per-row noise (multiplicative, on log2 scale) --------
  if (verbose) {
    cat("Step 3: Adding per-row noise (log-normal multiplicative) to create CXCL13_simulated_exp...\n")
    cat("  Per-row noise SD on log2 scale:", per_row_noise_sd_log2, "\n")
    cat("  Seed for row noise:", seed_row_noise, "\n\n")
  }
  
  set.seed(seed_row_noise)
  
  # Noise on log2 scale ~ N(0, per_row_noise_sd_log2^2)
  row_noise_log2 <- rnorm(n = nrow(dt), mean = 0, sd = per_row_noise_sd_log2)
  
  # Expression per row = clonotype-level expression * 2^(noise)
  dt[, CXCL13_simulated_exp := CXCL13_simulated_exp_raw * 2^row_noise_log2]
  
  if (verbose) {
    cat("  Summary of CXCL13_simulated_exp (row-level, linear scale):\n")
    print(summary(dt$CXCL13_simulated_exp))
    cat("\nPreview (head):\n")
    print(head(dt))
    cat("\n=== simulate_cxcl13_on_wells() finished ===\n\n")
  }
  
  return(dt)
}





############### Function calls #####################

file_path <- "C:\\Users\\EKN\\Desktop\\Master_Semester2\\Praktikum_DKFZ_immunology\\Datasets\\paired_clonotypes_PBMCplate_Tirtl.rds"

PBMC_tirtl_long2 <- generate_pbmc_tirtl_long(file_path, verbose = TRUE)

PBMC_tirtl_long_sim2 <- simulate_cxcl13_on_wells(
  PBMC_tirtl_long,
  subset_rows  = NULL,
  subset_wells = NULL,
  verbose      = TRUE
)


PBMC_tirtl_long_sim2_rows_AD <- simulate_cxcl13_on_wells(
  PBMC_tirtl_long2,
  subset_rows  = c("A", "B", "C", "D"),
  subset_wells = NULL,
  verbose      = TRUE
)


