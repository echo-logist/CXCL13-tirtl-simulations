## ============================================================
## 0. Setup
## ============================================================
library(data.table)
library(nnls)

## ============================================================
## 1. Seeded version of simulate_cxcl13_system()
##    (only change: measurement_seed argument + set.seed)
## ============================================================
simulate_cxcl13_system_seeded <- function(PBMC_tirtl_long,
                                          measurement_sd   = 5,
                                          measurement_seed = 2025,
                                          verbose          = TRUE) {
  PBMC_tirtl_long <- as.data.table(PBMC_tirtl_long)
  
  if (!"CXCL13_simulated_exp_raw" %in% names(PBMC_tirtl_long)) {
    stop("Column 'CXCL13_simulated_exp_raw' not found in PBMC_tirtl_long. ",
         "Run simulate_cxcl13_on_wells() first.")
  }
  
  if (verbose) {
    cat("=== simulate_cxcl13_system_seeded() ===\n")
    cat("Input rows:", nrow(PBMC_tirtl_long), "\n")
    cat("Number of unique wells:", length(unique(PBMC_tirtl_long$well)), "\n")
    cat("Number of unique clonotypes:", length(unique(PBMC_tirtl_long$clonotype_id)), "\n\n")
  }
  
  ## 1) Design matrix A (wells x clonotypes)
  if (verbose) cat("Step 1: Building design matrix A (wells x clonotypes)...\n")
  
  A_dt <- dcast(
    PBMC_tirtl_long,
    well ~ clonotype_id,
    value.var    = "rel_intensity",
    fun.aggregate = sum,
    fill         = 0
  )
  
  wells <- A_dt$well
  A     <- as.matrix(A_dt[, -"well"])
  colnames(A) <- paste0("clono", colnames(A))
  
  if (verbose) {
    cat("  Design matrix dimensions: ", nrow(A), " wells x ", ncol(A), " clonotypes\n", sep = "")
    cat("  First 5 wells:\n"); print(head(wells, 5))
    cat("  First 5 clonotype columns:\n"); print(head(colnames(A), 5))
    cat("\n")
  }
  
  ## 2) Extract true CXCL13 + reactive flags
  if (verbose) cat("Step 2: Extracting true per-clonotype CXCL13 values...\n")
  
  cx_per_clone <- unique(
    PBMC_tirtl_long[, .(clonotype_id, CXCL13_simulated_exp_raw)]
  )
  
  if ("is_reactive" %in% names(PBMC_tirtl_long)) {
    reactive_per_clone <- unique(
      PBMC_tirtl_long[, .(clonotype_id, is_reactive)]
    )
  } else {
    reactive_per_clone <- cx_per_clone[, .(clonotype_id, is_reactive = FALSE)]
  }
  
  clone_ids_in_A <- as.integer(sub("^clono", "", colnames(A)))
  
  cx_true <- cx_per_clone[
    match(clone_ids_in_A, clonotype_id),
    CXCL13_simulated_exp_raw
  ]
  
  is_reactive_ordered <- reactive_per_clone[
    match(clone_ids_in_A, clonotype_id),
    is_reactive
  ]
  
  if (verbose) {
    cat("  Length of cx_true:", length(cx_true), "\n")
    cat("  Summary of true CXCL13 per clonotype:\n")
    print(summary(cx_true))
    cat("  Reactive vs non-reactive (if available):\n")
    print(table(is_reactive_ordered, useNA = "ifany"))
    cat("\n")
  }
  
  ## 3) y_true using row-level noisy CXCL13_simulated_exp
  if (verbose) cat("Step 3: Computing y_true using CXCL13_simulated_exp...\n")
  
  if (!"CXCL13_simulated_exp" %in% names(PBMC_tirtl_long)) {
    stop("Column 'CXCL13_simulated_exp' not found in PBMC_tirtl_long. ",
         "Did you run simulate_cxcl13_on_wells() with per-row noise?")
  }
  
  PBMC_tirtl_long[, contrib := rel_intensity * CXCL13_simulated_exp]
  
  y_true_dt <- PBMC_tirtl_long[
    , .(y_true = sum(contrib, na.rm = TRUE)),
    by = well
  ]
  
  y_true <- y_true_dt[match(wells, well), y_true]
  
  if (verbose) {
    cat("  Summary of y_true:\n")
    print(summary(y_true))
    cat("\n")
  }
  
  ## 4) Add measurement noise
  if (verbose) cat("Step 4: Adding measurement noise to get y_obs...\n")
  set.seed(measurement_seed)
  y_obs <- y_true + rnorm(length(y_true), mean = 0, sd = measurement_sd)
  names(y_obs) <- wells
  
  if (verbose) {
    cat("  Measurement noise SD:", measurement_sd, "\n")
    cat("  Summary of y_obs:\n")
    print(summary(y_obs))
    cat("\n")
    cat("=== Simulation complete. Returning A, wells, cx_true, is_reactive, y_true, y_obs ===\n\n")
  }
  
  return(list(
    A           = A,
    wells       = wells,
    cx_true     = cx_true,
    is_reactive = is_reactive_ordered,
    y_true      = y_true,
    y_obs       = y_obs
  ))
}

## ============================================================
## 2. Fixed QC plot with rel_tol = 0.90 and smaller margins
## ============================================================
plot_cxcl13_true_vs_est <- function(cx_true,
                                    cx_hat,
                                    is_reactive = NULL,
                                    rel_tol    = 0.90,  # <- new default
                                    main       = "Comparison of Ground Truth and Estimated CXCL13 Levels per Clonotype") {
  ## Keep only finite values
  ok <- is.finite(cx_true) & is.finite(cx_hat)
  x  <- cx_true[ok]
  y  <- cx_hat[ok]
  
  if (!is.null(is_reactive)) {
    is_reactive <- is_reactive[ok]
  } else {
    is_reactive <- rep(FALSE, length(x))
  }
  
  ## Relative error, with safe denominator
  denom      <- pmax(abs(x), 1e-6)
  rel_err    <- abs(y - x) / denom
  within_tol <- rel_err <= rel_tol
  
  ## Colors for 4 categories
  cols <- ifelse(is_reactive & within_tol, "red3",
                 ifelse(is_reactive & !within_tol, "lightcoral",
                        ifelse(!is_reactive & within_tol, "black", "grey50")))
  
  ## Save old par and use smaller margins for multi-panel plotting
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mar = c(4, 4, 3, 6))
  
  ## Scatterplot
  plot(
    x, y,
    col      = cols,
    pch      = 16,
    xlab     = "True CXCL13 Expression per Clonotype",
    ylab     = "Estimated CXCL13 Expression per Clonotype (NNLS)",
    main     = main,
    cex.main = 0.9,
    cex.lab  = 0.9,
    cex.axis = 0.8
  )
  abline(0, 1, lwd = 2)
  
  ## Pearson r
  r_val <- suppressWarnings(cor(x, y))
  subtitle_text <- sprintf(
    "Pearson correlation (r) = %.3f   (n = %d)",
    r_val, length(x)
  )
  title(sub = subtitle_text, cex.sub = 0.8)
  
  ## Legend in right margin
  par(xpd = NA)
  usr <- par("usr")
  legend_x <- usr[2] + (usr[2] - usr[1]) * 0.05
  legend_y <- usr[4]
  
  legend(
    x = legend_x,
    y = legend_y,
    legend = c(
      "Reactive, correctly estimated",
      "Reactive",
      "Non-reactive, correctly estimated",
      "Non-reactive"
    ),
    col = c("red3", "lightcoral", "black", "grey50"),
    pch = 16,
    bty = "n",
    cex = 0.8
  )
  
  invisible(list(
    r          = r_val,
    rel_error  = rel_err,
    within_tol = within_tol
  ))
}

## ============================================================
## 3. Helper: run one simulation + fit + plot
##    (safe to length mismatches)
## ============================================================
run_and_plot_once <- function(PBMC_tirtl_long,
                              subset_rows           = c("A","B","C","D"),
                              per_row_noise_sd_log2 = 0,
                              measurement_sd        = 0,
                              main_title            = "",
                              verbose               = FALSE) {
  
  ## 1) Simulate clonotype + (optional) row noise (Noise 3 + Noise 1)
  sim_dt <- simulate_cxcl13_on_wells(
    PBMC_tirtl_long,
    subset_rows           = subset_rows,
    subset_wells          = NULL,
    per_row_noise_sd_log2 = per_row_noise_sd_log2,
    verbose               = verbose
  )
  
  ## 2) Simulate well-level system with measurement noise (Noise 2)
  sim_sys <- simulate_cxcl13_system_seeded(
    sim_dt,
    measurement_sd   = measurement_sd,
    measurement_seed = sample.int(.Machine$integer.max, 1),
    verbose          = verbose
  )
  
  ## 3) Compute depth vector r_w (for weights)
  depth_dt <- sim_dt[
    , .(r_w = sum(intensity, na.rm = TRUE)),
    by = well
  ]
  r_w_vec <- depth_dt[match(sim_sys$wells, well), r_w]
  
  ## 4) Fit NNLS
  fit <- fit_cxcl13_nnls(
    A              = sim_sys$A,
    y_well         = sim_sys$y_obs,
    use_background = TRUE,
    use_weights    = TRUE,
    r_w            = r_w_vec,
    cx_true        = sim_sys$cx_true,
    verbose        = verbose
  )
  
  ## 5) Align lengths before plotting (avoid recycling warnings)
  n_common <- min(length(sim_sys$cx_true), length(fit$cx_hat))
  cx_true_plot     <- sim_sys$cx_true[seq_len(n_common)]
  cx_hat_plot      <- fit$cx_hat[seq_len(n_common)]
  is_reactive_plot <- if (!is.null(sim_sys$is_reactive)) {
    sim_sys$is_reactive[seq_len(n_common)]
  } else {
    NULL
  }
  
  ## 6) QC plot (rel_tol default 0.90)
  plot_res <- plot_cxcl13_true_vs_est(
    cx_true     = cx_true_plot,
    cx_hat      = cx_hat_plot,
    is_reactive = is_reactive_plot,
    main        = main_title
  )
  
  invisible(list(
    sim_dt   = sim_dt,
    sim_sys  = sim_sys,
    fit      = fit,
    plot_res = plot_res
  ))
}

## ============================================================
## 4. PLOTS: three situations
## ============================================================

## 4.1 Only Noise 3 (no extra noise)
par(mfrow = c(1,1))
res_n3_only <- run_and_plot_once(
  PBMC_tirtl_long        = PBMC_tirtl_long2,
  per_row_noise_sd_log2  = 0,
  measurement_sd         = 0,
  main_title             = "Only x_c variability (No noise introduced)"
)
par(mfrow = c(1,1))

## 4.2 Noise 2 + Noise 3: vary measurement noise (Noise 2)
measurement_grid <- c(0, 0.5, 1, 2, 5)

par(mfrow = c(2, 3))   # up to 6 panels
noise2_plots <- lapply(measurement_grid, function(msd) {
  run_and_plot_once(
    PBMC_tirtl_long        = PBMC_tirtl_long2,
    per_row_noise_sd_log2  = 0,      # Noise 1 OFF
    measurement_sd         = msd,    # Noise 2 grid
    main_title             = paste0("Noise 2 (measurement SD) = ", msd)
  )
})
par(mfrow = c(1,1))

## 4.3 Noise 1 + Noise 3: vary row noise (Noise 1)
row_noise_grid <- c(0, 0.1, 0.3, 0.5)

par(mfrow = c(2, 2))
noise1_plots <- lapply(row_noise_grid, function(sdl2) {
  run_and_plot_once(
    PBMC_tirtl_long        = PBMC_tirtl_long2,
    per_row_noise_sd_log2  = sdl2,   # Noise 1 grid
    measurement_sd         = 0,      # Noise 2 OFF
    main_title             = paste0("Noise 1 (row SD log2) = ", sdl2)
  )
})
par(mfrow = c(1,1))

## ============================================================
## 5. run_one_sim(): returns correlation + recall
## ============================================================
run_one_sim <- function(PBMC_tirtl_long,
                        subset_rows           = c("A","B","C","D"),
                        per_row_noise_sd_log2 = 0,
                        measurement_sd        = 0,
                        reactive_fraction     = 0.01,
                        nonreact_mean_log2    = 0.3,
                        nonreact_sd_log2      = 0.4,
                        react_mean_log2       = 3.5,
                        react_sd_log2         = 0.7,
                        verbose               = FALSE) {
  
  ## 1) Simulate clonotype + row noise + x_c variability
  sim_dt <- simulate_cxcl13_on_wells(
    PBMC_tirtl_long,
    subset_rows           = subset_rows,
    subset_wells          = NULL,
    nonreact_mean_log2    = nonreact_mean_log2,
    nonreact_sd_log2      = nonreact_sd_log2,
    react_mean_log2       = react_mean_log2,
    react_sd_log2         = react_sd_log2,
    reactive_fraction     = reactive_fraction,
    per_row_noise_sd_log2 = per_row_noise_sd_log2,
    verbose               = verbose
  )
  
  ## 2) Simulate system with measurement noise
  sim_sys <- simulate_cxcl13_system_seeded(
    sim_dt,
    measurement_sd   = measurement_sd,
    measurement_seed = sample.int(.Machine$integer.max, 1),
    verbose          = verbose
  )
  
  ## 3) Depth vector
  depth_dt <- sim_dt[
    , .(r_w = sum(intensity, na.rm = TRUE)),
    by = well
  ]
  r_w_vec <- depth_dt[match(sim_sys$wells, well), r_w]
  
  ## 4) Fit NNLS
  fit <- fit_cxcl13_nnls(
    A              = sim_sys$A,
    y_well         = sim_sys$y_obs,
    use_background = TRUE,
    use_weights    = TRUE,
    r_w            = r_w_vec,
    cx_true        = sim_sys$cx_true,
    verbose        = verbose
  )
  
  ## 5) Metrics: correlation + recall
  cx_true <- sim_sys$cx_true
  cx_hat  <- fit$cx_hat
  is_reactive <- sim_sys$is_reactive
  
  ## enforce same length for safety
  n_common <- min(length(cx_true), length(cx_hat))
  cx_true <- cx_true[seq_len(n_common)]
  cx_hat  <- cx_hat[seq_len(n_common)]
  is_reactive <- is_reactive[seq_len(n_common)]
  
  ## Pearson correlation
  cor_val <- suppressWarnings(cor(cx_true, cx_hat, use = "complete.obs"))
  
  ## Recall: rank by cx_hat, top N = (# true reactive)
  n_reactive <- sum(is_reactive)
  if (n_reactive > 0 && all(is.finite(cx_hat))) {
    ord <- order(cx_hat, decreasing = TRUE)
    pred_reactive <- rep(FALSE, length(cx_hat))
    pred_reactive[ord[seq_len(n_reactive)]] <- TRUE
    
    tp     <- sum(pred_reactive & is_reactive)
    recall <- tp / n_reactive
  } else {
    recall <- NA_real_
  }
  
  return(list(
    cor                   = cor_val,
    recall                = recall,
    per_row_noise_sd_log2 = per_row_noise_sd_log2,
    measurement_sd        = measurement_sd
  ))
}

## ============================================================
## 6. Grid of noise values + QC summary
## ============================================================

set.seed(123)  # for reproducible seeds in measurement / system

row_noise_grid   <- c(0, 0.1, 0.3, 0.5)   # Noise 1
meas_noise_grid  <- c(0, 0.5, 1, 2, 5)    # Noise 2
n_reps_per_combo <- 50                    # adjust as needed

results_list <- list()
idx <- 1L

for (rn in row_noise_grid) {
  for (mn in meas_noise_grid) {
    for (rep_i in seq_len(n_reps_per_combo)) {
      
      sim_res <- run_one_sim(
        PBMC_tirtl_long        = PBMC_tirtl_long2,
        subset_rows            = c("A","B","C","D"),
        per_row_noise_sd_log2  = rn,
        measurement_sd         = mn,
        reactive_fraction      = 0.01,
        verbose                = FALSE
      )
      
      results_list[[idx]] <- data.table(
        per_row_noise_sd_log2 = sim_res$per_row_noise_sd_log2,
        measurement_sd        = sim_res$measurement_sd,
        cor                   = sim_res$cor,
        recall                = sim_res$recall
      )
      idx <- idx + 1L
    }
  }
}

results_dt <- rbindlist(results_list)

qc_summary <- results_dt[
  ,
  .(
    mean_cor             = mean(cor,    na.rm = TRUE),
    mean_recall          = mean(recall, na.rm = TRUE),
    prop_recall_ge_0_8   = mean(recall >= 0.8, na.rm = TRUE),
    prop_recall_ge_0_9   = mean(recall >= 0.9, na.rm = TRUE),
    n_runs               = .N
  ),
  by = .(per_row_noise_sd_log2, measurement_sd)
]

qc_summary[order(per_row_noise_sd_log2, measurement_sd)]



######################################################################33


















































####################################################


###############################################################################
## SIMULATION EXPERIMENTS:
## Effect of changing expression distributions on identifiability
##
## Includes:
##   1) Varying mean separation between reactive / non-reactive populations
##   2) Varying SDs (spread) of both populations
##
## Requirements:
##   - run_one_sim() already defined
##   - PBMC_tirtl_long2 loaded
###############################################################################

library(data.table)

###############################################################################
## 1. EXPERIMENT 1 — Effect of MEAN SEPARATION
##    (How far apart reactive vs non-reactive mean expression levels are)
###############################################################################

set.seed(1234)

# Fixed SDs
nonreact_sd_fixed <- 0.4
react_sd_fixed    <- 0.7

# Fix nonreactive mean and vary reactive mean
nonreact_mean_fixed <- 0.3
react_mean_grid     <- seq(0.3, 5, length.out = 60)

# Turn off Noise 1 and Noise 2
per_row_noise_sd_log2_fixed <- 0
measurement_sd_fixed        <- 0

n_reps <- 50  # repeat each setting

expr_gap_res_list <- list()
idx <- 1L

for (rm in react_mean_grid) {
  for (rep_i in seq_len(n_reps)) {
    
    sim_res <- run_one_sim(
      PBMC_tirtl_long        = PBMC_tirtl_long2,
      subset_rows            = c("A","B","C","D"),
      per_row_noise_sd_log2  = per_row_noise_sd_log2_fixed,
      measurement_sd         = measurement_sd_fixed,
      reactive_fraction      = 0.01,
      nonreact_mean_log2     = nonreact_mean_fixed,
      nonreact_sd_log2       = nonreact_sd_fixed,
      react_mean_log2        = rm,
      react_sd_log2          = react_sd_fixed,
      verbose                = FALSE
    )
    
    expr_gap_res_list[[idx]] <- data.table(
      nonreact_mean_log2 = nonreact_mean_fixed,
      react_mean_log2    = rm,
      nonreact_sd_log2   = nonreact_sd_fixed,
      react_sd_log2      = react_sd_fixed,
      cor                = sim_res$cor,
      recall             = sim_res$recall
    )
    idx <- idx + 1L
  }
}

expr_gap_dt <- rbindlist(expr_gap_res_list)

# Summary over mean separation
expr_gap_summary <- expr_gap_dt[
  ,
  .(
    mean_gap           = mean(react_mean_log2 - nonreact_mean_log2),
    mean_cor           = mean(cor,    na.rm = TRUE),
    mean_recall        = mean(recall, na.rm = TRUE),
    prop_recall_ge_0_5 = mean(recall >= 0.5, na.rm = TRUE),
    prop_recall_ge_0_8 = mean(recall >= 0.8, na.rm = TRUE)
  ),
  by = react_mean_log2
][order(react_mean_log2)]

print(expr_gap_summary)

# Optional quick plot
plot(
  x    = expr_gap_summary$mean_gap,
  y    = expr_gap_summary$mean_recall,
  type = "b",
  xlab = "Mean gap (reactive_mean_log2 − nonreactive_mean_log2)",
  ylab = "Mean recall",
  main = "Effect of Mean Separation on Recall"
)
abline(h = c(0.5, 0.8), lty = 2)


###############################################################################
## 2. EXPERIMENT 2 — Effect of SD (spread)
##    (How wide or tight the expression distributions are)
###############################################################################

set.seed(5678)

# Fix means and vary SDs
nonreact_mean_fixed <- 0.3
react_mean_fixed    <- 3.5

nonreact_sd_grid <- c(0.2, 0.4, 0.8)
react_sd_grid    <- c(0.2, 0.5, 1.0)

# Noise off again
per_row_noise_sd_log2_fixed <- 0
measurement_sd_fixed        <- 0

expr_sd_res_list <- list()
idx <- 1L

for (nsd in nonreact_sd_grid) {
  for (rsd in react_sd_grid) {
    for (rep_i in seq_len(n_reps)) {
      
      sim_res <- run_one_sim(
        PBMC_tirtl_long        = PBMC_tirtl_long2,
        subset_rows            = c("A","B","C","D"),
        per_row_noise_sd_log2  = per_row_noise_sd_log2_fixed,
        measurement_sd         = measurement_sd_fixed,
        reactive_fraction      = 0.01,
        nonreact_mean_log2     = nonreact_mean_fixed,
        nonreact_sd_log2       = nsd,
        react_mean_log2        = react_mean_fixed,
        react_sd_log2          = rsd,
        verbose                = FALSE
      )
      
      expr_sd_res_list[[idx]] <- data.table(
        nonreact_sd_log2   = nsd,
        react_sd_log2      = rsd,
        nonreact_mean_log2 = nonreact_mean_fixed,
        react_mean_log2    = react_mean_fixed,
        cor                = sim_res$cor,
        recall             = sim_res$recall
      )
      idx <- idx + 1L
    }
  }
}

expr_sd_dt <- rbindlist(expr_sd_res_list)

# Summary of effect of spreads
expr_sd_summary <- expr_sd_dt[
  ,
  .(
    mean_cor           = mean(cor,    na.rm = TRUE),
    mean_recall        = mean(recall, na.rm = TRUE),
    prop_recall_ge_0_5 = mean(recall >= 0.5, na.rm = TRUE),
    prop_recall_ge_0_8 = mean(recall >= 0.8, na.rm = TRUE)
  ),
  by = .(nonreact_sd_log2, react_sd_log2)
][order(nonreact_sd_log2, react_sd_log2)]

print(expr_sd_summary)

# Optional heatmap-like visualization
recall_mat <- tapply(
  expr_sd_summary$mean_recall,
  INDEX = list(expr_sd_summary$nonreact_sd_log2,
               expr_sd_summary$react_sd_log2),
  FUN   = identity
)

image(
  x = as.numeric(rownames(recall_mat)),
  y = as.numeric(colnames(recall_mat)),
  z = recall_mat,
  xlab = "Non-reactive SD (log2)",
  ylab = "Reactive SD (log2)",
  main = "Effect of SDs on Mean Recall"
)
contour(
  x = as.numeric(rownames(recall_mat)),
  y = as.numeric(colnames(recall_mat)),
  z = recall_mat,
  add = TRUE
)

