library(data.table)
library(nnls)  

### Function Definitions ###

### simulate
simulate_cxcl13_system <- function(PBMC_tirtl_long,
                                   measurement_sd = 5,
                                   verbose = TRUE) {
  library(data.table)
  PBMC_tirtl_long <- as.data.table(PBMC_tirtl_long)
  
  if (!"CXCL13_simulated_exp_raw" %in% names(PBMC_tirtl_long)) {
    stop("Column 'CXCL13_simulated_exp_raw' not found in PBMC_tirtl_long. ",
         "Run your clonotype-level CXCL13 simulation first.")
  }
  
  if (verbose) {
    cat("=== simulate_cxcl13_system() ===\n")
    cat("Input rows:", nrow(PBMC_tirtl_long), "\n")
    cat("Number of unique wells:", length(unique(PBMC_tirtl_long$well)), "\n")
    cat("Number of unique clonotypes:", length(unique(PBMC_tirtl_long$clonotype_id)), "\n\n")
  }
  
  # 1) Design matrix A (wells x clonotypes)
  if (verbose) cat("Step 1: Building design matrix A (wells x clonotypes)...\n")
  
  A_dt <- dcast(
    PBMC_tirtl_long,
    well ~ clonotype_id,
    value.var = "rel_intensity",
    fun.aggregate = sum,
    fill = 0
  )
  
  wells <- A_dt$well
  A <- as.matrix(A_dt[, -"well"])
  colnames(A) <- paste0("clono", colnames(A))
  
  if (verbose) {
    cat("  Design matrix dimensions: ", nrow(A), " wells x ", ncol(A), " clonotypes\n", sep = "")
    cat("  First 5 wells:\n"); print(head(wells, 5))
    cat("  First 5 clonotype columns:\n"); print(head(colnames(A), 5))
    cat("\n")
  }
  
  # 2) Extract true CXCL13 + reactive flags in same order as columns of A
  if (verbose) cat("Step 2: Extracting true per-clonotype CXCL13 values...\n")
  
  # one row per clonotype
  cx_per_clone <- unique(
    PBMC_tirtl_long[, .(clonotype_id, CXCL13_simulated_exp_raw)]
  )
  
  # if you also stored a reactive flag in PBMC_tirtl_long, pick it up; otherwise approximate
  if ("is_reactive" %in% names(PBMC_tirtl_long)) {
    reactive_per_clone <- unique(
      PBMC_tirtl_long[, .(clonotype_id, is_reactive)]
    )
  } else {
    # fallback: everything "non-reactive" if not present
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
  
  # 3) y_true and y_obs
  #if (verbose) cat("Step 3: Computing true well-level CXCL13 (y_true = A %*% cx_true)...\n")
  #y_true <- as.numeric(A %*% cx_true)
  
  # 3) y_true and y_obs  (NOW USING ROW-LEVEL NOISY CXCL13)
  if (verbose) cat("Step 3: Computing true well-level CXCL13 using row-level noisy CXCL13_simulated_exp...\n")
  
  if (!"CXCL13_simulated_exp" %in% names(PBMC_tirtl_long)) {
    stop("Column 'CXCL13_simulated_exp' not found in PBMC_tirtl_long. ",
         "Did you run simulate_cxcl13_on_wells() with per-row noise?")
  }
  
  # For each clonotype–well row: contribution = rel_intensity * noisy CXCL13
  PBMC_tirtl_long[, contrib := rel_intensity * CXCL13_simulated_exp]
  
  # Sum contributions per well to get y_true
  y_true_dt <- PBMC_tirtl_long[
    , .(y_true = sum(contrib, na.rm = TRUE)),
    by = well
  ]
  
  # Match the ordering of 'wells' (from the design matrix A)
  y_true <- y_true_dt[match(wells, well), y_true]
  
  if (verbose) {
    cat("  Summary of y_true (with row-level noise included):\n")
    print(summary(y_true))
    cat("\n")
  }
  
  
  if (verbose) {
    cat("  Summary of y_true:\n")
    print(summary(y_true))
    cat("\n")
  }
  
  if (verbose) cat("Step 4: Adding measurement noise to simulate observed y_obs...\n")
  set.seed(2025)
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








#### fit 

fit_cxcl13_nnls <- function(A,
                            y_well,
                            use_background   = TRUE,
                            use_weights      = FALSE,
                            r_w              = NULL,       # total TIRTL reads per well
                            weight_clip      = c(1e-3, 1e3),
                            cx_true          = NULL,
                            verbose          = TRUE,
                            # NEW: optional column filtering
                            filter_low_A     = TRUE,
                            A_min_total      = NULL,
                            # NEW: optional shrinkage on cx_hat
                            shrink_cx        = TRUE,
                            shrink_lambda    = 0) {
  if (verbose) {
    cat("=== fit_cxcl13_nnls() ===\n")
    cat("Input matrix A: ", nrow(A), " wells x ", ncol(A), " clonotypes\n", sep = "")
    cat("Length of y_well:", length(y_well), "\n\n")
  }
  
  if (nrow(A) != length(y_well)) {
    stop("Number of rows in A (", nrow(A),
         ") does not match length of y_well (", length(y_well), ").")
  }
  
  # -------- 0) Optional filtering of low-contribution clonotypes --------
  orig_p   <- ncol(A)
  keep_idx <- seq_len(orig_p)
  dropped  <- integer(0)
  
  if (filter_low_A) {
    col_sums <- colSums(A, na.rm = TRUE)
    
    # If user didn't specify a threshold, use 10th percentile of non-zero columns
    if (is.null(A_min_total)) {
      pos <- col_sums[col_sums > 0]
      if (length(pos) == 0) {
        warning("All column sums are zero; skipping filter_low_A step.")
      } else {
        thr <- as.numeric(stats::quantile(pos, 0.1))
      }
    } else {
      thr <- A_min_total
    }
    
    if (exists("thr")) {
      keep_idx <- which(col_sums >= thr)
      dropped  <- setdiff(seq_len(orig_p), keep_idx)
      
      if (length(keep_idx) == 0) {
        warning("filter_low_A removed all columns; reverting to full A.")
        keep_idx <- seq_len(orig_p)
        dropped  <- integer(0)
      } else {
        A <- A[, keep_idx, drop = FALSE]
        if (!is.null(cx_true) && length(cx_true) == orig_p) {
          cx_true <- cx_true[keep_idx]
        }
        if (verbose) {
          cat("Step 0: filter_low_A enabled.\n")
          cat("  Threshold (A_min_total):", thr, "\n")
          cat("  Kept", length(keep_idx), "of", orig_p, "clonotypes.\n\n")
        }
      }
    }
  }
  
  # -------- 1) Compute weights w_w from r_w (or fallback) --------
  if (use_weights) {
    if (verbose) cat("Step 1: Computing row weights based on TIRTL depth r_w...\n")
    
    if (!is.null(r_w)) {
      if (length(r_w) != nrow(A)) {
        stop("Length of r_w (", length(r_w),
             ") does not match number of wells / rows in A (", nrow(A), ").")
      }
      depth <- as.numeric(r_w)
      if (any(depth < 0, na.rm = TRUE)) {
        stop("r_w contains negative values; depths must be non-negative.")
      }
    } else {
      if (verbose) {
        cat("  NOTE: r_w not provided. Falling back to depth <- rowSums(A).\n",
            "  Because A usually contains relative intensities, this may be uninformative.\n", sep = "")
      }
      depth <- rowSums(A, na.rm = TRUE)
    }
    
    # avoid zeros when computing ratio
    depth_pos <- depth[depth > 0]
    if (length(depth_pos) == 0) {
      stop("All depths are zero or NA; cannot compute weights.")
    }
    med_depth <- median(depth_pos)
    
    depth_ratio <- depth / med_depth
    w <- sqrt(depth_ratio)
    
    # enforce 10^-3 <= w_w <= 10^3
    w[!is.finite(w)] <- 1
    w <- pmax(weight_clip[1], pmin(weight_clip[2], w))
    
    if (verbose) {
      cat("  Summary of depth (r_w):\n");     print(summary(depth))
      cat("  Median depth:", med_depth, "\n")
      cat("  Summary of depth_ratio:\n");     print(summary(depth_ratio))
      cat("  Summary of final weights (w_w):\n"); print(summary(w))
      cat("\n")
    }
    
    A_fit <- A * w
    y_fit <- y_well * w
  } else {
    if (verbose) cat("Step 1: No weights used (unweighted NNLS).\n\n")
    A_fit <- A
    y_fit <- y_well
    w     <- rep(1, length(y_well))
  }
  
  # -------- 2) Add background column if requested --------
  if (use_background) {
    if (verbose) cat("Step 2: Adding background column (ones) to A...\n")
    A_fit_bg <- cbind(A_fit, bg = 1)
  } else {
    if (verbose) cat("Step 2: No background term included.\n")
    A_fit_bg <- A_fit
  }
  
  if (verbose) {
    cat("  Fitting matrix dimensions: ", nrow(A_fit_bg), " x ", ncol(A_fit_bg), "\n\n", sep = "")
    cat("Step 3: Running NNLS...\n")
  }
  
  # -------- 3) Run NNLS --------
  fit_obj <- nnls(A_fit_bg, y_fit)
  coefs   <- coef(fit_obj)
  
  if (use_background) {
    cx_hat <- coefs[1:(ncol(A_fit_bg) - 1)]
    bg_hat <- coefs[ncol(A_fit_bg)]
  } else {
    cx_hat <- coefs
    bg_hat <- NULL
  }
  
  # save raw (unshrunken) estimates
  cx_hat_raw <- cx_hat
  
  # -------- 4) Fitted values and residuals (on original scale) --------
  if (verbose) cat("Step 4: Computing fitted values and residuals...\n")
  
  y_hat_weighted <- as.numeric(A_fit_bg %*% coefs)
  y_hat          <- y_hat_weighted / w      # back to unweighted scale
  residuals      <- y_well - y_hat
  
  # -------- 4b) Optional shrinkage of cx_hat --------
  if (shrink_cx && shrink_lambda > 0) {
    if (verbose) {
      cat("Step 4b: Applying shrinkage to cx_hat (lambda =", shrink_lambda, ")...\n")
    }
    cx_hat <- pmax(0, cx_hat_raw - shrink_lambda)
  }
  
  if (verbose) {
    cat("  Summary of estimated CXCL13 per clonotype (cx_hat):\n")
    print(summary(cx_hat))
    if (!is.null(bg_hat)) {
      cat("  Estimated background (bg_hat):", bg_hat, "\n")
    }
    cat("  Summary of fitted y_hat:\n")
    print(summary(y_hat))
    cat("  Summary of residuals (y_well - y_hat):\n")
    print(summary(residuals))
    cat("\n")
  }
  
  # -------- 5) Optional correlation with ground truth --------
  if (!is.null(cx_true)) {
    if (length(cx_true) == length(cx_hat_raw)) {
      cor_val <- cor(cx_true, cx_hat_raw, use = "complete.obs")
      if (verbose) {
        cat("Step 5: Pearson correlation between true and estimated CXCL13 (raw cx_hat):",
            round(cor_val, 3), "\n")
      }
    } else {
      warning("Length of cx_true (", length(cx_true),
              ") does not match length of cx_hat (", length(cx_hat_raw), "). Skipping correlation.")
    }
  }
  
  if (verbose) cat("=== Fitting complete. Returning estimates and diagnostics. ===\n\n")
  
  return(list(
    cx_hat      = cx_hat,      # possibly shrunk
    cx_hat_raw  = cx_hat_raw,  # original NNLS estimate
    bg_hat      = bg_hat,
    y_hat       = y_hat,
    residuals   = residuals,
    weights     = w,
    fit_obj     = fit_obj,
    kept_idx    = keep_idx,
    dropped_idx = dropped
  ))
}







###### QC plot #####
plot_cxcl13_true_vs_est <- function(cx_true,
cx_hat,
is_reactive = NULL,
rel_tol    = 0.20,
main       = "Comparison of Ground Truth and Estimated CXCL13 Levels per Clonotype") {
  
  # Keep only finite values
  ok <- is.finite(cx_true) & is.finite(cx_hat)
  x  <- cx_true[ok]
  y  <- cx_hat[ok]
  
  if (!is.null(is_reactive)) {
    is_reactive <- is_reactive[ok]
  } else {
    is_reactive <- rep(FALSE, length(x))
  }
  
  # Relative error
  denom    <- pmax(abs(x), 1e-6)
  rel_err  <- abs(y - x) / denom
  within_tol <- rel_err <= rel_tol
  
  # Colors for 4 categories
  cols <- ifelse(is_reactive & within_tol, "red3",
                 ifelse(is_reactive & !within_tol, "lightcoral",
                        ifelse(!is_reactive & within_tol, "black", "grey50")))
  
  # Save old par
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  # Increase right margin even further to prevent clipping
  par(mar = c(5, 5, 4, 15))   # bottom, left, top, RIGHT (bigger!)
  
  # Scatterplot
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
  
  # Pearson r
  r_val <- cor(x, y)
  subtitle_text <- sprintf(
    "Pearson correlation (r) = %.3f   (n = %d)",
    r_val, length(x)
  )
  title(sub = subtitle_text, cex.sub = 0.8)
  
  # Legend in right margin
  par(xpd = NA)
  
  usr <- par("usr")
  
  # Increase horizontal offset for legend placement
  legend_x <- usr[2] + (usr[2] - usr[1]) * 0.12   # 12% of plot width to right
  legend_y <- usr[4]
  
  legend(
    x = legend_x,
    y = legend_y,
    legend = c(
      #paste0("Reactive, within ", rel_tol * 100, "%"),
      #paste0("Reactive, outside ", rel_tol * 100, "%"),
      #paste0("Non-reactive, within ", rel_tol * 100, "%"),
      #paste0("Non-reactive, outside ", rel_tol * 100, "%")
      
      paste0("Reactive, correctly estimated"),
      paste0("Reactive"),
      paste0("Non-reactive, correctly estimated"),
      paste0("Non-reactive")
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





####################################################
################ function calls ######################
##################################################


# runs te simulation
sim <- simulate_cxcl13_system(PBMC_tirtl_long_sim2_rows_AD, measurement_sd = 0.1, verbose = TRUE)

# total TIRTL reads per well from the same subset used to build A
depth_dt <- PBMC_tirtl_long_sim2_rows_AD[
  , .(r_w = sum(intensity, na.rm = TRUE)),
  by = well]

# match ordering to rows of A / y_well
r_w_vec <- depth_dt[match(sim$wells, well), r_w]


# Then fit
fit <- fit_cxcl13_nnls(
  A             = sim$A,
  y_well        = sim$y_obs,
  use_background = TRUE,
  use_weights    = TRUE,
  r_w            = r_w_vec,
  cx_true        = sim$cx_true,
  verbose        = TRUE
)


plot_res <- plot_cxcl13_true_vs_est(
  cx_true     = sim$cx_true,
  cx_hat      = fit$cx_hat,
  is_reactive = sim$is_reactive,
  rel_tol     = 0.80
)

fit2 <- fit_cxcl13_nnls(
  A              = sim$A,
  y_well         = sim$y_obs,
  use_background = FALSE,
  use_weights    = FALSE,
  r_w            = r_w_vec,
  cx_true        = sim$cx_true,
  verbose        = TRUE,
  filter_low_A   = FALSE,
  shrink_cx      = FALSE,
  shrink_lambda  = 0.5
)


plot_res <- plot_cxcl13_true_vs_est(
  cx_true     = sim$cx_true,
  cx_hat      = fit2$cx_hat,
  is_reactive = sim$is_reactive,
  rel_tol     = 0.90
)


