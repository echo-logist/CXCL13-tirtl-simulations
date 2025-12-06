library(data.table)

# dt: PBMC_tirtl_long-style table (one row per clonotype–well)
# top_frac: fraction of clonotypes to keep (e.g. 0.1 = top 10%)
subset_top_clonotypes <- function(dt, top_frac = 0.2) {
  dt <- as.data.table(dt)
  
  # Total abundance per clonotype (global across wells)
  clono_rank <- dt[, .(total_rel = sum(rel_intensity, na.rm = TRUE)),
                   by = clonotype_id]
  setorder(clono_rank, -total_rel)  # descending
  
  n_clono <- nrow(clono_rank)
  n_keep  <- max(1, floor(top_frac * n_clono))
  
  keep_ids <- clono_rank$clonotype_id[seq_len(n_keep)]
  
  dt_kept <- dt[clonotype_id %in% keep_ids]
  
  dt_kept
}

run_one_sim_topfrac <- function(PBMC_base,
                                top_frac       = 0.2,
                                reactive_q     = 0.9,
                                measurement_sd = 0.6,
                                verbose        = FALSE,
                                seed_global    = 123) {
  set.seed(seed_global)
  
  # 1) Keep only top X% clonotypes globally by abundance
  dt_top <- subset_top_clonotypes(PBMC_base, top_frac = top_frac)
  
  # 2) Simulate CXCL13 per row + clonotype
  sim_rows <- simulate_cxcl13_on_wells(
    PBMC_tirtl_long = dt_top,
    subset_rows     = c("A","B","C","D"),  # or NULL if you want all rows
    subset_wells    = NULL,
    verbose         = verbose
  )
  
  # 3) Simulate well-level CXCL13 system with measurement noise
  sim <- simulate_cxcl13_system(
    PBMC_tirtl_long = sim_rows,
    measurement_sd  = measurement_sd,
    verbose         = verbose
  )
  
  # 4) (Optional) depth per well, not used for fit here but kept for completeness
  depth_dt <- sim_rows[
    , .(r_w = sum(intensity, na.rm = TRUE)),
    by = well
  ]
  r_w_vec <- depth_dt[match(sim$wells, well), r_w]
  
  # 5) Fit NNLS with your preferred configuration
  fit <- fit_cxcl13_nnls(
    A              = sim$A,
    y_well         = sim$y_obs,
    use_background = FALSE,
    use_weights    = FALSE,
    r_w            = r_w_vec,
    cx_true        = sim$cx_true,
    verbose        = verbose,
    filter_low_A   = FALSE,
    shrink_cx      = FALSE  # since this is currently best for you
  )
  
  cx_true     <- sim$cx_true
  cx_hat      <- fit$cx_hat
  is_reactive <- sim$is_reactive
  
  ok <- is.finite(cx_true) & is.finite(cx_hat)
  cx_true     <- cx_true[ok]
  cx_hat      <- cx_hat[ok]
  is_reactive <- is_reactive[ok]
  
  # Pearson correlation
  r_val <- cor(cx_true, cx_hat)
  
  # Recall of reactive clonotypes at a fixed quantile of cx_hat
  thr <- quantile(cx_hat, probs = reactive_q, na.rm = TRUE)
  pred_reactive   <- cx_hat >= thr
  recall_reactive <- mean(pred_reactive[is_reactive])
  
  # Rough measure of clonotypes per well after thinning
  mean_clonos_per_well <- mean(rowSums(sim$A > 0))
  
  data.table(
    top_frac             = top_frac,
    mean_clonos_per_well = mean_clonos_per_well,
    r_val                = r_val,
    recall_reactive      = recall_reactive
  )
}



top_fracs <- c(0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0)
n_reps    <- 10

results_topfrac <- rbindlist(lapply(top_fracs, function(tf) {
  rbindlist(lapply(1:n_reps, function(rep_id) {
    run_one_sim_topfrac(
      PBMC_base    = PBMC_tirtl_long2,
      top_frac     = tf,
      measurement_sd = 0.6,
      reactive_q   = 0.9,
      seed_global  = 1000 + rep_id,
      verbose      = FALSE
    )[, rep := rep_id]
  }))
}))

summary_topfrac <- results_topfrac[
  , .(
    mean_clonos_per_well = mean(mean_clonos_per_well),
    sd_clonos_per_well   = sd(mean_clonos_per_well),
    mean_r               = mean(r_val),
    sd_r                 = sd(r_val),
    mean_recall          = mean(recall_reactive),
    sd_recall            = sd(recall_reactive)
  ),
  by = top_frac
]

print(summary_topfrac)
