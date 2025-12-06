

library(data.table)
library(ggplot2)

###############################
# 1) Keep only 80 clonotypes  #
###############################

# Helper: keep top N clonotypes by total rel_intensity
subset_top_n_clonotypes <- function(dt, n_keep = 80) {
  dt <- as.data.table(dt)
  
  clono_rank <- dt[
    , .(total_rel = sum(rel_intensity, na.rm = TRUE)),
    by = clonotype_id
  ]
  setorder(clono_rank, -total_rel)  # descending
  
  n_keep <- min(n_keep, nrow(clono_rank))
  keep_ids <- clono_rank$clonotype_id[seq_len(n_keep)]
  
  dt_kept <- dt[clonotype_id %in% keep_ids]
  
  if (nrow(dt_kept) == 0) {
    stop("subset_top_n_clonotypes(): no rows left after subsetting; check n_keep or dt.")
  }
  
  dt_kept
}

# Apply to your full PBMC_tirtl_long2
PBMC_80 <- subset_top_n_clonotypes(PBMC_tirtl_long2, n_keep = 80)

cat("After subsetting to top 80 clonotypes globally:\n")
cat("  Rows:", nrow(PBMC_80), "\n")
cat("  Unique clonotypes:", length(unique(PBMC_80$clonotype_id)), "\n")
cat("  Unique wells:", length(unique(PBMC_80$well)), "\n\n")

##########################################
# 2) Simulate CXCL13 on these 80 clones  #
##########################################

# You can adjust these parameters if you want
nonreact_mean_log2 <- 0.0
nonreact_sd_log2   <- 0.1
react_mean_log2    <- 3.5
react_sd_log2      <- 0.7
reactive_fraction  <- 0.05
per_row_noise_sd_log2 <- 0.3
measurement_sd     <- 0.6

# Simulate clonotype-level + row-level CXCL13 on this reduced set
sim_rows_80 <- simulate_cxcl13_on_wells(
  PBMC_tirtl_long      = PBMC_80,
  subset_rows          = c("A","B","C","D"),   # or NULL if you want all plate rows
  subset_wells         = NULL,
  nonreact_mean_log2   = nonreact_mean_log2,
  nonreact_sd_log2     = nonreact_sd_log2,
  react_mean_log2      = react_mean_log2,
  react_sd_log2        = react_sd_log2,
  reactive_fraction    = reactive_fraction,
  per_row_noise_sd_log2 = per_row_noise_sd_log2,
  seed_clono_noise     = 777,
  seed_row_noise       = 12,
  verbose              = FALSE
)

# Simulate well-level system (A, cx_true, y_true, y_obs)
sim_80 <- simulate_cxcl13_system(
  PBMC_tirtl_long = sim_rows_80,
  measurement_sd  = measurement_sd,
  verbose         = FALSE
)

cat("Reduced system (80 clonotypes):\n")
cat("  Wells:", nrow(sim_80$A), "\n")
cat("  Clonotypes (columns in A):", ncol(sim_80$A), "\n\n")

#############################################
# 3) Correlation vs number of wells (80 cl) #
#############################################

# Helper: run one fit on a subset of wells and return correlation
run_one_well_subset_80 <- function(sim,
                                   n_wells,
                                   seed = 1,
                                   verbose = FALSE) {
  set.seed(seed)
  
  total_wells <- nrow(sim$A)
  if (n_wells > total_wells) {
    stop("Requested n_wells (", n_wells, ") > total number of wells (", total_wells, ").")
  }
  
  idx <- sort(sample(seq_len(total_wells), size = n_wells, replace = FALSE))
  
  A_sub   <- sim$A[idx, , drop = FALSE]
  y_sub   <- sim$y_obs[idx]
  cx_true <- sim$cx_true
  
  # Fit with your current "best" NNLS configuration
  fit_sub <- fit_cxcl13_nnls(
    A              = A_sub,
    y_well         = y_sub,
    use_background = FALSE,
    use_weights    = FALSE,
    r_w            = NULL,
    cx_true        = cx_true,
    verbose        = verbose,
    filter_low_A   = FALSE,
    shrink_cx      = FALSE
  )
  
  cx_hat <- fit_sub$cx_hat
  
  ok <- is.finite(cx_true) & is.finite(cx_hat)
  if (sum(ok) < 2) {
    r_val <- NA_real_
  } else {
    r_val <- cor(cx_true[ok], cx_hat[ok])
  }
  
  data.table(
    n_wells = n_wells,
    r_val   = r_val
  )
}

# Grid of well counts
n_total_80   <- nrow(sim_80$A)
n_well_grid  <- seq(10, n_total_80, by = 10)  # adjust step if you want
n_reps       <- 5  # repeats per point

res_wells_80 <- rbindlist(lapply(n_well_grid, function(nw) {
  rbindlist(lapply(1:n_reps, function(rep_id) {
    run_one_well_subset_80(
      sim   = sim_80,
      n_wells = nw,
      seed    = 2000 + nw + rep_id,
      verbose = FALSE
    )[, rep := rep_id]
  }))
}))

summary_wells_80 <- res_wells_80[
  , .(
    mean_r = mean(r_val, na.rm = TRUE),
    sd_r   = sd(r_val,   na.rm = TRUE)
  ),
  by = n_wells
]

cat("Correlation vs number of wells (80-clonotype system):\n")
print(summary_wells_80)

########################################
# 4) Plot in nice academic style       #
########################################

ggplot(summary_wells_80, aes(x = n_wells, y = mean_r)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = mean_r - sd_r, ymax = mean_r + sd_r),
              alpha = 0.2, linetype = 0) +
  labs(
    x = "Number of wells used in NNLS fit",
    y = "Pearson correlation (cx_true vs cx_hat)",
    title = "Effect of well count on CXCL13 deconvolution (80 clonotypes)"
  ) +
  theme_classic(base_size = 14)
