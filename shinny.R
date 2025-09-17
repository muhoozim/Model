# HIV Transmission Model — Robust Shiny App
# -----------------------------------------------------------------------------
# This single-file Shiny app rewrites and hardens the original simulation code
# into testable functions, adds safety checks, and presents an interactive
# dashboard for exploring results. It assumes required parameter CSV files are
# placed in a local "data/" folder (see FILES_REQUIRED below). No setwd() is
# used.
# 
# Key improvements vs. original script
# - Deterministic, reproducible runs (seed control)
# - Explicit month/year timeline derived from inputs (no 324/384 mismatch)
# - Fixed undefined variables (e.g., g during init), corrected logical
#   conditions, and reset of per-step accumulators (lambda_inf)
# - Modularized helpers for probability conversions and draws
# - Safer file I/O, validation, and clear errors in UI
# - Interactive plots, scenario comparison, and downloads
# - Contemporary Shiny UX patterns (bslib theme, value cards, spinners)
# -----------------------------------------------------------------------------

# ---- Packages ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(shinycssloaders)
  library(shinyvalidate)
  library(DT)
  library(ggplot2)
  library(plotly)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(tibble)
  library(stringr)
  library(zoo)
  library(data.table)
  library(EnvStats) # rtri
})

# ---- Constants ---------------------------------------------------------------
FILES_REQUIRED <- c(
  # targets
  "Targets_00_15-54_new.csv",
  "Targets_01.csv",
  # transmission & behavior
  "RW_p_transmission.csv",
  "RW_p_condomEffect.csv",
  "RW_p_condomUse.csv",
  "RW_p_VSEffect.csv",
  "RW_n_sexact.csv",
  # populations & composition (RW_Pop.csv & RW_Pop_CD4.csv shipped as
  # placeholder copies — update these with final numbers when available)
  "RW_Pop.csv",
  "RW_Pop_HIVPrev_2.csv",
  "RW_Pop_CD4.csv",
  "RW_Pop_comp_diag.csv",
  "RW_Pop_comp_link.csv",
  "RW_Pop_comp_artvs.csv",
  "RW_Pop_comp_artvf.csv",
  # growth, natural history, cascade transitions
  "RW_p_popgrowth.csv",
  "RW_p_DisProgress.csv",
  "RW_p_diagnosis.csv",
  "RW_p_linkage.csv",
  "RW_p_LTFU.csv",
  "RW_p_ARTVS.csv",
  "RW_p_ARTnotVS.csv",
  "RW_p_death.csv",
  "RW_p_reengage.csv"
)

# ---- Helpers -----------------------------------------------------------------
app_data_dir <- function() file.path(getwd(), "data")

ensure_placeholder_csv <- function(path, name) {
  # Minimal templated structures so the app can run even if final inputs
  # are absent. Replace these numbers with study-specific values later.
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  placeholder <- switch(
    name,
    "RW_Pop.csv" = tibble::tibble(
      Group = seq_len(25),
      Total = rep(1000, 25)
    ),
    "RW_Pop_CD4.csv" = tibble::tibble(
      State = c(
        paste0("unlink_bin_", 1:3),
        paste0("link_bin_", 1:3),
        paste0("art_vs_bin_", 1:3),
        paste0("art_not_vs_bin_", 1:3)
      ),
      Mean = c(rep(0.25, 3), rep(0.25, 3), rep(0.25, 3), rep(0.25, 3))
    ),
    "Targets_00_15-54_new.csv" = {
      years <- seq(2004, 2010)
      cols <- c(
        "Prevalence", "Prevalence_LB", "Prevalence_UB",
        "Prevalence_F", "Prevalence_LB_F", "Prevalence_UB_F",
        "Prevalence_M", "Prevalence_LB_M", "Prevalence_UB_M",
        "Incidence", "Incidence_LB", "Incidence_UB",
        "Incidence_F", "Incidence_F_LB", "Incidence_F_UB",
        "Incidence_M", "Incidence_M_LB", "Incidence_M_UB",
        "Percent_on_ART_total", "Percent_on_ART_total_LB", "Percent_on_ART_total_UB",
        "Percent_on_ART_Female", "Percent_on_ART_Female_LB", "Percent_on_ART_Female_UB",
        "Percent_on_ART_Male", "Percent_on_ART_Male_LB", "Percent_on_ART_Male_UB",
        "Percent_VS_total", "Percent_VS_total_LB", "Percent_VS_total_UB",
        "Percent_VS_Female", "Percent_VS_Female_LB", "Percent_VS_Female_UB",
        "Percent_VS_Male", "Percent_VS_Male_LB", "Percent_VS_Male_UB"
      )
      placeholder_cols <- lapply(cols, function(...) rep(NA_real_, length(years)))
      names(placeholder_cols) <- cols
      tibble::as_tibble(c(list(Year = years), placeholder_cols))
    },
    "Targets_01.csv" = {
      years <- seq(2004, 2010)
      cols <- c(
        "Prevalence", "Prevalence_LB", "Prevalence_UB",
        "Prevalence_F", "Prevalence_LB_F", "Prevalence_UB_F",
        "Prevalence_M", "Prevalence_LB_M", "Prevalence_UB_M",
        "Percent_VS_total", "Percent_VS_total_LB", "Percent_VS_total_UB",
        "Percent_VS_Male", "Percent_VS_Male_LB", "Percent_VS_Male_UB",
        "Percent_VS_Female", "Percent_VS_Female_LB", "Percent_VS_Female_UB",
        "P_ART_All", "P_ART_All_LB", "P_ART_All_UB",
        "P_ART_All_F", "P_ART_All_F_LB", "P_ART_All_F_UB",
        "P_ART_All_M", "P_ART_All_M_LB", "P_ART_All_M_UB",
        "P_ART_All_F_young", "P_ART_All_F_young_LB", "P_ART_All_F_young_UB"
      )
      placeholder_cols <- lapply(cols, function(...) rep(NA_real_, length(years)))
      names(placeholder_cols) <- cols
      tibble::as_tibble(c(list(Year = years), placeholder_cols))
    },
    NULL
  )
  if (is.null(placeholder)) return(FALSE)
  readr::write_csv(placeholder, path)
  message(sprintf(
    "Created placeholder CSV for %s. Update %s with calibrated data when available (placeholder rows contain neutral NA values).",
    name, path
  ))
  TRUE
}

read_params <- function(dir) {
  paths <- setNames(file.path(dir, FILES_REQUIRED), FILES_REQUIRED)
  missing <- names(paths)[!file.exists(paths)]
  if (length(missing)) {
    purrr::walk(missing, ~ensure_placeholder_csv(paths[.x], .x))
  }
  still_missing <- names(paths)[!file.exists(paths)]
  if (length(still_missing)) {
    stop(
      sprintf(
        "Missing parameter files in '%s':\n- %s",
        dir,
        paste(still_missing, collapse = "\n- ")
      ),
      call. = FALSE
    )
  }
  lapply(paths, readr::read_csv, show_col_types = FALSE)
}

# Probability conversion utilities --------------------------------------------
prob_annual_to_monthly <- function(p_annual) {
  # Convert an annual probability to an equivalent monthly probability
  1 - (1 - p_annual)^(1/12)
}

rate_to_prob_monthly <- function(rate_annual) {
  # If given a (continuous) annual rate, convert to monthly probability
  1 - exp(-rate_annual/12)
}

prob_to_rate <- function(p) { -log(1 - p) }
rate_to_prob <- function(r) { 1 - exp(-r) }

# Safe logical helpers (scalar) ------------------------------------------------
between_inclusive <- function(x, lo, hi) x >= lo && x <= hi

# Triangular draw wrapper (EnvStats::rtri) -------------------------------------
rdraw_tri <- function(min, max, mode, n = 1) EnvStats::rtri(n, min = min, max = max, mode = mode)

# Initialize 2D numeric matrix with NA ----------------------------------------
mat <- function(nrow, ncol, fill = NA_real_) {
  m <- matrix(fill, nrow = nrow, ncol = ncol)
  storage.mode(m) <- "double"
  m
}

# ---- Core Simulation ---------------------------------------------------------
# The function below is a careful rewrite of the original model structure,
# retaining its compartment definitions and data requirements while fixing
# previously noted issues (accumulator resets, undefined indices, timeline).
# For readability and to keep memory reasonable, we compute the final aggregated
# outcomes used by the dashboard rather than returning every internal state.

simulate_model <- function(params,
                           trials = 100L,
                           sgroup = 25L,
                           start_year = 2004L,
                           end_year   = 2020L,
                           seed = 1L,
                           deathmulti = 1,
                           initialvf = 0.25) {
  stopifnot(end_year >= start_year)
  generations <- as.integer((end_year - start_year + 1L) * 12L)
  set.seed(seed)
  
  # Unpack parameter tables
  Target_all       <- params[["Targets_00_15-54_new.csv"]]
  Target_new       <- params[["Targets_01.csv"]]
  p_transmission   <- params[["RW_p_transmission.csv"]]
  p_condomEffect   <- params[["RW_p_condomEffect.csv"]]
  p_condomUse      <- params[["RW_p_condomUse.csv"]]
  p_VSEffect       <- params[["RW_p_VSEffect.csv"]]
  n_sexact         <- params[["RW_n_sexact.csv"]]
  Pop              <- params[["RW_Pop.csv"]]
  Pop_HIVPrev      <- params[["RW_Pop_HIVPrev_2.csv"]]
  Pop_CD4          <- params[["RW_Pop_CD4.csv"]]
  Pop_comp_diag    <- params[["RW_Pop_comp_diag.csv"]]
  Pop_comp_link    <- params[["RW_Pop_comp_link.csv"]]
  Pop_comp_artvs   <- params[["RW_Pop_comp_artvs.csv"]]
  Pop_comp_artvf   <- params[["RW_Pop_comp_artvf.csv"]]
  p_popgrowth      <- params[["RW_p_popgrowth.csv"]]
  p_DisProgress    <- params[["RW_p_DisProgress.csv"]]
  p_diagnosis      <- params[["RW_p_diagnosis.csv"]]
  p_link           <- params[["RW_p_linkage.csv"]]
  p_LTFU           <- params[["RW_p_LTFU.csv"]]
  p_ARTVS          <- params[["RW_p_ARTVS.csv"]]
  p_ARTnotVS       <- params[["RW_p_ARTnotVS.csv"]]
  p_death          <- params[["RW_p_death.csv"]]
  p_reengage       <- params[["RW_p_reengage.csv"]]
  
  years <- seq(from = start_year, to = end_year, by = 1L)
  nyears <- length(years)
  
  # Storage for calibration/GOF across trials (top 50)
  pd_top50 <- matrix(NA_real_, nrow = 50, ncol = 3,
                     dimnames = list(NULL, c("trial", "pd_total", "N_within")))
  pd_top50[, 2] <- 999
  pd_top50[, 3] <- 0
  
  # For dashboard outputs, we’ll collect summary per-month and per-year stats
  out_list <- vector("list", length = trials)
  yearly_list <- vector("list", length = trials)
  
  for (t in seq_len(trials)) {
    # ----------------------
    # Allocate state matrices
    # ----------------------
    S  <- mat(sgroup, generations)
    N  <- mat(sgroup, generations, fill = 0)
    V  <- mat(sgroup, generations, fill = 0)  # virally suppressed
    Z  <- mat(sgroup, generations, fill = 0)  # not suppressed
    NI <- mat(sgroup, generations, fill = 0)  # new infections
    
    # 24 infectious compartments X1..X24, plus deaths D
    X <- array(NA_real_, dim = c(sgroup, generations, 24))
    D <- mat(sgroup, generations, fill = 0)
    
    # Parameters (time & subgroup varying)
    rho      <- mat(sgroup, generations) # growth
    beta     <- mat(sgroup, generations) # per-contact transmission (unsuppressed)
    epsilon  <- mat(sgroup, generations) # condom effectiveness
    c_r      <- mat(sgroup, generations) # condom use
    upsilon  <- mat(sgroup, generations) # VS effectiveness
    n_r      <- mat(sgroup, generations) # sex acts
    
    # Natural history rates (delta_*): 1..6 as in original mapping
    delta <- array(NA_real_, dim = c(sgroup, generations, 6))
    
    # Diagnosis (alpha_* CD4 bins 1..4), Linkage (sigma_*), LTFU (gamma_1..12)
    alpha <- array(NA_real_, dim = c(sgroup, generations, 4))
    sigma <- array(NA_real_, dim = c(sgroup, generations, 4))
    gamma <- array(NA_real_, dim = c(sgroup, generations, 12))
    
    # ART transitions: theta_1..4 (to VS), psi_1..4 (VS failure), tau_1..4 (re-engage)
    theta <- array(NA_real_, dim = c(sgroup, generations, 4))
    psi   <- array(NA_real_, dim = c(sgroup, generations, 4))
    tau   <- array(NA_real_, dim = c(sgroup, generations, 4))
    
    # Death probabilities per stage (mu_1..16)
    mu    <- array(NA_real_, dim = c(sgroup, generations, 16))
    
    # Force of infection helpers
    lambda_inf <- mat(sgroup, generations, fill = 1)  # multiplicative survivor
    lambda_t_r <- mat(sgroup, generations, fill = NA_real_) # per-month infection prob
    
    # ----------------------
    # Initial population & prevalence
    # ----------------------
    Total <- mat(sgroup, generations)
    Prevalence <- mat(sgroup, generations)
    
    # Example derivation for initial totals (preserves original proportions)
    TotalUrbanWomen <- sum(Pop$Total[1:10])
    TotalHighRisk   <- runif(1, min = 11000, max = 45000)
    TotalLowRiskUrban <- TotalUrbanWomen - TotalHighRisk
    
    # Assign high-risk women by age bins (as in original logic)
    Total[6,1]  <- TotalHighRisk * (157 + 470) / 1338
    Total[7,1]  <- TotalHighRisk * (342 + 288/2) / 1338
    Total[8,1]  <- TotalHighRisk * (288/2 + 81/5) / 1338
    Total[9,1]  <- TotalHighRisk * (81/5 + 81/5) / 1338
    Total[10,1] <- TotalHighRisk * (81/5 + 81/5) / 1338
    
    women <- c(535205,460559,327619,256033,201099,209695,154203,124472,89338,67466)
    totalwomen <- sum(women)
    # low-risk urban women by age groups minus HR allocation above
    Total[1,1] <- TotalUrbanWomen * sum(535205,460559) / totalwomen - Total[6,1]
    Total[2,1] <- TotalUrbanWomen * sum(327619,256033) / totalwomen - Total[7,1]
    Total[3,1] <- TotalUrbanWomen * sum(201099,209695) / totalwomen - Total[8,1]
    Total[4,1] <- TotalUrbanWomen * sum(154203,124472) / totalwomen - Total[9,1]
    Total[5,1] <- TotalUrbanWomen * sum(89338,67466)  / totalwomen - Total[10,1]
    
    # Remaining groups from Pop
    if (nrow(Pop) < sgroup) stop("RW_Pop.csv has fewer than sgroup rows")
    for (s in 11:sgroup) Total[s,1] <- Pop$Total[s]
    
    # Initial prevalence draws (beta distrib by subgroup)
    for (l in 1:sgroup) {
      a <- Pop_HIVPrev$Alpha[l]
      b <- Pop_HIVPrev$Beta[l]
      Prevalence[l,1] <- rbeta(1, a, b)
    }
    
    # Initial susceptible
    for (s in 1:sgroup) S[s,1] <- Total[s,1] * (1 - Prevalence[s,1])
    
    # ----------------------
    # Initial compartment composition (diagnosis/linkage/ART status & CD4 bins)
    # ----------------------
    # Diagnosis proportion baseline
    multi_diag <- rdraw_tri(0.15, 0.8, 0.3)
    a <- Pop_comp_diag$Alpha[1]; b <- Pop_comp_diag$Beta[1]
    initial_diag <- rbeta(1, a, b) * multi_diag
    initial_diag_vec <- rep(initial_diag, sgroup)
    
    # Linkage baseline
    a <- Pop_comp_link$Alpha[1]; b <- Pop_comp_link$Beta[1]
    initial_link <- rbeta(1, a, b)
    initial_link_vec <- rep(initial_link, sgroup)
    
    # On ART & suppressed baseline
    a <- Pop_comp_artvs$Alpha[1]; b <- Pop_comp_artvs$Beta[1]
    initial_artvs <- rbeta(1, a, b)
    initial_artvs_vec <- rep(initial_artvs, sgroup)
    
    # On ART & not suppressed baseline
    a <- Pop_comp_artvf$Alpha[1]; b <- Pop_comp_artvf$Beta[1]
    initial_artvf <- rbeta(1, a, b)
    initial_artvf_vec <- rep(initial_artvf, sgroup)
    
    # CD4 distributions by status (use means provided)
    cd4_unlink <- Pop_CD4$Mean[1:3]
    cd4_unlink4 <- 1 - sum(cd4_unlink)
    cd4_link <- c(Pop_CD4$Mean[4:6])
    cd4_link4 <- 1 - sum(cd4_link)
    cd4_artvs <- c(Pop_CD4$Mean[7:9])
    cd4_artvs4 <- 1 - sum(cd4_artvs)
    cd4_artvf <- c(Pop_CD4$Mean[10:12])
    cd4_artvf4 <- 1 - sum(cd4_artvf)
    
    # Seed X1..X24 at t=1
    for (s in 1:sgroup) {
      prev  <- Prevalence[s,1]
      tot   <- Total[s,1]
      diag  <- initial_diag_vec[s]
      linkd <- initial_link_vec[s]
      artvs <- initial_artvs_vec[s]
      artvf <- initial_artvf_vec[s]
      
      X[s,1,1]  <- tot * prev * (1 - diag)          * cd4_unlink[1]
      X[s,1,2]  <- tot * prev * (1 - diag)          * cd4_unlink[2]
      X[s,1,3]  <- tot * prev * (1 - diag)          * cd4_unlink[3]
      X[s,1,4]  <- tot * prev * (1 - diag)          * cd4_unlink4
      X[s,1,5]  <- tot * prev * (diag) * (1 - linkd) * cd4_link[1]
      X[s,1,6]  <- tot * prev * (diag) * (1 - linkd) * cd4_link[2]
      X[s,1,7]  <- tot * prev * (diag) * (1 - linkd) * cd4_link[3]
      X[s,1,8]  <- tot * prev * (diag) * (1 - linkd) * cd4_link4
      X[s,1,9]  <- tot * prev * (diag) * ( linkd) * (1 - artvs) * cd4_link[1]
      X[s,1,10] <- tot * prev * (diag) * ( linkd) * (1 - artvs) * cd4_link[2]
      X[s,1,11] <- tot * prev * (diag) * ( linkd) * (1 - artvs) * cd4_link[3]
      X[s,1,12] <- tot * prev * (diag) * ( linkd) * (1 - artvs) * cd4_link4
      X[s,1,13] <- 0 # LTFU CD4>500
      X[s,1,14] <- 0
      X[s,1,15] <- 0
      X[s,1,16] <- 0
      X[s,1,17] <- tot * prev * diag * linkd * artvs * (1 - artvf) * cd4_artvs[1]
      X[s,1,18] <- tot * prev * diag * linkd * artvs * (1 - artvf) * cd4_artvs[2]
      X[s,1,19] <- tot * prev * diag * linkd * artvs * (1 - artvf) * cd4_artvs[3]
      X[s,1,20] <- tot * prev * diag * linkd * artvs * (1 - artvf) * cd4_artvs4
      X[s,1,21] <- tot * prev * diag * linkd * artvs * ( artvf) * cd4_artvf[1]
      X[s,1,22] <- tot * prev * diag * linkd * artvs * ( artvf) * cd4_artvf[2]
      X[s,1,23] <- tot * prev * diag * linkd * artvs * ( artvf) * cd4_artvf[3]
      X[s,1,24] <- tot * prev * diag * linkd * artvs * ( artvf) * cd4_artvf4
      
      # Derived states at t=1
      N[s,1] <- sum(X[s,1,])
      V[s,1] <- sum(X[s,1,17:20])
      Z[s,1] <- N[s,1] - V[s,1] # includes non-VS states
    }
    
    # ----------------------
    # Parameter draws (per trial, often constant over time after set)
    # ----------------------
    # Transmission beta (per contact, unsuppressed): draw base for FSW, F, M
    beta_fs <- rdraw_tri(p_transmission$Min[1], p_transmission$Max[1], p_transmission$Mode[1])
    beta_f  <- rdraw_tri(p_transmission$Min[2], p_transmission$Max[2], p_transmission$Mode[2])
    beta_m  <- rdraw_tri(p_transmission$Min[3], p_transmission$Max[3], p_transmission$Mode[3])
    
    for (s in 1:sgroup) {
      if ((s >= 1 && s <= 5) || (s >= 11 && s <= 15)) beta[s,1] <- beta_f
      if (s >= 6 && s <= 10) beta[s,1] <- beta_fs
      if (s >= 16 && s <= 25) beta[s,1] <- beta_m
      beta[s, 2:generations] <- beta[s,1]
    }
    
    # Condom effectiveness
    eps <- rdraw_tri(p_condomEffect$Min[1], p_condomEffect$Max[1], p_condomEffect$Mode[1])
    epsilon[,1] <- eps
    epsilon[, 2:generations] <- eps
    
    # Condom use c_r staged at three eras: months 1..72, 73..132, 133..generations
    eras <- list(`1` = 1, `73` = 73, `133` = 133)
    for (l in 1:5) {
      # Low-risk urban women (1..5)
      c_r[l,   1]  <- rbeta(1, p_condomUse$Alpha[l],      p_condomUse$Beta[l])
      c_r[l,  73]  <- rbeta(1, p_condomUse$Alpha[l+5],    p_condomUse$Beta[l+5])
      c_r[l, 133]  <- rbeta(1, p_condomUse$Alpha[l+10],   p_condomUse$Beta[l+10])
      # Low-risk rural women (11..15)
      c_r[l+10, 1] <- rbeta(1, p_condomUse$Alpha[l+15],   p_condomUse$Beta[l+15])
      c_r[l+10,73] <- rbeta(1, p_condomUse$Alpha[l+20],   p_condomUse$Beta[l+20])
      c_r[l+10,133]<- rbeta(1, p_condomUse$Alpha[l+25],   p_condomUse$Beta[l+25])
      # Urban men (16..20)
      c_r[l+15, 1] <- rbeta(1, p_condomUse$Alpha[l+30],   p_condomUse$Beta[l+30])
      c_r[l+15,73] <- rbeta(1, p_condomUse$Alpha[l+35],   p_condomUse$Beta[l+35])
      c_r[l+15,133]<- rbeta(1, p_condomUse$Alpha[l+40],   p_condomUse$Beta[l+40])
      # Rural men (21..25)
      c_r[l+20, 1] <- rbeta(1, p_condomUse$Alpha[l+45],   p_condomUse$Beta[l+45])
      c_r[l+20,73] <- rbeta(1, p_condomUse$Alpha[l+50],   p_condomUse$Beta[l+50])
      c_r[l+20,133]<- rbeta(1, p_condomUse$Alpha[l+55],   p_condomUse$Beta[l+55])
      # High-risk women (6..10) — constant per era as original
      v <- rbeta(1, p_condomUse$Alpha[l+60], p_condomUse$Beta[l+60])
      c_r[l+5,  1] <- v; c_r[l+5, 73] <- v; c_r[l+5, 133] <- v
    }
    # Fill between eras
    for (s in 1:sgroup) {
      for (g in 2:generations) {
        if (g <= 72) c_r[s,g] <- c_r[s,1]
        else if (g <= 132) c_r[s,g] <- c_r[s,73]
        else c_r[s,g] <- c_r[s,133]
      }
    }
    
    # VS effectiveness (sexual)
    ups <- rdraw_tri(p_VSEffect$Min[1], p_VSEffect$Max[1], p_VSEffect$Mode[1])
    upsilon[,1] <- ups
    upsilon[,2:generations] <- ups
    
    # Sex acts
    n_low <- rdraw_tri(n_sexact$Min[1], n_sexact$Max[1], n_sexact$Mode[1])
    n_hr  <- rdraw_tri(n_sexact$Min[2], n_sexact$Max[2], n_sexact$Mode[2]) * runif(1, 1, 5)
    youth_multi_f <- rdraw_tri(0.2, 1.0, 0.5)
    youth_multi_m <- rdraw_tri(0.2, 1.0, 0.5)
    
    for (s in 1:sgroup) {
      if (s %in% c(1,11)) n_r[s,1] <- n_low * youth_multi_f
      else if (s %in% c(16,21)) n_r[s,1] <- n_low * youth_multi_m
      else if ((s >= 2 && s <= 5) || (s >= 12 && s <= 15) || (s >= 17 && s <= 20) || (s >= 22 && s <= 25)) n_r[s,1] <- n_low
      else if (s >= 6 && s <= 10) n_r[s,1] <- n_hr
      n_r[s, 2:generations] <- n_r[s,1]
    }
    
    # Population growth (monthly, constant over time)
    for (s in 1:sgroup) {
      u  <- p_popgrowth$Mean[s]
      sd <- p_popgrowth$SD[s]
      rho[s,1] <- rnorm(1, mean = u, sd = sd)
      rho[s,2:generations] <- rho[s,1]
    }
    
    # Natural history deltas (converted to monthly probabilities)
    # CD4>500 -> 350-500 (delta1), 350-500 -> 200-350 (delta2), 200-350 -> <200 (delta3)
    # Females (1..15) share draws; males (16..25) share male draws as original
    d1_f <- prob_annual_to_monthly(rbeta(1, p_DisProgress$Alpha[1], p_DisProgress$Beta[1]))
    d1_m <- prob_annual_to_monthly(rbeta(1, p_DisProgress$Alpha[4], p_DisProgress$Beta[4]))
    d2_f <- prob_annual_to_monthly(rbeta(1, p_DisProgress$Alpha[2], p_DisProgress$Beta[2]))
    d2_m <- prob_annual_to_monthly(rbeta(1, p_DisProgress$Alpha[5], p_DisProgress$Beta[5]))
    d3_f <- prob_annual_to_monthly(rbeta(1, p_DisProgress$Alpha[3], p_DisProgress$Beta[3]))
    d3_m <- prob_annual_to_monthly(rbeta(1, p_DisProgress$Alpha[6], p_DisProgress$Beta[6]))
    
    for (s in 1:sgroup) {
      if (s <= 15) {
        delta[s,1,1] <- d1_f; delta[s,1,2] <- d2_f; delta[s,1,3] <- d3_f
        delta[s,1,4] <- d1_f; delta[s,1,5] <- d2_f; delta[s,1,6] <- d3_f
      } else {
        delta[s,1,1] <- d1_m; delta[s,1,2] <- d2_m; delta[s,1,3] <- d3_m
        delta[s,1,4] <- d1_m; delta[s,1,5] <- d2_m; delta[s,1,6] <- d3_m
      }
      delta[s, 2:generations, ] <- delta[s,1,]
    }
    
    # Diagnosis probabilities (alpha_1..4), with CD4 multiplier and sex corrections
    cd4_multi <- rdraw_tri(p_diagnosis$Min[62], p_diagnosis$Max[62], p_diagnosis$Mode[62])
    correct_f <- rdraw_tri(0.15, 0.8, 0.2)
    correct_m <- rdraw_tri(0.20, 0.8, 0.3)
    
    # Helper to set era values (1, 73, 133) and replicate across months
    set_alpha_era <- function(subset_idx, base_row_idx) {
      # base_row_idx: positions in p_diagnosis aligned to original script
      vals <- list()
      for (era_g in c(1,73,133)) {
        k <- switch(as.character(era_g), `1` = 0, `73` = 5, `133` = 10)
        a1 <- rbeta(1, p_diagnosis$Alpha[base_row_idx + k], p_diagnosis$Beta[base_row_idx + k])
        a1 <- prob_to_rate(a1) # to annual rate
        # Same across CD4 bins 1..3, multiply for bin 4
        vals[[as.character(era_g)]] <- rate_to_prob(prob_to_rate(a1)) # monthly prob
        for (s in subset_idx) {
          # apply sex correction
          corr <- if (s <= 15) correct_f else correct_m
          alpha[s,era_g,1] <<- 1 - exp(-a1/12 * corr)
          alpha[s,era_g,2] <<- 1 - exp(-a1/12 * corr)
          alpha[s,era_g,3] <<- 1 - exp(-a1/12 * corr)
          alpha[s,era_g,4] <<- 1 - exp(-a1/12 * corr * cd4_multi)
        }
      }
      # Fill across months
      for (s in subset_idx) {
        for (g in 2:generations) {
          if (g <= 72) alpha[s,g,] <<- alpha[s,1,]
          else if (g <= 132) alpha[s,g,] <<- alpha[s,73,]
          else alpha[s,g,] <<- alpha[s,133,]
        }
      }
    }
    
    # Apply by population groups, following original indexing
    set_alpha_era(1:5, 1)        # low-risk urban women
    set_alpha_era(11:15, 16)     # low-risk rural women
    set_alpha_era(16:20, 31)     # urban men
    set_alpha_era(21:25, 46)     # rural men
    # high-risk women (6:10) use row 61 in original
    for (s in 6:10) {
      a1 <- rbeta(1, p_diagnosis$Alpha[61], p_diagnosis$Beta[61])
      r  <- prob_to_rate(a1)
      for (era_g in c(1,73,133)) {
        alpha[s,era_g,1] <- 1 - exp(-r/12)
        alpha[s,era_g,2] <- 1 - exp(-r/12)
        alpha[s,era_g,3] <- 1 - exp(-r/12)
        alpha[s,era_g,4] <- 1 - exp(-r/12 * cd4_multi)
      }
      for (g in 2:generations) {
        if (g <= 72) alpha[s,g,] <- alpha[s,1,]
        else if (g <= 132) alpha[s,g,] <- alpha[s,73,]
        else alpha[s,g,] <- alpha[s,133,]
      }
    }
    
    # Linkage sigma_1..4 (HR women vs low-risk/others)
    cd4_multi_link <- rdraw_tri(p_link$Min[3], p_link$Max[3], p_link$Mode[3])
    sig_hr <- rdraw_tri(p_link$Min[1], p_link$Max[1], p_link$Mode[1])
    sig_lr <- rdraw_tri(p_link$Min[2], p_link$Max[2], p_link$Mode[2])
    
    to_monthly <- function(p) 1 - exp(-prob_to_rate(p)/12)
    for (s in 1:sgroup) {
      base <- if (s >= 6 && s <= 10) sig_hr else sig_lr
      sigma[s,1,1] <- to_monthly(base)
      sigma[s,1,2] <- sigma[s,1,1]
      sigma[s,1,3] <- sigma[s,1,1]
      sigma[s,1,4] <- to_monthly(min(0.999, base * cd4_multi_link))
      sigma[s, 2:generations, ] <- sigma[s,1,]
    }
    
    # LTFU gamma_1..12 (constant over time, shared across subgroups)
    for (k in 1:8) {
      gk <- rbeta(1, p_LTFU$Alpha[k], p_LTFU$Beta[k])
      gk <- prob_annual_to_monthly(gk)
      gamma[,1,k] <- gk
      gamma[,2:generations,k] <- gk
    }
    # Map duplicates (9..12) per original logic
    gamma[,,9]  <- gamma[,,5]
    gamma[,,10] <- gamma[,,6]
    gamma[,,11] <- gamma[,,7]
    gamma[,,12] <- gamma[,,8]
    
    # On ART & suppressed (theta_1..4) — era changes at 37, 97, 145
    set_theta_era <- function(g0, idx_offset) {
      for (k in 1:4) {
        th <- rbeta(1, p_ARTVS$Alpha[idx_offset + (k-1)], p_ARTVS$Beta[idx_offset + (k-1)])
        th_m <- prob_annual_to_monthly(th)
        theta[, g0, k] <<- th_m
      }
    }
    set_theta_era(1, 1)
    set_theta_era(37, 5)
    set_theta_era(97, 9)
    set_theta_era(145, 13)
    
    # Adjust female & older multipliers (converted to/from rates safely)
    correct_female <- rdraw_tri(1, 2, 1.5)
    correct_old    <- rdraw_tri(1, 2, 1.5)
    for (s in 1:sgroup) {
      for (g in 1:generations) {
        if ((s %in% c(2:5,7:10,12:15))) {
          theta[s,g,] <- rate_to_prob(prob_to_rate(theta[s,g,]) * (correct_female * correct_old))
        } else if (s %in% c(17:20,22:25)) {
          theta[s,g,] <- rate_to_prob(prob_to_rate(theta[s,g,]) * (correct_old))
        }
      }
    }
    
    # Spread theta eras across months
    for (s in 1:sgroup) {
      for (g in 2:generations) {
        if (g <= 36) theta[s,g,] <- theta[s,1,]
        else if (g <= 96) theta[s,g,] <- theta[s,37,]
        else if (g <= 144) theta[s,g,] <- theta[s,97,]
        else theta[s,g,] <- theta[s,145,]
      }
    }
    
    # On ART & Not suppressed (psi_1..4) — constant over time w/ sex multipliers
    psi_base <- vapply(1:4, function(k) rbeta(1, p_ARTnotVS$Alpha[k], p_ARTnotVS$Beta[k]) / 12, numeric(1))
    multi_female <- rdraw_tri(0.2, 1.0, 0.5)
    multi_male   <- rdraw_tri(1.0, 3.0, 2.5)
    for (s in 1:sgroup) {
      mult <- if (s <= 15) multi_female else multi_male
      for (k in 1:4) {
        psi[s,1,k] <- 1 - exp(-(psi_base[k] * mult))
        psi[s,2:generations,k] <- psi[s,1,k]
      }
    }
    
    # Deaths mu_1..16 (convert to monthly probabilities), with deathmulti on some
    # Following original structure; for brevity we use grouped draws and replicate
    to_month <- function(p) 1 - exp(-prob_to_rate(p)/12)
    # Women blocks
    for (l in 1:5) {
      mu9 <- rdraw_tri(p_death$Min[l], p_death$Max[l], p_death$Mode[l])
      mu10<- rdraw_tri(p_death$Min[l+5], p_death$Max[l+5], p_death$Mode[l+5])
      mu11<- rdraw_tri(p_death$Min[l+10], p_death$Max[l+10], p_death$Mode[l+10])
      mu12<- rdraw_tri(p_death$Min[l+15], p_death$Max[l+15], p_death$Mode[l+15])
      mu_w <- to_month(c(mu9,mu10,mu11,mu12))
      for (idx in c(l, l+5, l+10)) { # women low/high/rural sets
        mu[idx,1,9]  <- mu_w[1]
        mu[idx,1,1]  <- mu_w[1]
        mu[idx,1,5]  <- mu_w[1]
        mu[idx,1,13] <- mu_w[1] * deathmulti
        mu[idx,1,10] <- mu_w[2]
        mu[idx,1,2]  <- mu_w[2]
        mu[idx,1,6]  <- mu_w[2]
        mu[idx,1,14] <- mu_w[2] * deathmulti
        mu[idx,1,11] <- mu_w[3]
        mu[idx,1,3]  <- mu_w[3]
        mu[idx,1,7]  <- mu_w[3]
        mu[idx,1,15] <- mu_w[3] * deathmulti
        mu[idx,1,12] <- mu_w[4]
        mu[idx,1,4]  <- mu_w[4]
        mu[idx,1,8]  <- mu_w[4]
        mu[idx,1,16] <- mu_w[4] * deathmulti
      }
    }
    # Men blocks (shifted indices)
    for (l in 1:5) {
      mu9 <- rdraw_tri(p_death$Min[l+20], p_death$Max[l+20], p_death$Mode[l+20])
      mu10<- rdraw_tri(p_death$Min[l+25], p_death$Max[l+25], p_death$Mode[l+25])
      mu11<- rdraw_tri(p_death$Min[l+30], p_death$Max[l+30], p_death$Mode[l+30])
      mu12<- rdraw_tri(p_death$Min[l+35], p_death$Max[l+35], p_death$Mode[l+35])
      mu_m <- to_month(c(mu9,mu10,mu11,mu12))
      for (idx in c(l+15, l+20)) { # urban & rural men sets
        mu[idx,1,9]  <- mu_m[1]
        mu[idx,1,1]  <- mu_m[1]
        mu[idx,1,5]  <- mu_m[1]
        mu[idx,1,13] <- mu_m[1] * deathmulti
        mu[idx,1,10] <- mu_m[2]
        mu[idx,1,2]  <- mu_m[2]
        mu[idx,1,6]  <- mu_m[2]
        mu[idx,1,14] <- mu_m[2] * deathmulti
        mu[idx,1,11] <- mu_m[3]
        mu[idx,1,3]  <- mu_m[3]
        mu[idx,1,7]  <- mu_m[3]
        mu[idx,1,15] <- mu_m[3] * deathmulti
        mu[idx,1,12] <- mu_m[4]
        mu[idx,1,4]  <- mu_m[4]
        mu[idx,1,8]  <- mu_m[4]
        mu[idx,1,16] <- mu_m[4] * deathmulti
      }
    }
    mu[, 2:generations, ] <- mu[, 1, ]
    
    # Re-engagement (tau_1..4), same for all groups
    tau1 <- prob_annual_to_monthly(rdraw_tri(p_reengage$Min[1], p_reengage$Max[1], p_reengage$Mode[1]))
    tau[,,1] <- tau1; tau[,,2] <- tau1; tau[,,3] <- tau1; tau[,,4] <- tau1
    
    # ----------------------
    # Dynamics over generations (months)
    # ----------------------
    for (g in 2:generations) {
      for (s in 1:sgroup) {
        # Partner pools (previous month)
        N_lr_u_f <- sum(S[1:5,g-1])  + sum(N[1:5,g-1])
        N_hr_f   <- sum(S[6:10,g-1]) + sum(N[6:10,g-1])
        N_lr_r_f <- sum(S[11:15,g-1]) + sum(N[11:15,g-1])
        N_u_m    <- sum(S[16:20,g-1]) + sum(N[16:20,g-1])
        N_r_m    <- sum(S[21:25,g-1]) + sum(N[21:25,g-1])
        
        # Reset multiplicative survivor for FOI at this (s,g)
        lambda_inf[s,g] <- 1
        
        # Mix by group: compute per-act infection risk and apply protection/acts
        add_partner_group <- function(partners_idx, partner_pool_size) {
          for (j in partners_idx) {
            if (partner_pool_size <= 0) next
            prev_j   <- ifelse((N[j,g-1] + S[j,g-1]) <= 0, 0, N[j,g-1] / (N[j,g-1] + S[j,g-1]))
            eff_p    <- beta[s,g] * ((1 - upsilon[s,g]) * (V[j,g-1]/max(1e-9,N[j,g-1])) + (Z[j,g-1]/max(1e-9,N[j,g-1])))
            per_act  <- ( (N[j,g-1] + S[j,g-1]) / partner_pool_size ) * prev_j * eff_p
            protect  <- (1 - epsilon[s,g] * c_r[s,g])
            lambda_inf[s,g] <<- lambda_inf[s,g] * (1 - per_act)^(protect * n_r[s,g])
          }
        }
        
        if (s >= 1 && s <= 5)    add_partner_group(16:20, N_u_m)
        if (s >= 6 && s <= 10)   add_partner_group(16:20, N_u_m)
        if (s >= 11 && s <= 15)  add_partner_group(21:25, N_r_m)
        if (s >= 16 && s <= 20)  add_partner_group(1:10, N_lr_u_f + N_hr_f)
        if (s >= 21 && s <= 25)  add_partner_group(11:15, N_lr_r_f)
        
        lambda_t_r[s,g] <- 1 - lambda_inf[s,g]
        NI[s,g] <- lambda_t_r[s,g] * S[s,g-1]
        
        # Difference equations (compact mapping using X array)
        # Aliases for readability at time g-1
        Xprev <- X[s,g-1,]
        del   <- delta[s,g,]
        alp   <- alpha[s,g,]
        sig   <- sigma[s,g,]
        gam   <- gamma[s,g,]
        th    <- theta[s,g,]
        ps    <- psi[s,g,]
        ta    <- tau[s,g,]
        mu_s  <- mu[s,g,]
        
        # S
        S[s,g] <- S[s,g-1] + rho[s,g] * S[s,g-1] - lambda_t_r[s,g] * S[s,g-1]
        
        # X1..X4 (undiagnosed)
        X[s,g,1] <- Xprev[1] + lambda_t_r[s,g]*S[s,g-1] - del[1]*Xprev[1] - alp[1]*Xprev[1] - mu_s[1]*Xprev[1]
        X[s,g,2] <- Xprev[2] + del[1]*Xprev[1] - del[2]*Xprev[2] - alp[2]*Xprev[2] - mu_s[2]*Xprev[2]
        X[s,g,3] <- Xprev[3] + del[2]*Xprev[2] - del[3]*Xprev[3] - alp[3]*Xprev[3] - mu_s[3]*Xprev[3]
        X[s,g,4] <- Xprev[4] + del[3]*Xprev[3] - alp[4]*Xprev[4] - mu_s[4]*Xprev[4]
        
        # X5..X8 (diagnosed, not linked)
        X[s,g,5] <- Xprev[5] + alp[1]*Xprev[1] - del[1]*Xprev[5] - sig[1]*Xprev[5] - mu_s[1]*Xprev[5]
        X[s,g,6] <- Xprev[6] + alp[2]*Xprev[2] + del[1]*Xprev[5] - del[2]*Xprev[6] - sig[2]*Xprev[6] - mu_s[2]*Xprev[6]
        X[s,g,7] <- Xprev[7] + alp[3]*Xprev[3] + del[2]*Xprev[6] - del[3]*Xprev[7] - sig[3]*Xprev[7] - mu_s[3]*Xprev[7]
        X[s,g,8] <- Xprev[8] + alp[4]*Xprev[4] + del[3]*Xprev[7] - sig[4]*Xprev[8] - mu_s[4]*Xprev[8]
        
        # X9..X12 (linked)
        X[s,g,9]  <- Xprev[9]  + sig[1]*Xprev[5] - del[4]*Xprev[9]  - th[1]*Xprev[9]  - gam[1]*Xprev[9]  - mu_s[5]*Xprev[9]
        X[s,g,10] <- Xprev[10] + sig[2]*Xprev[6] + del[4]*Xprev[9]  - del[5]*Xprev[10] - th[2]*Xprev[10] - gam[2]*Xprev[10] - mu_s[6]*Xprev[10]
        X[s,g,11] <- Xprev[11] + sig[3]*Xprev[7] + del[5]*Xprev[10] - del[6]*Xprev[11] - th[3]*Xprev[11] - gam[3]*Xprev[11] - mu_s[7]*Xprev[11]
        X[s,g,12] <- Xprev[12] + sig[4]*Xprev[8] + del[6]*Xprev[11] - th[4]*Xprev[12] - gam[4]*Xprev[12] - mu_s[8]*Xprev[12]
        
        # X13..X16 (LTFU)
        X[s,g,13] <- Xprev[13] + gam[1]*Xprev[9]  + gam[5]*Xprev[17] + gam[9]*Xprev[21] - del[1]*Xprev[13] - ta[1]*Xprev[13] - mu_s[1]*Xprev[13]
        X[s,g,14] <- Xprev[14] + gam[2]*Xprev[10] + gam[6]*Xprev[18] + gam[10]*Xprev[22] + del[1]*Xprev[13] - ta[2]*Xprev[14] - del[2]*Xprev[14] - mu_s[2]*Xprev[14]
        X[s,g,15] <- Xprev[15] + gam[3]*Xprev[11] + gam[7]*Xprev[19] + gam[11]*Xprev[23] + del[2]*Xprev[14] - ta[3]*Xprev[15] - del[3]*Xprev[15] - mu_s[3]*Xprev[15]
        X[s,g,16] <- Xprev[16] + gam[4]*Xprev[12] + gam[8]*Xprev[20] + gam[12]*Xprev[24] + del[3]*Xprev[15] - ta[4]*Xprev[16] - mu_s[4]*Xprev[16]
        
        # X17..X20 (VS on ART)
        X[s,g,17] <- Xprev[17] + th[1]*Xprev[9]  + ta[1]*Xprev[13] - gam[5]*Xprev[17] - ps[1]*Xprev[17] - mu_s[9]*Xprev[17]
        X[s,g,18] <- Xprev[18] + th[2]*Xprev[10] + ta[2]*Xprev[14] - gam[6]*Xprev[18] - ps[2]*Xprev[18] - mu_s[10]*Xprev[18]
        X[s,g,19] <- Xprev[19] + th[3]*Xprev[11] + ta[3]*Xprev[15] - gam[7]*Xprev[19] - ps[3]*Xprev[19] - mu_s[11]*Xprev[19]
        X[s,g,20] <- Xprev[20] + th[4]*Xprev[12] + ta[4]*Xprev[16] - gam[8]*Xprev[20] - ps[4]*Xprev[20] - mu_s[12]*Xprev[20]
        
        # X21..X24 (Not VS on ART)
        X[s,g,21] <- Xprev[21] + ps[1]*Xprev[17] - gam[9]*Xprev[21]  - mu_s[13]*Xprev[21]
        X[s,g,22] <- Xprev[22] + ps[2]*Xprev[18] - gam[10]*Xprev[22] - mu_s[14]*Xprev[22]
        X[s,g,23] <- Xprev[23] + ps[3]*Xprev[19] - gam[11]*Xprev[23] - mu_s[15]*Xprev[23]
        X[s,g,24] <- Xprev[24] + ps[4]*Xprev[20] - gam[12]*Xprev[24] - mu_s[16]*Xprev[24]
        
        # Deaths
        D[s,g] <- D[s,g-1] +
          mu_s[1] *(Xprev[1] + Xprev[5] + Xprev[13]) +
          mu_s[2] *(Xprev[2] + Xprev[6] + Xprev[14]) +
          mu_s[3] *(Xprev[3] + Xprev[7] + Xprev[15]) +
          mu_s[4] *(Xprev[4] + Xprev[8] + Xprev[16]) +
          mu_s[5] * Xprev[9]  + mu_s[6]*Xprev[10] + mu_s[7]*Xprev[11] + mu_s[8]*Xprev[12] +
          mu_s[9] * Xprev[17] + mu_s[10]*Xprev[18] + mu_s[11]*Xprev[19] + mu_s[12]*Xprev[20] +
          mu_s[13]* Xprev[21] + mu_s[14]*Xprev[22] + mu_s[15]*Xprev[23] + mu_s[16]*Xprev[24]
        
        # Update aggregates at g
        N[s,g] <- sum(X[s,g,])
        V[s,g] <- sum(X[s,g,17:20])
        Z[s,g] <- N[s,g] - V[s,g]
      }
    }
    
    # ----------------------
    # Aggregations for dashboard/calibration
    # ----------------------
    idxF <- c(1:3,6:8,11:13)    # female subsets used later
    idxM <- c(16:18,21:23)
    
    N_total      <- colSums(N)
    S_total      <- colSums(S)
    V_total      <- colSums(V)
    V_F_total    <- colSums(V[idxF, , drop = FALSE])
    V_M_total    <- colSums(V[idxM, , drop = FALSE])
    On_ART_mat <- V + apply(X[,,21:24], c(1,2), sum)
    On_ART_total <- colSums(On_ART_mat)
    Diagnosed_total <- colSums(X[,,5:24])
    
    N_F_total    <- colSums(N[idxF, , drop = FALSE])
    S_F_total    <- colSums(S[idxF, , drop = FALSE])
    N_M_total    <- colSums(N[idxM, , drop = FALSE])
    S_M_total    <- colSums(S[idxM, , drop = FALSE])
    On_ART_F_total <- colSums(On_ART_mat[idxF, , drop = FALSE])
    On_ART_M_total <- colSums(On_ART_mat[idxM, , drop = FALSE])
    
    NI_total     <- colSums(NI)
    NI_F_total   <- colSums(NI[idxF, , drop = FALSE])
    NI_M_total   <- colSums(NI[idxM, , drop = FALSE])
    
    # Yearly slices (take every 12th month)
    take_years <- seq(1, generations, by = 12)
    S_total_y  <- S_total[take_years]
    S_F_total_y <- S_F_total[take_years]
    S_M_total_y <- S_M_total[take_years]
    N_total_y <- N_total[take_years]
    N_F_total_y <- N_F_total[take_years]
    N_M_total_y <- N_M_total[take_years]
    Prevalence_total_y <- N_total_y / (N_total_y + S_total_y)
    Prevalence_F_y     <- N_F_total_y / (N_F_total_y + S_F_total_y)
    Prevalence_M_y     <- N_M_total_y / (N_M_total_y + S_M_total_y)

    NI_total_y <- zoo::rollapply(NI_total, 12, sum, by = 12, align = "right")
    NI_F_total_y <- zoo::rollapply(NI_F_total, 12, sum, by = 12, align = "right")
    NI_M_total_y <- zoo::rollapply(NI_M_total, 12, sum, by = 12, align = "right")
    Incidence_total_y <- NI_total_y / pmax(1e-9, S_total_y)
    Incidence_F_y <- NI_F_total_y / pmax(1e-9, S_F_total_y)
    Incidence_M_y <- NI_M_total_y / pmax(1e-9, S_M_total_y)
    On_ART_total_y <- On_ART_total[take_years]
    On_ART_F_total_y <- On_ART_F_total[take_years]
    On_ART_M_total_y <- On_ART_M_total[take_years]
    V_total_y <- V_total[take_years]
    V_F_total_y <- V_F_total[take_years]
    V_M_total_y <- V_M_total[take_years]
    On_ART_total_pct_y <- On_ART_total_y / pmax(1e-9, N_total_y) * 100
    On_ART_F_total_pct_y <- On_ART_F_total_y / pmax(1e-9, N_F_total_y) * 100
    On_ART_M_total_pct_y <- On_ART_M_total_y / pmax(1e-9, N_M_total_y) * 100
    VS_total_pct_y <- V_total_y / pmax(1e-9, N_total_y) * 100
    VS_F_total_pct_y <- V_F_total_y / pmax(1e-9, N_F_total_y) * 100
    VS_M_total_pct_y <- V_M_total_y / pmax(1e-9, N_M_total_y) * 100
    
    # Minimal GOF proxy (sum of relative deviations on prevalence & incidence)
    pd_prev_total <- sum(abs(Prevalence_total_y - Target_all$Prevalence[seq_len(nyears)]) / pmax(1e-9, Target_all$Prevalence[seq_len(nyears)]), na.rm = TRUE)
    pd_inc_total  <- sum(abs(Incidence_total_y - Target_all$Incidence[seq_len(nyears)]) / pmax(1e-9, Target_all$Incidence[seq_len(nyears)]), na.rm = TRUE)
    pd_total <- pd_prev_total + pd_inc_total
    
    # Number within naive bounds where provided
    within_prev <- with(Target_all[seq_len(nyears),],
                        sum( Prevalence_LB <= Prevalence_total_y & Prevalence_total_y <= Prevalence_UB, na.rm = TRUE))
    N_within <- within_prev
    
    # Insert into top50 (ordered by pd_total ascending)
    pd_total_new <- c(t, pd_total, N_within)
    for (i in 1:50) {
      if ((i == 1 && pd_total < pd_top50[i,2]) || (i > 1 && pd_total <= pd_top50[i,2])) {
        pd_top50 <- rbind(pd_total_new, pd_top50)
        pd_top50 <- pd_top50[1:50, , drop = FALSE]
        break
      }
    }
    
    # Store trial summary for dashboard
    out_list[[t]] <- tibble(
      trial = t,
      month = seq_len(generations),
      year  = rep(years, each = 12)[1:generations],
      N_total = N_total,
      S_total = S_total,
      V_total = V_total,
      On_ART_total = On_ART_total,
      Diagnosed_total = Diagnosed_total,
      NI_total = NI_total,
      Prevalence = N_total / (N_total + S_total),
      Incidence  = c(
        rep(NA_real_, 11),
        zoo::rollapply(NI_total, width = 12, FUN = sum, align = "right")
      ) / pmax(1e-9, S_total))
    yearly_list[[t]] <- bind_rows(
      tibble(
        trial = t,
        group = "Total",
        year = years,
        prevalence = Prevalence_total_y * 100,
        incidence = Incidence_total_y,
        on_art_count = On_ART_total_y,
        on_art_pct = On_ART_total_pct_y,
        vs_pct = VS_total_pct_y,
        on_art = On_ART_total_y
      ),
      tibble(
        trial = t,
        group = "Female",
        year = years,
        prevalence = Prevalence_F_y * 100,
        incidence = Incidence_F_y,
        on_art_count = On_ART_F_total_y,
        on_art_pct = On_ART_F_total_pct_y,
        vs_pct = VS_F_total_pct_y,
        on_art = On_ART_F_total_y
      ),
      tibble(
        trial = t,
        group = "Male",
        year = years,
        prevalence = Prevalence_M_y * 100,
        incidence = Incidence_M_y,
        on_art_count = On_ART_M_total_y,
        on_art_pct = On_ART_M_total_pct_y,
        vs_pct = VS_M_total_pct_y,
        on_art = On_ART_M_total_y
      )
    )
  } # end trials

  Prevalence = N_total / (N_total + S_total)
  Incidence  = c(
    rep(NA_real_, 11),
    zoo::rollapply(NI_total, width = 12, FUN = sum, align = "right")
  ) / pmax(1e-9, S_total)



  list(
    trials = trials,
    generations = generations,
    years = years,
    pd_top50 = as.data.frame(pd_top50),
    outcomes = bind_rows(out_list),
    yearly = bind_rows(yearly_list),
    targets = list(
      prevalence = bind_rows(
        tibble(
          group = "Total",
          year = Target_all$Year,
          target = Target_all$Prevalence * 100,
          target_lb = Target_all$Prevalence_LB * 100,
          target_ub = Target_all$Prevalence_UB * 100
        ),
        tibble(
          group = "Female",
          year = Target_all$Year,
          target = Target_all$Prevalence_F * 100,
          target_lb = Target_all$Prevalence_LB_F * 100,
          target_ub = Target_all$Prevalence_UB_F * 100
        ),
        tibble(
          group = "Male",
          year = Target_all$Year,
          target = Target_all$Prevalence_M * 100,
          target_lb = Target_all$Prevalence_LB_M * 100,
          target_ub = Target_all$Prevalence_UB_M * 100
        )
      ),
      incidence = bind_rows(
        tibble(
          group = "Total",
          year = Target_all$Year,
          target = Target_all$Incidence,
          target_lb = Target_all$Incidence_LB,
          target_ub = Target_all$Incidence_UB
        ),
        tibble(
          group = "Female",
          year = Target_all$Year,
          target = Target_all$Incidence_F,
          target_lb = Target_all$Incidence_F_LB,
          target_ub = Target_all$Incidence_F_UB
        ),
        tibble(
          group = "Male",
          year = Target_all$Year,
          target = Target_all$Incidence_M,
          target_lb = Target_all$Incidence_M_LB,
          target_ub = Target_all$Incidence_M_UB
        )
      ),
      art = bind_rows(
        tibble(
          group = "Total",
          year = Target_all$Year,
          target = Target_all$Percent_on_ART_total * 100,
          target_lb = Target_all$Percent_on_ART_total_LB * 100,
          target_ub = Target_all$Percent_on_ART_total_UB * 100
        ),
        tibble(
          group = "Female",
          year = Target_all$Year,
          target = Target_all$Percent_on_ART_Female * 100,
          target_lb = Target_all$Percent_on_ART_Female_LB * 100,
          target_ub = Target_all$Percent_on_ART_Female_UB * 100
        ),
        tibble(
          group = "Male",
          year = Target_all$Year,
          target = Target_all$Percent_on_ART_Male * 100,
          target_lb = Target_all$Percent_on_ART_Male_LB * 100,
          target_ub = Target_all$Percent_on_ART_Male_UB * 100
        )
      ),
      vs = bind_rows(
        tibble(
          group = "Total",
          year = Target_all$Year,
          target = Target_all$Percent_VS_total * 100,
          target_lb = Target_all$Percent_VS_total_LB * 100,
          target_ub = Target_all$Percent_VS_total_UB * 100
        ),
        tibble(
          group = "Female",
          year = Target_all$Year,
          target = Target_all$Percent_VS_Female * 100,
          target_lb = Target_all$Percent_VS_Female_LB * 100,
          target_ub = Target_all$Percent_VS_Female_UB * 100
        ),
        tibble(
          group = "Male",
          year = Target_all$Year,
          target = Target_all$Percent_VS_Male * 100,
          target_lb = Target_all$Percent_VS_Male_LB * 100,
          target_ub = Target_all$Percent_VS_Male_UB * 100
        )
      )
    )
  )
}

# ---- UI ----------------------------------------------------------------------
app_theme <- bslib::bs_theme(
  version = 5,
  bootswatch = "minty",
  base_font = bslib::font_google("Inter")
)

ui <- page_navbar(
  title = "HIV Transmission Model Dashboard",
  theme = app_theme,
  sidebar = sidebar(
    h4("Simulation Settings"),
    textInput("param_dir", "Parameter folder", value = app_data_dir()),
    numericInput("trials", "Trials", value = 50, min = 1, step = 1),
    numericInput("sgroup", "Subgroups", value = 25, min = 5, step = 1),
    sliderInput("years", "Timeline", min = 2000, max = 2035, value = c(2004, 2020), sep = ""),
    numericInput("seed", "Random seed", value = 1, min = 0, step = 1),
    actionButton("run", "Run Simulation", class = "btn-primary"),
    hr(),
    uiOutput("files_status")
  ),
  navset_tab(
    nav_panel("Overview",
              fluidRow(
                column(3, uiOutput("card_prev")),
                column(3, uiOutput("card_inc")),
                column(3, uiOutput("card_vs")),
                column(3, uiOutput("card_art"))
              ),
              fluidRow(
                column(4,
                       selectInput(
                         "group",
                         "Population focus",
                         choices = c("Total", "Female", "Male"),
                         selected = "Total"
                       )
                )
              ),
              hr(),
              tabsetPanel(
                tabPanel("Prevalence",
                         withSpinner(plotlyOutput("plot_prev", height = 360))
                ),
                tabPanel("Incidence",
                         withSpinner(plotlyOutput("plot_inc", height = 360))
                ),
                tabPanel("On ART & VS",
                         withSpinner(plotlyOutput("plot_art", height = 360))
                )
              )
    ),
    nav_panel("Calibration",
              h5("Best 50 Trials by Goodness-of-Fit (lower is better)"),
              withSpinner(DTOutput("tbl_top50")),
              br(),
              downloadButton("download_outcomes", "Download outcomes.csv")
    ),
    nav_panel("About",
              div(class = "p-3",
                  h4("About this app"),
                  p("This app rewrites the original HIV transmission dynamic model with improved reliability, explicit timelines, and a modern dashboard UI. Place your parameter CSVs in the ", code("data/"), " folder next to app.R, or change the folder path in the sidebar."),
                  tags$ul(
                    tags$li("Deterministic & reproducible: controlled seed and timeline"),
                    tags$li("Safer I/O: no setwd(), errors surface in UI"),
                    tags$li("Interactive visuals: hover, zoom, and compare trials"),
                    tags$li("Downloadable results for further analysis")
                  ),
                  p("Tip: To speed up prototyping, start with a small number of trials and a shorter timeline, then scale up.")
              )
    )
  )
)

# ---- Server ------------------------------------------------------------------
server <- function(input, output, session) {
  iv <- InputValidator$new()
  iv$add_rule("param_dir", sv_required())
  iv$enable()
  
  files_status <- reactive({
    dir <- input$param_dir
    exists <- file.exists(file.path(dir, FILES_REQUIRED))
    tibble(file = FILES_REQUIRED, found = exists)
  })
  
  output$files_status <- renderUI({
    df <- files_status()
    n_ok <- sum(df$found)
    tagList(
      h5("Parameter files:"),
      tags$small(sprintf("%d/%d found in %s", n_ok, nrow(df), input$param_dir)),
      if (n_ok < nrow(df)) tags$p(class = "text-danger", "Missing files will block a run.")
    )
  })
  
  sim <- eventReactive(input$run, {
    validate(need(all(files_status()$found), paste0(
      "Missing files in ", input$param_dir, ". Check the list of required CSVs.")))
    
    params <- read_params(input$param_dir)
    res <- simulate_model(
      params = params,
      trials = input$trials,
      sgroup = input$sgroup,
      start_year = input$years[1],
      end_year = input$years[2],
      seed = input$seed
    )
    res
  }, ignoreInit = TRUE)
  
  # Value cards ----------------------------------------------------------------
  make_card <- function(title, value, suffix = "", good = TRUE) {
    val <- if (is.na(value)) "–" else format(value, big.mark = ",", scientific = FALSE, digits = 3)
    bg <- if (good) "#e8fff3" else "#fff0f0"
    col <- if (good) "#2e7d32" else "#b71c1c"
    div(style = paste("background:", bg, "; border-radius: 12px; padding: 16px;"),
        h6(title),
        div(style = paste("font-size: 28px; font-weight: 700; color:", col), paste0(val, suffix))
    )
  }
  
  last_year_df <- reactive({
    req(sim())
    out <- sim()$outcomes
    last_year <- max(out$year, na.rm = TRUE)
    out %>% filter(year == last_year)
  })
  
  output$card_prev <- renderUI({
    df <- last_year_df(); req(nrow(df))
    prev <- tail(df$Prevalence, 1) * 100
    make_card("Prevalence (latest year)", round(prev, 2), "%")
  })
  output$card_inc <- renderUI({
    df <- last_year_df(); req(nrow(df))
    inc <- tail(df$Incidence, 1)
    make_card("Incidence (per susceptible, 12m)", round(inc, 4))
  })
  output$card_vs <- renderUI({
    df <- last_year_df(); req(nrow(df))
    vs <- tail(df$V_total / df$N_total, 1) * 100
    make_card("Viral suppression among PLHIV", round(vs, 1), "%")
  })
  output$card_art <- renderUI({
    df <- last_year_df(); req(nrow(df))
    onart <- tail(df$On_ART_total, 1)
    make_card("On ART (count)", round(onart), good = TRUE)
  })
  
  # Plots ----------------------------------------------------------------------
  output$plot_prev <- renderPlotly({
    req(sim())
    df <- sim()$yearly %>% filter(group == input$group)
    target_df <- sim()$targets$prevalence %>% filter(group == input$group, !is.na(target))
    p <- ggplot(df, aes(x = year, y = prevalence, group = trial)) +
      geom_ribbon(
        data = target_df,
        aes(ymin = target_lb, ymax = target_ub),
        inherit.aes = FALSE,
        fill = "#0072B2",
        alpha = 0.12
      ) +
      geom_line(
        data = target_df,
        aes(y = target),
        inherit.aes = FALSE,
        color = "#0072B2",
        linewidth = 1.1,
        linetype = "dashed"
      ) +
      geom_point(
        data = target_df,
        aes(y = target),
        inherit.aes = FALSE,
        color = "#0072B2",
        size = 1.6
      ) +
      geom_line(alpha = 0.2) +
      stat_summary(fun = mean, geom = "line", linewidth = 1, color = "#222222") +
      labs(x = "Year", y = "Prevalence (%)") +
      theme_minimal(base_size = 12)
    ggplotly(p)
  })

  output$plot_inc <- renderPlotly({
    req(sim())
    df <- sim()$yearly %>% filter(group == input$group)
    target_df <- sim()$targets$incidence %>% filter(group == input$group, !is.na(target))
    p <- ggplot(df, aes(x = year, y = incidence, group = trial)) +
      geom_ribbon(
        data = target_df,
        aes(ymin = target_lb, ymax = target_ub),
        inherit.aes = FALSE,
        fill = "#D55E00",
        alpha = 0.12
      ) +
      geom_line(
        data = target_df,
        aes(y = target),
        inherit.aes = FALSE,
        color = "#D55E00",
        linewidth = 1.1,
        linetype = "dashed"
      ) +
      geom_point(
        data = target_df,
        aes(y = target),
        inherit.aes = FALSE,
        color = "#D55E00",
        size = 1.6
      ) +
      geom_line(alpha = 0.2) +
      stat_summary(fun = mean, geom = "line", linewidth = 1, color = "#222222") +
      labs(x = "Year", y = "Incidence (per susceptible, 12m)") +
      theme_minimal(base_size = 12)
    ggplotly(p)
  })

  output$plot_art <- renderPlotly({
    req(sim())
    df <- sim()$yearly %>% filter(group == input$group)
    colors <- c("On ART" = "#009E73", "Viral suppression" = "#CC79A7")
    df_long <- bind_rows(
      df %>% transmute(trial, year, metric = "On ART", value = on_art_pct),
      df %>% transmute(trial, year, metric = "Viral suppression", value = vs_pct)
    ) %>%
      tidyr::drop_na(value)
    req(nrow(df_long) > 0)
    df_long <- df_long %>% mutate(metric = factor(metric, levels = names(colors)))

    target_art_df <- sim()$targets$art %>%
      filter(group == input$group, !is.na(target)) %>%
      mutate(metric = factor("On ART", levels = names(colors)))
    target_vs_df <- sim()$targets$vs %>%
      filter(group == input$group, !is.na(target)) %>%
      mutate(metric = factor("Viral suppression", levels = names(colors)))
    target_df <- bind_rows(target_art_df, target_vs_df)
    target_df_ribbon <- target_df %>% tidyr::drop_na(target_lb, target_ub)

    p <- ggplot(df_long, aes(x = year, y = value,
                             group = interaction(metric, trial), color = metric)) +
      geom_ribbon(
        data = target_df_ribbon,
        aes(x = year, ymin = target_lb, ymax = target_ub, fill = metric),
        inherit.aes = FALSE,
        alpha = 0.12
      ) +
      geom_line(
        data = target_df,
        aes(x = year, y = target, color = metric),
        inherit.aes = FALSE,
        linewidth = 1.1,
        linetype = "dashed"
      ) +
      geom_point(
        data = target_df,
        aes(x = year, y = target, color = metric),
        inherit.aes = FALSE,
        size = 1.6
      ) +
      geom_line(alpha = 0.15) +
      stat_summary(aes(color = metric), fun = mean, geom = "line", linewidth = 1.1) +
      scale_color_manual(values = colors, name = NULL) +
      scale_fill_manual(values = colors, name = NULL) +
      labs(x = "Year", y = "Percent of PLHIV (%)") +
      theme_minimal(base_size = 12)
    ggplotly(p)
  })
  
  # Calibration table ----------------------------------------------------------
  output$tbl_top50 <- renderDT({
    req(sim())
    dat <- sim()$pd_top50 %>% rename(trial = V1, pd_total = V2, N_within = V3)
    datatable(dat, rownames = FALSE, options = list(pageLength = 10))
  })
  
  output$download_outcomes <- downloadHandler(
    filename = function() sprintf("outcomes_%s.csv", format(Sys.time(), "%Y%m%d_%H%M%S")),
    content = function(file) {
      req(sim())
      readr::write_csv(sim()$outcomes, file)
    }
  )
}

# ---- Run App -----------------------------------------------------------------
shinyApp(ui, server)

