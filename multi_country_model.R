#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(glue)
  library(rlang)
})

get_script_path <- function() {
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  }
  if (!is.null(attr(body(sys.function()), "srcref"))) {
    return(dirname(normalizePath(attr(body(sys.function()), "srcref")[[1]]@srcfile$filename)))
  }
  return(getwd())
}

default_data_root <- function() {
  script_dir <- get_script_path()
  if (basename(script_dir) == "") {
    return(file.path(getwd(), "data"))
  }
  candidate <- file.path(script_dir, "data")
  if (dir.exists(candidate)) {
    return(candidate)
  }
  file.path(getwd(), "data")
}

load_country_inputs <- function(country, data_root) {
  country_path <- file.path(data_root, country)
  if (!dir.exists(country_path)) {
    stop(glue("No data directory found for '{country}'. Checked {country_path}"))
  }
  disease_path <- file.path(country_path, "disease_burden.csv")
  prevention_path <- file.path(country_path, "prevention.csv")
  targets_path <- file.path(country_path, "targets.csv")

  required_files <- c(disease_path, prevention_path, targets_path)
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0) {
    stop(glue("Missing required files for {country}: {paste(basename(missing_files), collapse = ', ')}"))
  }

  disease <- read_csv(disease_path, show_col_types = FALSE)
  prevention <- read_csv(prevention_path, show_col_types = FALSE)
  targets <- read_csv(targets_path, show_col_types = FALSE)

  list(
    disease = disease,
    prevention = prevention,
    targets = targets,
    country_path = country_path
  )
}

prepare_country_data <- function(inputs) {
  merged <- inputs$disease %>%
    inner_join(inputs$prevention, by = "year") %>%
    arrange(year)
  list(
    merged = merged,
    targets = inputs$targets,
    country_path = inputs$country_path
  )
}

simulate_metrics <- function(data, calibration) {
  calibration <- as.list(calibration)
  effect_transmission <- exp(calibration$transmission_multiplier * (1 - data$condom_use))
  effect_art <- exp(-calibration$art_effect * data$art_coverage)
  effect_prep <- exp(-calibration$prep_effect * data$prep_coverage)

  prevalence <- pmin(pmax(data$baseline_prevalence * effect_transmission * effect_art, 1e-4), 0.99)
  incidence <- pmax(data$baseline_incidence * effect_transmission * effect_art * effect_prep, 1e-6)
  mortality <- pmax(data$baseline_mortality * exp(-0.5 * calibration$art_effect * data$art_coverage), 1e-6)
  art_coverage <- pmin(pmax(data$art_coverage * (1 + calibration$art_scale), 0), 0.99)

  tibble(
    year = data$year,
    prevalence = prevalence,
    incidence = incidence,
    mortality = mortality,
    art_coverage = art_coverage
  )
}

objective_function <- function(params, data, targets) {
  calibration <- list(
    transmission_multiplier = params[[1]],
    art_effect = params[[2]],
    prep_effect = params[[3]],
    art_scale = params[[4]]
  )
  predictions <- simulate_metrics(data, calibration)

  predictions_long <- predictions %>%
    select(year, prevalence, incidence, art_coverage) %>%
    pivot_longer(cols = -year, names_to = "indicator", values_to = "model_value")

  merged <- targets %>%
    semi_join(predictions_long, by = c("indicator", "year")) %>%
    inner_join(predictions_long, by = c("indicator", "year"))

  if (nrow(merged) == 0) {
    return(1e6)
  }

  scaling <- ifelse(!is.na(merged$upper) & !is.na(merged$lower) & merged$upper != merged$lower,
                    merged$upper - merged$lower,
                    pmax(abs(merged$value), 1e-6))

  mean(((merged$model_value - merged$value) / scaling) ^ 2)
}

calibrate_country <- function(data, targets) {
  start <- c(transmission_multiplier = 0.4, art_effect = 1.0, prep_effect = 0.5, art_scale = 0.1)
  lower <- c(-1.5, 0.01, 0.01, -0.5)
  upper <- c(2.5, 4.0, 2.0, 0.5)

  opt <- optim(
    par = start,
    fn = objective_function,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    data = data,
    targets = targets
  )

  list(
    calibration = list(
      transmission_multiplier = opt$par[[1]],
      art_effect = opt$par[[2]],
      prep_effect = opt$par[[3]],
      art_scale = opt$par[[4]]
    ),
    value = opt$value,
    convergence = opt$convergence,
    message = opt$message
  )
}

render_plots <- function(predictions, targets, country, output_dir) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  indicator_labels <- c(
    prevalence = "HIV prevalence",
    incidence = "HIV incidence",
    art_coverage = "ART coverage"
  )

  for (indicator in names(indicator_labels)) {
    indicator_targets <- targets %>% filter(indicator == .env$indicator)
    indicator_predictions <- predictions %>% select(year, !!sym(indicator))

    plot_data <- indicator_predictions %>%
      rename(model = !!sym(indicator))

    p <- ggplot(plot_data, aes(x = year)) +
      geom_line(aes(y = model, color = "Model"), linewidth = 1) +
      geom_point(data = indicator_targets, aes(y = value, color = "Target"), size = 2) +
      geom_errorbar(data = indicator_targets,
                    aes(ymin = lower, ymax = upper, color = "Target"), width = 0.2) +
      scale_color_manual(values = c(Model = "#1f78b4", Target = "#e31a1c"), name = NULL) +
      labs(title = glue("{country}: {indicator_labels[[indicator]]}"),
           x = "Year",
           y = "Rate") +
      theme_minimal() +
      theme(legend.position = "bottom")

    filename <- file.path(output_dir, glue("{country}_{indicator}.png"))
    ggsave(filename, plot = p, width = 7, height = 4, dpi = 300)
  }
}

produce_summary <- function(predictions, targets) {
  predictions_long <- predictions %>%
    pivot_longer(cols = -year, names_to = "indicator", values_to = "model_value")

  targets_summary <- targets %>%
    group_by(indicator) %>%
    summarise(
      observed_mean = mean(value, na.rm = TRUE),
      observed_latest = value[which.max(year)],
      .groups = "drop"
    )

  model_summary <- predictions_long %>%
    group_by(indicator) %>%
    summarise(
      model_mean = mean(model_value, na.rm = TRUE),
      model_latest = model_value[which.max(year)],
      .groups = "drop"
    )

  full_join(targets_summary, model_summary, by = "indicator")
}

main <- function() {
  option_list <- list(
    make_option(c("-c", "--country"), type = "character", default = "rwanda",
                help = "Country to calibrate. Must match a folder name inside the data directory."),
    make_option(c("-d", "--data-root"), type = "character", default = default_data_root(),
                help = "Root directory that contains per-country subfolders with input CSV files."),
    make_option(c("-o", "--output-dir"), type = "character", default = "outputs",
                help = "Directory where plots and summaries will be written."),
    make_option(c("--print-summary"), action = "store_true", default = FALSE,
                help = "Print the calibration summary table to the console.")
  )

  parser <- OptionParser(option_list = option_list)
  args <- parse_args(parser)

  inputs <- load_country_inputs(args$country, args$`data-root`)
  prepared <- prepare_country_data(inputs)

  calibration <- calibrate_country(prepared$merged, prepared$targets)
  if (calibration$convergence != 0) {
    warning(glue("Calibration did not converge cleanly (code {calibration$convergence}): {calibration$message}"))
  }

  predictions <- simulate_metrics(prepared$merged, calibration$calibration)

  output_dir <- file.path(args$`output-dir`, args$country)
  render_plots(predictions, prepared$targets, args$country, output_dir)

  summary_tbl <- produce_summary(predictions, prepared$targets)
  summary_path <- file.path(output_dir, glue("{args$country}_summary.csv"))
  write_csv(summary_tbl, summary_path)

  if (isTRUE(args$`print-summary`)) {
    print(summary_tbl)
  }

  message(glue("Calibration complete for {args$country}."))
  message(glue("Plots and summary saved to {normalizePath(output_dir)}"))
}

if (identical(environment(), globalenv())) {
  main()
}
