################################################################################
# Simplified Bayesian analysis for scratching and tick removal
#
# Purpose
# -------
# This script analyzes each tick species / worksheet separately.
# For each worksheet, it fits TWO straightforward Bayesian models:
#
#   1. Scratch_bouts
#      - Outcome type: count
#      - Model: negative-binomial regression
#      - Why: scratch bouts are counts and may be overdispersed
#        relative to a simple Poisson model.
#
#   2. Ticks_remaining
#      - Outcome type: count out of a known starting number of ticks
#      - Model: beta-binomial regression with 15 trials
#      - Why: ticks remaining is bounded from 0 to 15, not an unrestricted count.
#
# If a worksheet contains repeated measurements from the same Guinea_Pig, the model
# includes a Guinea_Pig random intercept:
#
#      outcome ~ Infestation + (1 | Guinea_Pig)
#
# If Guinea_Pig is not repeated in a worksheet, the model automatically drops the
# random intercept:
#
#      outcome ~ Infestation
#
# This is simpler and easier to explain than a fully multivariate Bayesian MANOVA.
# It still directly answers the biological questions:
#
#   - Does scratching differ by infestation?
#   - Do ticks remaining out of 15 differ by infestation?
#   - Which infestation groups differ from each other?
#
# Output
# ------
# The script creates an output folder with:
#
#   - predicted_means_all_datasets.csv
#       Expected mean scratch bouts and expected ticks remaining by Infestation.
#
#   - pairwise_differences_all_datasets.csv
#       Pairwise contrasts between infestation levels, on the original outcome scale.
#       For example: infestation B minus infestation A.
#
#   - model_diagnostics_all_datasets.csv
#       Basic MCMC diagnostics: Rhat, ESS, divergences, treedepth hits.
#
#   - posterior_predictive_check_*.png
#       Plots to check whether the fitted models reproduce the observed data.
#
#   - fitted_model_*.rds
#       Saved brms model objects, so results can be inspected later.
#
################################################################################


################################################################################
# 1. User settings
################################################################################

# Excel workbook to analyze.
# Put this R script in the same folder as the workbook, or use a full path.
input_file <- paste(here::here(), RawDat, RawDatSub, "Correlating Scratch and Removal.xlsx", sep = "/")

# Starting number of ticks per animal.
# You said this is usually 15. If a future experiment has row-specific starting
# tick counts, add a Ticks_initial column to the workbook and modify the code below.
tick_start_count <- 15L

# Folder where all results will be saved.
output_dir <- "bayesian_simple_output_ticks15"

# MCMC settings.
# For final reporting, keep 4 chains and at least 3000 to 4000 total iterations.
# For a quick test run, you can temporarily reduce iter to 1000.
chains <- 4
iter <- 4000
warmup <- 1000
seed <- 20260518

# Bayesian models sometimes need a stricter setting to avoid divergent transitions.
# If the diagnostics report divergences, increase adapt_delta to 0.99.
adapt_delta <- 0.95
max_treedepth <- 12

# Biomedical interpretation thresholds for pairwise differences.
# These are not statistical assumptions. They define what you consider a negligible
# effect size on the original scale.
#
# Example: if the scratch difference is between -1 and +1 bout, we call it
# practically small.
rope_scratch_bouts <- 1

# Example: if the tick difference is between -1 and +1 tick out of 15, we call it
# practically small.
rope_ticks_remaining <- 1

# If TRUE, the script will reuse saved fitted models if they already exist.
# Set this to FALSE if you change the model formula, priors, or data and want to refit.
reuse_saved_models <- TRUE


################################################################################
# 2. Packages
################################################################################

# Install packages once if needed:
# install.packages(c("readxl", "dplyr", "tidyr", "ggplot2", "brms",
#                    "posterior", "bayesplot"))
#
# Recommended Stan backend for brms:
# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/",
#                                        getOption("repos")))
# cmdstanr::install_cmdstan()

required_packages <- c(
  "readxl",
  "dplyr",
  "tidyr",
  "ggplot2",
  "brms",
  "posterior",
  "bayesplot",
  "tibble"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Package '", pkg, "' is required but not installed.\n",
      "Install it with install.packages('", pkg, "') and rerun the script."
    )
  }
}

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(brms)
library(posterior)
library(bayesplot)
library(tibble)

# Use cmdstanr if available; otherwise brms will use rstan if it is installed.
stan_backend <- if (requireNamespace("cmdstanr", quietly = TRUE)) "cmdstanr" else "rstan"

# Create output folders.
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "models"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "plots"), showWarnings = FALSE, recursive = TRUE)


################################################################################
# 3. Helper functions for data cleaning and validation
################################################################################

# This helper makes column matching robust to capitalization and punctuation.
# For example, "Guinea_Pig", "guinea pig", and "Guinea.Pig" are all treated as
# the same name.
clean_key <- function(x) {
  tolower(gsub("[^a-z0-9]", "", x))
}

# Find a required column in a data frame, even if capitalization/punctuation differs.
find_column <- function(dat, target_name) {
  target_key <- clean_key(target_name)
  name_keys <- clean_key(names(dat))
  hit <- which(name_keys == target_key)

  if (length(hit) == 0) {
    stop(
      "Required column '", target_name, "' was not found. Columns present are: ",
      paste(names(dat), collapse = ", ")
    )
  }

  if (length(hit) > 1) {
    stop(
      "Column '", target_name, "' matched more than one column: ",
      paste(names(dat)[hit], collapse = ", ")
    )
  }

  names(dat)[hit]
}

# Read one worksheet and standardize its columns.
read_one_sheet <- function(input_file, sheet_name, tick_start_count) {
  raw <- readxl::read_excel(input_file, sheet = sheet_name)

  # Locate columns. This allows minor spelling/capitalization differences.
  col_infestation <- find_column(raw, "Infestation")
  col_gp <- find_column(raw, "Guinea_Pig")
  col_scratch <- find_column(raw, "Scratch_bouts")
  col_ticks <- find_column(raw, "Ticks_remaining")

  dat <- raw %>%
    transmute(
      Dataset = sheet_name,
      Guinea_Pig = as.factor(.data[[col_gp]]),
      Infestation_raw = .data[[col_infestation]],
      Scratch_bouts = as.numeric(.data[[col_scratch]]),
      Ticks_remaining = as.numeric(.data[[col_ticks]])
    )

  # Preserve the natural infestation order.
  # If Infestation is numeric, use numeric order; otherwise use order of appearance.
  if (is.numeric(dat$Infestation_raw)) {
    infestation_levels <- sort(unique(dat$Infestation_raw))
  } else {
    infestation_levels <- unique(dat$Infestation_raw)
  }

  dat <- dat %>%
    mutate(
      Infestation = factor(Infestation_raw, levels = infestation_levels),
      Ticks_initial = tick_start_count
    ) %>%
    select(Dataset, Guinea_Pig, Infestation,
           Scratch_bouts, Ticks_remaining, Ticks_initial)

  # Remove rows that are incomplete for the variables used here.
  # It is better to be explicit than to let the model silently drop rows.
  before_n <- nrow(dat)
  dat <- dat %>%
    filter(
      !is.na(Guinea_Pig),
      !is.na(Infestation),
      !is.na(Scratch_bouts),
      !is.na(Ticks_remaining)
    )
  after_n <- nrow(dat)

  if (after_n < before_n) {
    message(sheet_name, ": removed ", before_n - after_n, " incomplete rows.")
  }

  validate_one_sheet(dat, sheet_name, tick_start_count)

  dat
}

# Check the most important data requirements before fitting models.
validate_one_sheet <- function(dat, sheet_name, tick_start_count) {
  if (nrow(dat) == 0) {
    stop(sheet_name, ": no usable rows after removing missing values.")
  }

  if (nlevels(dat$Infestation) < 2) {
    stop(sheet_name, ": Infestation must have at least two levels.")
  }

  # Scratch_bouts should be a non-negative count.
  if (any(dat$Scratch_bouts < 0)) {
    stop(sheet_name, ": Scratch_bouts contains negative values.")
  }

  if (any(abs(dat$Scratch_bouts - round(dat$Scratch_bouts)) > 1e-8)) {
    warning(sheet_name, ": Scratch_bouts contains non-integer values. A count model assumes counts.")
  }

  # Ticks_remaining should be an integer count between 0 and tick_start_count.
  if (any(dat$Ticks_remaining < 0 | dat$Ticks_remaining > tick_start_count)) {
    stop(
      sheet_name,
      ": Ticks_remaining must be between 0 and ", tick_start_count,
      " for the beta-binomial model."
    )
  }

  if (any(abs(dat$Ticks_remaining - round(dat$Ticks_remaining)) > 1e-8)) {
    warning(sheet_name, ": Ticks_remaining contains non-integer values. A beta-binomial model assumes counts.")
  }

  # Check for duplicate rows within an animal and infestation level.
  # Duplicates are not automatically wrong, but they can indicate accidental data entry.
  duplicates <- dat %>%
    count(Guinea_Pig, Infestation) %>%
    filter(n > 1)

  if (nrow(duplicates) > 0) {
    warning(
      sheet_name,
      ": found duplicate Guinea_Pig x Infestation rows. ",
      "If these are technical replicates, consider averaging or modeling replicate structure."
    )
  }
}

# Determine whether Guinea_Pig is repeated in this worksheet.
# If at least one Guinea_Pig occurs in more than one row, we use a random intercept.
has_repeated_guinea_pigs <- function(dat) {
  max(table(dat$Guinea_Pig)) > 1
}


################################################################################
# 4. Model definitions
################################################################################

# Build the model formula for Scratch_bouts.
# The negative-binomial model is used for counts that may be overdispersed.
make_scratch_formula <- function(repeated) {
  if (repeated) {
    brms::bf(Scratch_bouts ~ Infestation + (1 | Guinea_Pig))
  } else {
    brms::bf(Scratch_bouts ~ Infestation)
  }
}

# Build the model formula for Ticks_remaining.
# The beta-binomial model treats Ticks_remaining as successes out of 15 trials.
# Here "success" means "tick remained".
make_ticks_formula <- function(repeated) {
  if (repeated) {
    brms::bf(Ticks_remaining | trials(Ticks_initial) ~ Infestation + (1 | Guinea_Pig))
  } else {
    brms::bf(Ticks_remaining | trials(Ticks_initial) ~ Infestation)
  }
}

# Priors are intentionally simple and weakly informative.
# They stabilize small-sample models without forcing the result.
#
# For Scratch_bouts:
#   - Intercept is on the log-count scale.
#   - normal(log(5), 2) allows a very wide baseline mean range.
#   - Infestation coefficients are log ratios.
#   - normal(0, 1.5) allows large fold-changes but discourages extreme estimates.
#
# For Ticks_remaining:
#   - Intercept is on the logit probability scale.
#   - normal(0, 2) is centered around 50% remaining but very wide.
#   - Infestation coefficients are log odds ratios.
#
# For repeated datasets:
#   - Random-effect SD prior regularizes animal-to-animal variability.
make_priors <- function(outcome, repeated) {
  if (outcome == "Scratch_bouts") {
    priors <- c(
      brms::prior(normal(log(5), 2), class = "Intercept"),
      brms::prior(normal(0, 1.5), class = "b")
    )
  } else if (outcome == "Ticks_remaining") {
    priors <- c(
      brms::prior(normal(0, 2), class = "Intercept"),
      brms::prior(normal(0, 1.5), class = "b")
    )
  } else {
    stop("Unknown outcome: ", outcome)
  }

  if (repeated) {
    priors <- c(priors, brms::prior(exponential(1), class = "sd"))
  }

  priors
}


################################################################################
# 5. Fit one Bayesian model
################################################################################

fit_one_model <- function(dat, dataset_name, outcome, formula, family, priors) {
  model_file <- file.path(
    output_dir,
    "models",
    paste0("fitted_model_", dataset_name, "_", outcome, ".rds")
  )

  if (reuse_saved_models && file.exists(model_file)) {
    message("Reusing saved model: ", model_file)
    return(readRDS(model_file))
  }

  message("Fitting ", dataset_name, " - ", outcome, " ...")

  fit <- brms::brm(
    formula = formula,
    data = dat,
    family = family,
    prior = priors,
    chains = chains,
    iter = iter,
    warmup = warmup,
    seed = seed,
    backend = stan_backend,
    control = list(
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth
    ),
    refresh = 500
  )

  saveRDS(fit, model_file)
  fit
}


################################################################################
# 6. Summaries that biomedical readers usually want
################################################################################

# Summarize posterior draws using median and 95% credible interval.
summarize_draws_simple <- function(x) {
  tibble::tibble(
    mean = mean(x),
    median = stats::median(x),
    lower_95 = stats::quantile(x, 0.025),
    upper_95 = stats::quantile(x, 0.975),
    prob_greater_0 = mean(x > 0),
    prob_less_0 = mean(x < 0)
  )
}

# Create expected mean predictions for each Infestation level.
#
# Important:
#   re_formula = NA means predictions are population-average estimates.
#   In plain language: the estimates are not for one particular guinea pig;
#   they are for a typical animal from the studied population.
predicted_means_by_infestation <- function(fit, dat, outcome) {
  infestation_levels <- levels(dat$Infestation)

  newdata <- data.frame(
    Infestation = factor(infestation_levels, levels = infestation_levels)
  )

  # Include a placeholder Guinea_Pig level. If the model does not use it,
  # brms ignores this extra column. If the model does use it, it is ignored
  # because re_formula = NA removes animal-specific random effects.
  newdata$Guinea_Pig <- levels(dat$Guinea_Pig)[1]

  # The tick model needs the number of trials.
  if (outcome == "Ticks_remaining") {
    newdata$Ticks_initial <- tick_start_count
  }

  epred <- brms::posterior_epred(
    fit,
    newdata = newdata,
    re_formula = NA
  )

  # epred is a matrix: posterior draws x infestation levels.
  result <- lapply(seq_along(infestation_levels), function(j) {
    draws_j <- epred[, j]

    tibble::tibble(
      Infestation = infestation_levels[j],
      outcome = outcome,
      expected_mean = mean(draws_j),
      median = stats::median(draws_j),
      lower_95 = stats::quantile(draws_j, 0.025),
      upper_95 = stats::quantile(draws_j, 0.975)
    )
  }) %>%
    bind_rows()

  # For ticks, also provide the proportion remaining out of 15.
  if (outcome == "Ticks_remaining") {
    result <- result %>%
      mutate(
        expected_proportion_remaining = expected_mean / tick_start_count,
        lower_95_proportion = lower_95 / tick_start_count,
        upper_95_proportion = upper_95 / tick_start_count
      )
  }

  result
}

# Pairwise comparisons between Infestation levels.
#
# The comparison is always:
#
#   second level - first level
#
# Example:
#   If the contrast is "Infestation_2 - Infestation_1" and the median difference
#   is -3 ticks, then Infestation_2 has about 3 fewer ticks remaining than
#   Infestation_1.
pairwise_differences <- function(fit, dat, outcome, rope) {
  infestation_levels <- levels(dat$Infestation)

  newdata <- data.frame(
    Infestation = factor(infestation_levels, levels = infestation_levels)
  )

  # Include a placeholder Guinea_Pig level. See the comment in
  # predicted_means_by_infestation() for why this is safe.
  newdata$Guinea_Pig <- levels(dat$Guinea_Pig)[1]

  if (outcome == "Ticks_remaining") {
    newdata$Ticks_initial <- tick_start_count
  }

  epred <- brms::posterior_epred(
    fit,
    newdata = newdata,
    re_formula = NA
  )

  comparisons <- utils::combn(infestation_levels, 2, simplify = FALSE)

  result <- lapply(comparisons, function(pair) {
    level_1 <- pair[1]
    level_2 <- pair[2]

    i <- which(infestation_levels == level_1)
    j <- which(infestation_levels == level_2)

    diff_draws <- epred[, j] - epred[, i]

    out <- summarize_draws_simple(diff_draws) %>%
      mutate(
        outcome = outcome,
        level_1 = level_1,
        level_2 = level_2,
        contrast = paste0(level_2, " - ", level_1),
        scale = "original outcome scale",
        rope = rope,
        prob_practically_small = mean(abs(diff_draws) <= rope),
        .before = 1
      )

    # For ticks, add the difference as a proportion of the 15 starting ticks.
    if (outcome == "Ticks_remaining") {
      prop_diff_draws <- diff_draws / tick_start_count
      prop_out <- summarize_draws_simple(prop_diff_draws) %>%
        mutate(
          outcome = outcome,
          level_1 = level_1,
          level_2 = level_2,
          contrast = paste0(level_2, " - ", level_1),
          scale = "proportion of 15 ticks",
          rope = rope / tick_start_count,
          prob_practically_small = mean(abs(prop_diff_draws) <= rope / tick_start_count),
          .before = 1
        )

      out <- bind_rows(out, prop_out)
    }

    out
  }) %>%
    bind_rows()

  result
}


################################################################################
# 7. Diagnostics
################################################################################

# Good Bayesian practice requires checking computation and model fit.
#
# Key diagnostic targets:
#   - Rhat should usually be very close to 1.00.
#   - Effective sample sizes should not be tiny.
#   - Divergent transitions should be 0.
#   - Posterior predictive checks should show that simulated data look like the
#     observed data.
extract_diagnostics <- function(fit, dataset_name, outcome) {
  draws_summary <- posterior::summarise_draws(posterior::as_draws_df(fit))

  max_rhat <- suppressWarnings(max(draws_summary$rhat, na.rm = TRUE))
  min_ess_bulk <- suppressWarnings(min(draws_summary$ess_bulk, na.rm = TRUE))
  min_ess_tail <- suppressWarnings(min(draws_summary$ess_tail, na.rm = TRUE))

  nuts <- brms::nuts_params(fit)

  n_divergent <- nuts %>%
    filter(Parameter == "divergent__") %>%
    dplyr::summarize(n = sum(Value == 1)) %>%
    pull(n)

  n_max_treedepth <- nuts %>%
    filter(Parameter == "treedepth__") %>%
    dplyr::summarize(n = sum(Value >= max_treedepth)) %>%
    pull(n)

  tibble::tibble(
    Dataset = dataset_name,
    outcome = outcome,
    max_Rhat = max_rhat,
    min_ESS_bulk = min_ess_bulk,
    min_ESS_tail = min_ess_tail,
    divergent_transitions = n_divergent,
    max_treedepth_hits = n_max_treedepth,
    suggested_action = dplyr::case_when(
      n_divergent > 0 ~ "Divergences detected: refit with adapt_delta = 0.99 and inspect priors/model.",
      max_rhat > 1.01 ~ "Rhat > 1.01: increase iterations and inspect trace plots.",
      min_ess_bulk < 400 ~ "Low ESS: increase iterations.",
      TRUE ~ "No major MCMC warning from these checks."
    )
  )
}

# Save a posterior predictive check plot.
# The plot compares your observed data with fake datasets simulated from the model.
# If the model is reasonable, the observed pattern should look compatible with the
# simulated patterns.
save_pp_check_plot <- function(fit, dataset_name, outcome) {
  plot_file <- file.path(
    output_dir,
    "plots",
    paste0("posterior_predictive_check_", dataset_name, "_", outcome, ".png")
  )

  p <- brms::pp_check(fit, type = "bars", ndraws = 100) +
    ggplot2::ggtitle(paste(dataset_name, outcome, "posterior predictive check"))

  ggplot2::ggsave(
    filename = plot_file,
    plot = p,
    width = 8,
    height = 5,
    dpi = 300
  )

  plot_file
}


################################################################################
# 8. Run analysis for one worksheet
################################################################################

analyze_one_dataset <- function(dat) {
  dataset_name <- unique(dat$Dataset)
  repeated <- has_repeated_guinea_pigs(dat)

  message("\n==============================")
  message("Dataset: ", dataset_name)
  message("Rows: ", nrow(dat))
  message("Infestation levels: ", paste(levels(dat$Infestation), collapse = ", "))
  message("Repeated Guinea_Pig measurements: ", repeated)
  message("==============================\n")

  # ------------------------------
  # Scratch_bouts model
  # ------------------------------
  scratch_formula <- make_scratch_formula(repeated)
  scratch_priors <- make_priors("Scratch_bouts", repeated)

  fit_scratch <- fit_one_model(
    dat = dat,
    dataset_name = dataset_name,
    outcome = "Scratch_bouts",
    formula = scratch_formula,
    family = brms::negbinomial(),
    priors = scratch_priors
  )

  scratch_means <- predicted_means_by_infestation(
    fit = fit_scratch,
    dat = dat,
    outcome = "Scratch_bouts"
  )

  scratch_pairwise <- pairwise_differences(
    fit = fit_scratch,
    dat = dat,
    outcome = "Scratch_bouts",
    rope = rope_scratch_bouts
  )

  scratch_diag <- extract_diagnostics(
    fit = fit_scratch,
    dataset_name = dataset_name,
    outcome = "Scratch_bouts"
  )

  scratch_pp_file <- save_pp_check_plot(
    fit = fit_scratch,
    dataset_name = dataset_name,
    outcome = "Scratch_bouts"
  )

  # ------------------------------
  # Ticks_remaining model
  # ------------------------------
  ticks_formula <- make_ticks_formula(repeated)
  ticks_priors <- make_priors("Ticks_remaining", repeated)

  fit_ticks <- fit_one_model(
    dat = dat,
    dataset_name = dataset_name,
    outcome = "Ticks_remaining",
    formula = ticks_formula,
    family = brms::beta_binomial(),
    priors = ticks_priors
  )

  ticks_means <- predicted_means_by_infestation(
    fit = fit_ticks,
    dat = dat,
    outcome = "Ticks_remaining"
  )

  ticks_pairwise <- pairwise_differences(
    fit = fit_ticks,
    dat = dat,
    outcome = "Ticks_remaining",
    rope = rope_ticks_remaining
  )

  ticks_diag <- extract_diagnostics(
    fit = fit_ticks,
    dataset_name = dataset_name,
    outcome = "Ticks_remaining"
  )

  ticks_pp_file <- save_pp_check_plot(
    fit = fit_ticks,
    dataset_name = dataset_name,
    outcome = "Ticks_remaining"
  )

  list(
    predicted_means = bind_rows(scratch_means, ticks_means) %>%
      mutate(Dataset = dataset_name, .before = 1),
    pairwise = bind_rows(scratch_pairwise, ticks_pairwise) %>%
      mutate(Dataset = dataset_name, .before = 1),
    diagnostics = bind_rows(scratch_diag, ticks_diag),
    pp_check_files = c(scratch_pp_file, ticks_pp_file)
  )
}


################################################################################
# 9. Main analysis
################################################################################

# Read every worksheet in the Excel file.
sheets <- readxl::excel_sheets(input_file)
message("Worksheets found: ", paste(sheets, collapse = ", "))

all_data <- lapply(sheets, function(s) {
  read_one_sheet(input_file, sheet_name = s, tick_start_count = tick_start_count)
})

# Analyze each worksheet separately.
results <- lapply(all_data, analyze_one_dataset)

# Combine outputs across all worksheets.
predicted_means_all <- bind_rows(lapply(results, `[[`, "predicted_means"))
pairwise_all <- bind_rows(lapply(results, `[[`, "pairwise"))
diagnostics_all <- bind_rows(lapply(results, `[[`, "diagnostics"))

# Save CSV files.
readr_needed <- requireNamespace("readr", quietly = TRUE)
if (readr_needed) {
  readr::write_csv(
    predicted_means_all,
    file.path(output_dir, "predicted_means_all_datasets.csv")
  )
  readr::write_csv(
    pairwise_all,
    file.path(output_dir, "pairwise_differences_all_datasets.csv")
  )
  readr::write_csv(
    diagnostics_all,
    file.path(output_dir, "model_diagnostics_all_datasets.csv")
  )
} else {
  utils::write.csv(
    predicted_means_all,
    file.path(output_dir, "predicted_means_all_datasets.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    pairwise_all,
    file.path(output_dir, "pairwise_differences_all_datasets.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    diagnostics_all,
    file.path(output_dir, "model_diagnostics_all_datasets.csv"),
    row.names = FALSE
  )
}

# Write a short interpretation guide for collaborators.
interpretation_text <- c(
  "Simplified Bayesian analysis guide",
  "==================================",
  "",
  "Models fitted:",
  "1. Scratch_bouts: negative-binomial model.",
  "2. Ticks_remaining: beta-binomial model with 15 starting ticks.",
  "3. Repeated Guinea_Pig measurements receive a random intercept.",
  "",
  "How to read predicted_means_all_datasets.csv:",
  "- expected_mean is the estimated mean outcome for each Infestation level.",
  "- lower_95 and upper_95 are the 95% credible interval.",
  "- For Ticks_remaining, values are ticks remaining out of 15.",
  "",
  "How to read pairwise_differences_all_datasets.csv:",
  "- contrast is level_2 minus level_1.",
  "- A positive difference means level_2 has a higher expected outcome than level_1.",
  "- A negative difference means level_2 has a lower expected outcome than level_1.",
  "- prob_greater_0 is the posterior probability that the difference is positive.",
  "- prob_less_0 is the posterior probability that the difference is negative.",
  "- prob_practically_small is the probability that the difference is within the ROPE.",
  "",
  "Main assumptions:",
  "- Animals are independent experimental units.",
  "- Infestation is coded correctly and ordered as intended.",
  "- Scratch_bouts is a non-negative count.",
  "- Ticks_remaining is a count between 0 and 15.",
  "- Guinea_Pig identifies repeated measurements correctly.",
  "- The posterior predictive check plots should look compatible with the observed data.",
  "- MCMC diagnostics should show no divergences and Rhat near 1.00.",
  "",
  "If diagnostics are poor:",
  "- Increase adapt_delta to 0.99 and rerun.",
  "- Increase iter to 6000 or 8000 if ESS is low.",
  "- Check for data entry issues, outliers, or accidental duplicate rows.",
  "- Consider whether a zero-inflated model is needed if there are many zeros."
)

writeLines(
  interpretation_text,
  con = file.path(output_dir, "README_interpretation_guide.txt")
)

message("\nAnalysis complete.")
message("Results saved to: ", normalizePath(output_dir))
message("Key files:")
message("  - predicted_means_all_datasets.csv")
message("  - pairwise_differences_all_datasets.csv")
message("  - model_diagnostics_all_datasets.csv")
message("  - README_interpretation_guide.txt")
message("  - posterior predictive check plots in the plots/ folder")

################################################################################
# End of script
################################################################################
