# Project title: Forest growth strengthens diversity effects on
# host-parasitoid interactions
# Script title: Analysis helper functions
# Author: Massimo Martini
#
# Script description: This file is sourced by Design_models.R,
# Maximal_models.R, Maximal_SEM.R, Parsimonious_SEM.R, and Network_SEM.R. It
# contains shared functions for model checks and comparisons, residual
# diagnostics, planned
# inference and sensitivity summaries, and SEM standardization and R2
# calculations. Keeping these operations here applies the same checks and
# reporting rules throughout the analysis workflow. Sourcing this file does not
# itself read data or write analysis outputs.
#
# Function index
#
# Candidate-model comparison and validation
# - validate_candidate_models: check a candidate-model registry.
# - model_is_valid: apply the shared numerical fit rule.
# - compare_candidate_models: tabulate fit and information criteria.
# - summarize_fixed_terms: collect selected fixed-effect estimates.
# - run_candidate_diagnostics: run the shared DHARMa residual checks.
# - test_interaction_blocks: test hierarchy-preserving predictor blocks.
# - test_within_plot_temporal_autocorrelation: test residual time structure.
# - evaluate_maximal_structure_candidates: assemble a structure-selection audit.
# - validate_maximal_model: audit one selected maximal model.
#
# SEM inference and effect standardization
# - .sem_assert_scalar_character: validate one text identifier.
# - sem_term_wald_test: test a term at declared reference conditions.
# - .sem_apply_at: apply conditioning values to evaluation data.
# - .sem_weighted_mean: calculate a validated weighted mean.
# - .sem_predict: make model-class-aware predictions.
# - .sem_average_slope: calculate an average numerical derivative.
# - standardize_response_path: response-scale path standardization.
# - standardize_oe_path: Menard observed-empirical standardization.
# - standardize_relevant_range_path: relevant-range standardization.
# - .sem_node_row: retrieve one node-scale definition.
# - .sem_range_to_model_scale: transform biological range endpoints.
# - prepare_sem_relevant_ranges: validate and prepare the range registry.
# - standardize_sem_path: calculate all effect summaries for one path.
# - standardize_sem_edges: apply path standardization to an edge registry.
# - build_sem_component_results: assemble one complete SEM component result.
#
# Network-SEM focal inference and export
# - .sem_required_successes: convert a success threshold to a replicate count.
# - bootstrap_focal_dbi: bootstrap the bounded focal-richness component.
# - export_focal_dbi_bootstrap: save focal-bootstrap results and fit audits.
# - sem_focal_dbi_moments: calculate exact DBI count means and variances.
# - sem_predict_focal_dbi: supply DBI predictions to the SEM workflow.
# - sem_focal_dbi_inference_adapter: turn bootstrap draws into term inference.
# - sem_network_evaluation: enforce feasible network-richness contrasts.
# - sem_network_component_r2: calculate R2 for network component families.
# - network_sem_path_table: summarize the network SEM's directed paths.
# - export_network_sem: save a network-SEM handoff and readable result tables.
#
# Component-model R2 and cross-validation
# - .sem_weighted_variance: calculate unbiased weighted variance.
# - .sem_random_design_column: reconstruct a random-effect design column.
# - .sem_random_effect_variance: average random-effect variance by observation.
# - .sem_dispersion_prediction: obtain response-scale dispersion values.
# - .sem_binomial_trials: recover binomial trial totals.
# - .sem_conditional_response_variance: calculate family-specific variance.
# - sem_component_nakagawa_r2: calculate marginal and conditional R2.
# - make_grouped_sem_folds: assign intact plots to cross-validation folds.
# - sem_component_predictive_r2: calculate grouped predictive R2.
#
# Parsimonious-SEM reduction, testing, and export
# - compare_component_reduction: compare a reduced component with its parent.
# - run_component_reduction_plan: replay the recorded reduction decisions.
# - audit_component_missing_paths: test provisional removed paths.
# - assemble_reduced_sem_graph: reconstruct the finalized pathway graph.
# - derive_sem_basis_set: declare the graph's independence claims.
# - audit_sem_basis_set: test the final independence claims.
# - sem_fisher_c: combine missing-path tests into Fisher's C.
# - build_sem_path_table: summarize retained pathways.
# - export_parsimonious_sem_handoff: save RDS and CSV results.
#
# Shared refit validation
# - validate_sem_refit: check convergence and coefficient validity after refit.
#
# ------------------------------------------------------------------------------
# Shared model-comparison and diagnostic helpers
# ------------------------------------------------------------------------------
# These functions are used by Design_models.R and Maximal_models.R to compare
# fitted candidates and summarize model diagnostics.

# Validate a candidate-model list before it enters the shared workflow.
# `models` must be a non-empty, uniquely named list of glmmTMB fits.
# The function stops on invalid input and otherwise returns TRUE invisibly.
# Details: Confirms that the input is a non-empty, uniquely named list and
# that every entry is a fitted glmmTMB object before later comparisons run.
# Returns: TRUE invisibly; invalid registries stop with an informative error.
validate_candidate_models <- function(models) {
  if (!is.list(models) || length(models) == 0L) {
    stop("`models` must be a non-empty list of fitted glmmTMB models.")
  }

  model_names <- names(models)
  if (is.null(model_names) ||
      any(!nzchar(model_names)) ||
      anyDuplicated(model_names)) {
    stop("`models` must have non-empty, unique names.")
  }

  is_glmmTMB <- vapply(models, inherits, logical(1), what = "glmmTMB")
  if (any(!is_glmmTMB)) {
    invalid <- paste(
      sprintf(
        "%s <%s>",
        model_names[!is_glmmTMB],
        vapply(
          models[!is_glmmTMB],
          function(x) paste(class(x), collapse = "/"),
          character(1)
        )
      ),
      collapse = ", "
    )
    stop("Every candidate must be a fitted glmmTMB model. Invalid: ", invalid)
  }

  invisible(TRUE)
}


# Apply the numerical-validity rule used throughout the design analysis.
# A fit is valid only if it is a glmmTMB model, the optimizer reports
# convergence, and the Hessian is positive definite. Returns one logical value.
# Details: Centralizes the minimum numerical acceptance rule used throughout
# the workflow; it does not assess residual behavior or biological adequacy.
# Returns: One TRUE or FALSE value.
model_is_valid <- function(fit_object) {
  inherits(fit_object, "glmmTMB") &&
    isTRUE(fit_object$fit$convergence == 0L) &&
    isTRUE(fit_object$sdr$pdHess)
}


# Compare a named set of candidate models fitted to the same observations.
# The returned tibble records convergence, Hessian status, sample size,
# parameter count, log likelihood, AIC, AICc, BIC, and delta information
# criteria. Invalid fits retain their status but receive missing fit measures.
# Details: Computes likelihood and information criteria only for numerically
# valid fits, includes the small-sample AIC correction, and ranks valid models.
# Returns: One tibble row per candidate, ordered by validity and delta AIC.
compare_candidate_models <- function(models) {
  validate_candidate_models(models)

  comparison <- dplyr::bind_rows(
    Map(function(fit_object, model_name) {
      converged <- isTRUE(fit_object$fit$convergence == 0L)
      pd_hessian <- isTRUE(fit_object$sdr$pdHess)
      valid <- converged && pd_hessian
      model_loglik <- suppressWarnings(stats::logLik(fit_object))
      n <- stats::nobs(fit_object)
      k <- attr(model_loglik, "df")
      model_aic <- if (valid) stats::AIC(fit_object) else NA_real_
      model_bic <- if (valid) stats::BIC(fit_object) else NA_real_
      model_aicc <- if (valid && n > k + 1) {
        model_aic + (2 * k * (k + 1)) / (n - k - 1)
      } else {
        NA_real_
      }
      optimizer_message <- if (length(fit_object$fit$message) > 0L) {
        paste(fit_object$fit$message, collapse = "; ")
      } else {
        NA_character_
      }

      tibble::tibble(
        model = model_name,
        converged = converged,
        pd_hessian = pd_hessian,
        valid = valid,
        n = n,
        parameters = k,
        logLik = if (valid) as.numeric(model_loglik) else NA_real_,
        AIC = model_aic,
        AICc = model_aicc,
        BIC = model_bic,
        optimizer_message = optimizer_message
      )
    }, models, names(models))
  )

  valid_aic <- comparison$AIC[comparison$valid & is.finite(comparison$AIC)]
  best_aic <- if (length(valid_aic) > 0L) min(valid_aic) else NA_real_
  valid_aicc <- comparison$AICc[
    comparison$valid & is.finite(comparison$AICc)
  ]
  best_aicc <- if (length(valid_aicc) > 0L) min(valid_aicc) else NA_real_

  comparison |>
    dplyr::mutate(
      delta_AIC = ifelse(valid, AIC - best_aic, NA_real_),
      delta_AICc = ifelse(valid, AICc - best_aicc, NA_real_)
    ) |>
    dplyr::arrange(dplyr::desc(valid), delta_AIC)
}


# Extract selected conditional fixed effects from each candidate model.
# `terms` is a character vector of fitted coefficient names. The returned
# tibble contains estimates, standard errors, tests, 95% Wald intervals, and
# the candidate-model identifier.
# Details: Uses conditional-model coefficients and keeps only the requested
# fitted term names, allowing the same biological effects to be compared.
# Returns: A stacked tibble of estimates, tests, intervals, and model names.
summarize_fixed_terms <- function(models, terms) {
  validate_candidate_models(models)

  dplyr::bind_rows(
    lapply(names(models), function(nm) {
      broom.mixed::tidy(
        models[[nm]],
        effects = "fixed",
        component = "cond",
        conf.int = TRUE
      ) |>
        dplyr::filter(term %in% terms) |>
        dplyr::mutate(model = nm, .before = 1)
    })
  )
}


# Run a common set of DHARMa residual checks for each valid candidate model.
# `n_sim` controls the number of simulated residual datasets and `seed` fixes
# those simulations. Models with AR(1) or OU covariance use estimated rotation.
# The result contains a compact test-summary tibble and the simulated-residual
# objects. The bootstrap outlier test retains minor Monte Carlo variation.
# Details: Skips invalid fits, uses a reproducible model-specific seed, and
# tests uniformity, dispersion, outliers, zero inflation, and quantiles.
# Returns: A list containing a compact summary and each DHARMa residual object.
run_candidate_diagnostics <- function(models, n_sim = 1000L, seed = 123L) {
  validate_candidate_models(models)

  residuals <- setNames(
    lapply(seq_along(models), function(i) {
      model <- models[[i]]
      if (!model_is_valid(model)) {
        return(NULL)
      }

      has_temporal_covariance <- grepl(
        "\\b(ar1|ou)\\s*\\(",
        paste(deparse(stats::formula(model)), collapse = " "),
        perl = TRUE
      )

      DHARMa::simulateResiduals(
        fittedModel = model,
        n = n_sim,
        plot = FALSE,
        seed = seed + i - 1L,
        rotation = if (has_temporal_covariance) "estimated" else NULL
      )
    }),
    names(models)
  )

  test_summary <- dplyr::bind_rows(
    lapply(names(residuals), function(nm) {
      simulated <- residuals[[nm]]
      if (is.null(simulated)) {
        return(tibble::tibble(
          model = nm,
          uniformity_p = NA_real_,
          dispersion = NA_real_,
          dispersion_p = NA_real_,
          outlier_p = NA_real_,
          zero_ratio = NA_real_,
          zero_p = NA_real_,
          quantile_p = NA_real_
        ))
      }

      uniformity_test <- DHARMa::testUniformity(simulated, plot = FALSE)
      dispersion_test <- DHARMa::testDispersion(simulated, plot = FALSE)
      outlier_test <- DHARMa::testOutliers(
        simulated,
        type = if (isTRUE(simulated$integerResponse)) "bootstrap" else "binomial",
        nBoot = min(n_sim, 1000L),
        plot = FALSE
      )
      zero_test <- DHARMa::testZeroInflation(simulated, plot = FALSE)
      quantile_test <- suppressWarnings(
        DHARMa::testQuantiles(simulated, plot = FALSE)
      )
      scalar_or_na <- function(x) {
        if (length(x) == 0L) NA_real_ else as.numeric(x[[1]])
      }

      uniformity_p <- scalar_or_na(uniformity_test$p.value)
      dispersion_estimate <- scalar_or_na(dispersion_test$statistic)
      dispersion_p <- scalar_or_na(dispersion_test$p.value)
      outlier_p <- scalar_or_na(outlier_test$p.value)
      zero_ratio <- scalar_or_na(zero_test$statistic)
      zero_p <- scalar_or_na(zero_test$p.value)
      quantile_p <- scalar_or_na(quantile_test$p.value)

      tibble::tibble(
        model = nm,
        uniformity_p = uniformity_p,
        dispersion = dispersion_estimate,
        dispersion_p = dispersion_p,
        outlier_p = outlier_p,
        zero_ratio = zero_ratio,
        zero_p = zero_p,
        quantile_p = quantile_p
      )
    })
  )

  list(summary = test_summary, residuals = residuals)
}

# ------------------------------------------------------------------------------
# Maximal-model inference and validation
# ------------------------------------------------------------------------------

# Test interacting predictors as hierarchy-preserving likelihood-ratio blocks.
# For each eligible main effect, the reduced model removes that effect and all
# interactions containing it. `exclude` protects required design terms.
# Returns dropped terms, reduced fits, comparisons, and a P-value summary.
# Details: For each eligible main effect, refits a model without that effect
# and all interactions containing it, then applies an optional P adjustment.
# Returns: An invisible list of dropped terms, refits, LRTs, and a summary.
test_interaction_blocks <- function(model,
                                    exclude = character(0),
                                    alpha = 0.05,
                                    p_adjust_method = "none") {
  validate_candidate_models(list(maximal = model))
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be one number between 0 and 1.")
  }
  p_adjust_method <- match.arg(p_adjust_method, stats::p.adjust.methods)

  conditional_formula <- stats::formula(model, component = "cond")
  fixed_formula <- stats::formula(
    model,
    component = "cond",
    fixed.only = TRUE
  )

  model_terms <- stats::terms(fixed_formula)
  term_labels <- attr(model_terms, "term.labels")
  term_order <- attr(model_terms, "order")
  main_effects <- term_labels[term_order == 1L]
  interactions <- term_labels[term_order > 1L]

  interaction_parts <- strsplit(interactions, ":", fixed = TRUE)
  interacting_effects <- main_effects[vapply(
    main_effects,
    function(term) any(vapply(
      interaction_parts,
      function(parts) term %in% parts,
      logical(1)
    )),
    logical(1)
  )]
  interacting_effects <- setdiff(interacting_effects, exclude)

  if (length(interacting_effects) == 0L) {
    stop("No eligible fixed effects with interactions were found.")
  }

  dropped_terms <- setNames(
    lapply(interacting_effects, function(term) {
      related_interactions <- interactions[vapply(
        interaction_parts,
        function(parts) term %in% parts,
        logical(1)
      )]
      c(term, related_interactions)
    }),
    interacting_effects
  )

  reduced_models <- setNames(
    vector("list", length(dropped_terms)),
    names(dropped_terms)
  )
  comparisons <- setNames(
    vector("list", length(dropped_terms)),
    names(dropped_terms)
  )
  summary_rows <- setNames(
    vector("list", length(dropped_terms)),
    names(dropped_terms)
  )

  for (term in names(dropped_terms)) {
    terms_to_drop <- dropped_terms[[term]]
    reduction_formula <- stats::as.formula(
      paste(
        ". ~ .",
        paste0("- (", terms_to_drop, ")", collapse = " ")
      ),
      env = environment(conditional_formula)
    )

    reduced_model <- update(model, reduction_formula)
    if (stats::nobs(reduced_model) != stats::nobs(model)) {
      stop(
        "Removing `", term, "` changed the number of observations. ",
        "Fit the maximal model to an explicit complete-case dataset."
      )
    }

    comparison <- stats::anova(reduced_model, model, test = "Chisq")
    reduced_models[[term]] <- reduced_model
    comparisons[[term]] <- comparison

    p_column <- grep("^Pr\\(", names(comparison), value = TRUE)
    if (length(p_column) != 1L) {
      stop("Could not identify the likelihood-ratio P-value column.")
    }
    comparison_row <- nrow(comparison)
    summary_rows[[term]] <- tibble::tibble(
      predictor = term,
      df = as.numeric(comparison[["Chi Df"]][comparison_row]),
      likelihood_ratio = as.numeric(comparison[["Chisq"]][comparison_row]),
      p_value = as.numeric(comparison[[p_column]][comparison_row])
    )

    cat("\nPredictor block:", term, "\n")
    cat("Dropped terms:", paste(terms_to_drop, collapse = ", "), "\n")
    print(comparison)
  }

  summary_table <- dplyr::bind_rows(summary_rows) |>
    dplyr::mutate(
      p_adjusted = stats::p.adjust(p_value, method = p_adjust_method),
      significantly_improves_fit = p_adjusted < alpha
    )

  significant <- summary_table$predictor[
    summary_table$significantly_improves_fit %in% TRUE
  ]
  not_significant <- summary_table$predictor[
    summary_table$significantly_improves_fit %in% FALSE
  ]
  not_classified <- summary_table$predictor[
    is.na(summary_table$significantly_improves_fit)
  ]
  collapse_or_none <- function(x) {
    if (length(x) == 0L) "none" else paste(x, collapse = ", ")
  }

  cat(
    "\nSummary (alpha = ", alpha,
    "; P adjustment = ", p_adjust_method, ")\n",
    sep = ""
  )
  cat(
    "Significantly improve model fit: ",
    collapse_or_none(significant),
    "\n",
    sep = ""
  )
  cat(
    "Do not significantly improve model fit: ",
    collapse_or_none(not_significant),
    "\n",
    sep = ""
  )
  if (length(not_classified) > 0L) {
    cat(
      "Could not be classified: ",
      collapse_or_none(not_classified),
      "\n",
      sep = ""
    )
  }

  invisible(list(
    dropped_terms = dropped_terms,
    reduced_models = reduced_models,
    comparisons = comparisons,
    summary = summary_table
  ))
}
# Test residual temporal autocorrelation separately within each plot.
# Plots require at least three unique observed years. The result contains all
# plot-level tests, Holm-adjusted P values, and a descriptive Fisher combination
# across plots. `model_data` must follow the fitted-model row order.
# Details: Orders observations by year within site-by-plot groups, omits groups
# unsuitable for testing, and reports Holm and descriptive Fisher summaries.
# Returns: A list with the overall summary and one row per plot-level test.
test_within_plot_temporal_autocorrelation <- function(
    simulated_residuals,
    model,
    seed = 123L,
    model_data = NULL) {
  if (is.null(model_data)) {
    model_data <- stats::model.frame(model)
  }

  if (nrow(model_data) != length(simulated_residuals$scaledResiduals)) {
    stop(
      "`model_data` rows must match the simulated residuals."
    )
  }

  required_columns <- c("site", "plot_id", "year")

  if (!all(required_columns %in% names(model_data))) {
    return(list(
      summary = tibble::tibble(
        plots_tested = 0L,
        minimum_p = NA_real_,
        minimum_holm_p = NA_real_,
        fisher_p = NA_real_,
        note = "site, plot_id, or year was absent from the model frame"
      ),
      by_plot = tibble::tibble()
    ))
  }

  plot_group <- interaction(
    model_data$site,
    model_data$plot_id,
    drop = TRUE
  )
  numeric_year <- suppressWarnings(as.numeric(as.character(model_data$year)))
  plot_indices <- split(seq_len(nrow(model_data)), plot_group)

  temporal_rows <- lapply(seq_along(plot_indices), function(i) {
    plot_name <- names(plot_indices)[[i]]
    index <- plot_indices[[i]]
    index <- index[order(numeric_year[index])]

    if (length(index) < 3L || anyNA(numeric_year[index]) ||
        anyDuplicated(numeric_year[index])) {
      return(tibble::tibble(
        plot = plot_name,
        observations = length(index),
        statistic = NA_real_,
        p_value = NA_real_,
        note = "fewer than three observations, missing years, or duplicate years"
      ))
    }

    test_result <- tryCatch({
      plot_residuals <- DHARMa::recalculateResiduals(
        simulated_residuals,
        sel = index,
        seed = seed + i - 1L
      )
      DHARMa::testTemporalAutocorrelation(
        plot_residuals,
        time = numeric_year[index],
        plot = FALSE
      )
    }, error = function(e) e)

    if (inherits(test_result, "error")) {
      tibble::tibble(
        plot = plot_name,
        observations = length(index),
        statistic = NA_real_,
        p_value = NA_real_,
        note = conditionMessage(test_result)
      )
    } else {
      tibble::tibble(
        plot = plot_name,
        observations = length(index),
        statistic = as.numeric(test_result$statistic[[1]]),
        p_value = as.numeric(test_result$p.value[[1]]),
        note = NA_character_
      )
    }
  })

  by_plot <- dplyr::bind_rows(temporal_rows) |>
    dplyr::mutate(
      p_holm = stats::p.adjust(p_value, method = "holm")
    )
  valid_p <- by_plot$p_value[is.finite(by_plot$p_value)]

  temporal_summary <- tibble::tibble(
    plots_tested = length(valid_p),
    minimum_p = if (length(valid_p) > 0L) min(valid_p) else NA_real_,
    minimum_holm_p = if (length(valid_p) > 0L) {
      min(by_plot$p_holm, na.rm = TRUE)
    } else {
      NA_real_
    },
    fisher_p = if (length(valid_p) > 0L) {
      stats::pchisq(
        -2 * sum(log(pmax(valid_p, .Machine$double.xmin))),
        df = 2 * length(valid_p),
        lower.tail = FALSE
      )
    } else {
      NA_real_
    },
    note = NA_character_
  )

  list(summary = temporal_summary, by_plot = by_plot)
}


# Evaluate prespecified nuisance-structure candidates fitted to the same data.
# The function combines information criteria, DHARMa checks, within-plot
# temporal tests, singularity checks, variance components, and fixed effects.
# Returns one audit list while leaving model selection to the analyst.
# Details: Collects fit criteria, residual and temporal checks, singularity,
# variance components, and coefficients without making the selection itself.
# Returns: One named audit list containing candidates and all comparison tables.
evaluate_maximal_structure_candidates <- function(
    models,
    seed,
    n_sim = 1000L,
    temporal_data = NULL) {
  comparison <- compare_candidate_models(models)
  diagnostics <- run_candidate_diagnostics(
    models,
    n_sim = n_sim,
    seed = seed
  )

  temporal <- stats::setNames(
    lapply(seq_along(models), function(i) {
      model_name <- names(models)[[i]]
      residuals <- diagnostics$residuals[[model_name]]

      if (is.null(residuals)) {
        return(list(
          summary = tibble::tibble(
            plots_tested = NA_integer_,
            minimum_p = NA_real_,
            minimum_holm_p = NA_real_,
            fisher_p = NA_real_,
            note = "Model did not converge with a positive-definite Hessian."
          ),
          by_plot = tibble::tibble()
        ))
      }

      test_within_plot_temporal_autocorrelation(
        residuals,
        models[[model_name]],
        seed = seed + 1000L + i - 1L,
        model_data = temporal_data
      )
    }),
    names(models)
  )

  temporal_summary <- dplyr::bind_rows(
    lapply(names(temporal), function(model_name) {
      dplyr::mutate(
        temporal[[model_name]]$summary,
        model = model_name,
        .before = 1
      )
    })
  )

  singularity <- stats::setNames(
    lapply(models, function(model) {
      tryCatch(
        performance::check_singularity(model, check = "terms"),
        error = function(e) list(error = conditionMessage(e))
      )
    }),
    names(models)
  )

  fixed_terms <- unique(unlist(lapply(models, function(model) {
    names(glmmTMB::fixef(model)$cond)
  })))

  list(
    candidates = models,
    comparison = comparison,
    diagnostics = diagnostics,
    temporal = temporal,
    temporal_summary = temporal_summary,
    singularity = singularity,
    variance_components = lapply(models, VarCorr),
    fixed_effects = summarize_fixed_terms(models, fixed_terms)
  )
}


# Validate one selected maximal model before inferential comparisons.
# Records numerical fit status, fixed-gradient size, AIC/AICc, singularity,
# DHARMa tests, collinearity, and within-plot temporal autocorrelation.
# Stops if convergence or the Hessian check fails; otherwise returns the audit.
# Details: Stops before inference if convergence or Hessian checks fail; valid
# fits receive DHARMa, collinearity, singularity, gradient, and temporal checks.
# Returns: An invisible audit list and concise console summaries.
validate_maximal_model <- function(
    model,
    model_name,
    n_sim = 1000L,
    seed = 123L,
    temporal_data = NULL) {
  model_loglik <- suppressWarnings(stats::logLik(model))
  parameters <- attr(model_loglik, "df")
  n <- stats::nobs(model)
  model_aic <- if (model_is_valid(model)) stats::AIC(model) else NA_real_
  model_aicc <- if (is.finite(model_aic) && n > parameters + 1L) {
    model_aic + (2 * parameters * (parameters + 1L)) /
      (n - parameters - 1L)
  } else {
    NA_real_
  }
  fixed_gradient <- model$sdr$gradient.fixed
  converged <- isTRUE(model$fit$convergence == 0L)
  pd_hessian <- isTRUE(model$sdr$pdHess)
  optimizer_message <- if (length(model$fit$message) > 0L) {
    paste(model$fit$message, collapse = "; ")
  } else {
    NA_character_
  }

  singular <- tryCatch(
    as.logical(performance::check_singularity(model))[[1]],
    error = function(e) NA
  )

  fit_status <- tibble::tibble(
    model = model_name,
    converged = converged,
    pd_hessian = pd_hessian,
    singular = singular,
    n = n,
    parameters = parameters,
    max_abs_fixed_gradient = if (length(fixed_gradient) > 0L) {
      max(abs(fixed_gradient))
    } else {
      NA_real_
    },
    AIC = model_aic,
    AICc = model_aicc,
    optimizer_message = optimizer_message
  )

  cat("\nMaximal-model fit status:\n")
  print(fit_status)

  if (!model_is_valid(model)) {
    stop(
      "`", model_name,
      "` did not converge with a positive-definite Hessian; ",
      "inferential tests were not run."
    )
  }

  diagnostic_set <- stats::setNames(list(model), model_name)
  dharma <- run_candidate_diagnostics(
    diagnostic_set,
    n_sim = n_sim,
    seed = seed
  )

  collinearity <- tryCatch(
    as.data.frame(performance::check_collinearity(model)),
    error = function(e) data.frame(error = conditionMessage(e))
  )

  temporal <- test_within_plot_temporal_autocorrelation(
    dharma$residuals[[model_name]],
    model,
    seed = seed + 1000L,
    model_data = temporal_data
  )

  cat("\nDHARMa diagnostic summary:\n")
  print(dharma$summary)
  cat("\nCollinearity summary:\n")
  print(collinearity)
  cat("\nWithin-plot temporal-autocorrelation summary:\n")
  print(temporal$summary)

  invisible(list(
    fit_status = fit_status,
    dharma = dharma,
    collinearity = collinearity,
    temporal = temporal
  ))
}

# ------------------------------------------------------------------------------
# Maximal-SEM inference, standardization, and component-R2 helpers
# ------------------------------------------------------------------------------
# These functions are used by Maximal_SEM.R and Parsimonious_SEM.R. Internal
# helpers begin with `.sem_`; user-facing functions assemble tests,
# standardizations, and R2.

# Check identifiers before they are used to select models, data, or variables.
# Details: Prevents missing, empty, or vector-valued names from reaching model,
# dataset, node, or variable lookups.
# Returns: The validated value invisibly.
.sem_assert_scalar_character <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop("`", name, "` must be one non-empty character value.")
  }
  invisible(x)
}

# Test a fitted term at declared moderator reference values using the exact
# fixed-effect covariance matrix. For interaction terms, this reduces to the
# fitted interaction-coefficient Wald test.
# Details: Forms the appropriate fixed-effect gradient, incorporates interacting
# moderators at declared reference values, and uses the full covariance matrix.
# Returns: A one-row tibble with estimate, uncertainty, z test, and conditioning.
sem_term_wald_test <- function(model, term, term_type, reference_values) {
  coefficients <- glmmTMB::fixef(model)$cond
  covariance <- as.matrix(stats::vcov(model)$cond)
  if (!term %in% names(coefficients)) {
    stop("Term `", term, "` is absent from the fitted component model.")
  }
  covariance <- covariance[
    names(coefficients), names(coefficients), drop = FALSE
  ]
  gradient <- stats::setNames(numeric(length(coefficients)), names(coefficients))
  gradient[[term]] <- 1
  conditioning <- character()

  if (identical(term_type, "direct_path")) {
    interaction_terms <- names(coefficients)[vapply(
      strsplit(names(coefficients), ":", fixed = TRUE),
      function(parts) length(parts) > 1L && term %in% parts,
      logical(1)
    )]
    for (interaction_term in interaction_terms) {
      moderators <- setdiff(
        strsplit(interaction_term, ":", fixed = TRUE)[[1]],
        term
      )
      if (!all(moderators %in% names(reference_values))) {
        stop("Missing reference value for `", interaction_term, "`.")
      }
      moderator_values <- reference_values[moderators]
      gradient[[interaction_term]] <- prod(moderator_values)
      conditioning <- c(
        conditioning,
        paste0(moderators, "=", format(moderator_values, digits = 8))
      )
    }
  }

  estimate <- sum(gradient * coefficients)
  standard_error <- sqrt(drop(gradient %*% covariance %*% gradient))
  statistic <- estimate / standard_error
  tibble::tibble(
    inference_estimand = dplyr::if_else(
      term_type == "direct_path",
      "reference_link_slope",
      "interaction_coefficient"
    ),
    inference_method = "Wald z from exact fixed-effect covariance",
    reference_conditioning = if (length(conditioning) == 0L) {
      if (term_type == "direct_path") "no interacting moderators" else NA_character_
    } else {
      paste(unique(conditioning), collapse = "; ")
    },
    inference_estimate = estimate,
    inference_standard_error = standard_error,
    inference_conf_low = estimate - stats::qnorm(0.975) * standard_error,
    inference_conf_high = estimate + stats::qnorm(0.975) * standard_error,
    inference_test_statistic = statistic,
    inference_p_value = 2 * stats::pnorm(-abs(statistic))
  )
}

# Apply declared conditioning values while preserving every evaluation row.
# Details: Replaces named variables with scalar or row-specific conditioning
# values while preventing the focal predictor from also being fixed.
# Returns: A modified copy of the evaluation data with row order preserved.
.sem_apply_at <- function(data, at = list(), predictor = NULL) {
  if (is.null(at)) at <- list()
  if (!is.list(at) || (length(at) > 0L && is.null(names(at)))) {
    stop("`at` must be a named list of conditioning values.")
  }
  if (!is.null(predictor) && predictor %in% names(at)) {
    stop("The focal predictor cannot also be fixed through `at`.")
  }

  out <- data
  for (variable_name in names(at)) {
    if (!variable_name %in% names(out)) {
      stop("Conditioning variable `", variable_name, "` is absent from `data`.")
    }
    value <- at[[variable_name]]
    if (!length(value) %in% c(1L, nrow(out))) {
      stop(
        "Conditioning value for `", variable_name,
        "` must have length one or nrow(data)."
      )
    }
    out[[variable_name]] <- if (length(value) == 1L) {
      rep(value, nrow(out))
    } else {
      value
    }
  }
  out
}

# Calculate a weighted mean after validating values and non-negative weights.
# Details: Uses equal weights when none are supplied and rejects missing,
# infinite, negative, mismatched, or zero-total weights.
# Returns: One finite weighted mean.
.sem_weighted_mean <- function(x, weights = NULL) {
  if (is.null(weights)) weights <- rep(1, length(x))
  if (!is.numeric(weights) || length(weights) != length(x)) {
    stop("`weights` must be NULL or a numeric vector matching the data rows.")
  }
  keep <- is.finite(x) & is.finite(weights) & weights >= 0
  if (!all(keep)) {
    stop("Effects or averaging weights contain missing, infinite, or negative values.")
  }
  if (sum(weights) <= 0) stop("Averaging weights must have a positive sum.")
  stats::weighted.mean(x, weights)
}

# Fixed-effect predictions for the model classes currently used by the SEMs.
# A custom function can support another class; it must accept arguments
# (model, newdata, type, include_random) and return one numeric value per row.
# Details: Standardizes fixed- versus conditional-prediction behavior across
# supported model classes and permits an explicit custom prediction adapter.
# Returns: One finite numeric prediction for every evaluation row.
.sem_predict <- function(
    model,
    newdata,
    type = c("response", "link"),
    include_random = FALSE,
    predict_function = NULL) {
  type <- match.arg(type)

  if (!is.null(predict_function)) {
    predicted <- predict_function(
      model = model,
      newdata = newdata,
      type = type,
      include_random = include_random
    )
  } else if (inherits(model, "glmmTMB")) {
    predicted <- stats::predict(
      model,
      newdata = newdata,
      type = type,
      re.form = if (isTRUE(include_random)) NULL else NA,
      allow.new.levels = TRUE
    )
  } else if (inherits(model, "merMod")) {
    predicted <- stats::predict(
      model,
      newdata = newdata,
      type = type,
      re.form = if (isTRUE(include_random)) NULL else NA,
      allow.new.levels = TRUE
    )
  } else if (inherits(model, "glm")) {
    predicted <- stats::predict(model, newdata = newdata, type = type)
  } else if (inherits(model, "lm")) {
    if (type == "link") type <- "response"
    predicted <- stats::predict(model, newdata = newdata, type = type)
  } else {
    stop(
      "Unsupported model class. Supply `predict_function` so fixed- and ",
      "random-effect prediction semantics are explicit."
    )
  }

  predicted <- as.numeric(predicted)
  if (length(predicted) != nrow(newdata) || any(!is.finite(predicted))) {
    stop("Predictions must be finite and have one value per evaluation row.")
  }
  predicted
}

# Numerical derivative of either the exact response-scale conditional mean or
# the linear predictor. Near declared bounds this automatically becomes a
# one-sided derivative. Otherwise it is a central difference.
# Details: Perturbs the focal predictor within declared bounds, predicts at both
# values, and averages row-level finite-difference slopes with optional weights.
# Returns: A list containing the mean slope, row slopes, step, and row count.
.sem_average_slope <- function(
    model,
    data,
    predictor,
    type,
    scale_reference,
    at = list(),
    weights = NULL,
    step = NULL,
    bounds = c(-Inf, Inf),
    include_random = FALSE,
    predict_function = NULL,
    evaluation_function = NULL) {
  .sem_assert_scalar_character(predictor, "predictor")
  if (!predictor %in% names(data) || !is.numeric(data[[predictor]])) {
    stop("`predictor` must name a numeric column in `data`.")
  }
  if (!is.numeric(scale_reference) || length(scale_reference) != 1L ||
      !is.finite(scale_reference) || scale_reference <= 0) {
    stop("`scale_reference` must be one positive finite number.")
  }
  if (!is.numeric(bounds) || length(bounds) != 2L ||
      is.na(bounds[1]) || is.na(bounds[2]) || bounds[1] >= bounds[2]) {
    stop("`bounds` must contain an ordered lower and upper bound.")
  }

  evaluation_data <- .sem_apply_at(data, at = at, predictor = predictor)
  x <- evaluation_data[[predictor]]
  if (any(!is.finite(x))) stop("The focal predictor contains non-finite values.")

  if (is.null(step)) {
    step <- .Machine$double.eps^(1 / 3) * max(scale_reference, 1)
  }
  if (!is.numeric(step) || length(step) != 1L || !is.finite(step) || step <= 0) {
    stop("`step` must be one positive finite number.")
  }

  x_low <- pmax(x - step, bounds[1])
  x_high <- pmin(x + step, bounds[2])
  denominator <- x_high - x_low
  if (any(!is.finite(denominator)) || any(denominator <= 0)) {
    stop("The derivative step and predictor bounds do not define valid contrasts.")
  }

  low_data <- evaluation_data
  high_data <- evaluation_data
  low_data[[predictor]] <- x_low
  high_data[[predictor]] <- x_high
  if (!is.null(evaluation_function)) {
    evaluation <- evaluation_function(
      evaluation_data, predictor, low_data, high_data, mode = "slope"
    )
    low_data <- evaluation$low
    high_data <- evaluation$high
    if (nrow(low_data) != nrow(data) || nrow(high_data) != nrow(data)) {
      stop("Slope evaluation must preserve every native-sample row.")
    }
    denominator <- high_data[[predictor]] - low_data[[predictor]]
    if (any(!is.finite(denominator)) || any(denominator <= 0)) {
      stop("Adjusted slope contrasts must have positive finite widths.")
    }
    if (!is.null(evaluation$step)) step <- evaluation$step
  }

  low_prediction <- .sem_predict(
    model, low_data, type, include_random, predict_function
  )
  high_prediction <- .sem_predict(
    model, high_data, type, include_random, predict_function
  )
  individual_slope <- (high_prediction - low_prediction) / denominator

  list(
    average_slope = .sem_weighted_mean(individual_slope, weights),
    individual_slope = individual_slope,
    step = step,
    n_evaluation = nrow(evaluation_data)
  )
}

# Standardize an average response-scale partial effect by the node SDs.
# Details: Calculates the response-scale average partial effect and multiplies
# it by the source-node SD divided by the destination-node SD.
# Returns: A one-row data frame with the derivative and standardized effect.
standardize_response_path <- function(
    model,
    data,
    predictor,
    sd_x_reference,
    sd_y_reference,
    at = list(),
    average_weights = NULL,
    derivative_step = NULL,
    predictor_bounds = c(-Inf, Inf),
    include_random = FALSE,
    predict_function = NULL,
    evaluation_function = NULL) {
  if (!is.numeric(sd_y_reference) || length(sd_y_reference) != 1L ||
      !is.finite(sd_y_reference) || sd_y_reference <= 0) {
    stop("`sd_y_reference` must be one positive finite number.")
  }

  slope <- .sem_average_slope(
    model = model,
    data = data,
    predictor = predictor,
    type = "response",
    scale_reference = sd_x_reference,
    at = at,
    weights = average_weights,
    step = derivative_step,
    bounds = predictor_bounds,
    include_random = include_random,
    predict_function = predict_function,
    evaluation_function = evaluation_function
  )

  data.frame(
    n_evaluation = slope$n_evaluation,
    derivative_step = slope$step,
    average_partial_effect = slope$average_slope,
    sd_x_reference = sd_x_reference,
    sd_y_reference = sd_y_reference,
    standardized_response_effect =
      slope$average_slope * sd_x_reference / sd_y_reference,
    stringsAsFactors = FALSE
  )
}

# Calculate Menard observed-empirical standardization on the link scale.
# Details: Combines the average link-scale slope with the fitted component's
# observed-empirical R2 and linear-predictor variance.
# Returns: A one-row data frame with OE quantities and the standardized effect.
standardize_oe_path <- function(
    model,
    data,
    predictor,
    response,
    sd_x_reference,
    at = list(),
    average_weights = NULL,
    derivative_step = NULL,
    predictor_bounds = c(-Inf, Inf),
    include_random = FALSE,
    predict_function = NULL,
    evaluation_function = NULL) {
  .sem_assert_scalar_character(response, "response")
  if (!response %in% names(data) || !is.numeric(data[[response]])) {
    stop("`response` must name a numeric column in `data`.")
  }

  # Menard's R2 and implied response SD describe the complete fitted component
  # on its native sample. Conditioning values affect the focal link-scale slope
  # but do not replace the data used to calculate this full-model R2.
  eta_hat <- .sem_predict(
    model, data, "link", include_random, predict_function
  )
  mu_hat <- .sem_predict(
    model, data, "response", include_random, predict_function
  )
  y <- data[[response]]
  keep <- is.finite(y) & is.finite(eta_hat) & is.finite(mu_hat)
  if (sum(keep) < 3L || stats::sd(y[keep]) <= 0 || stats::sd(mu_hat[keep]) <= 0) {
    stop("The native response and fitted values must vary for OE standardization.")
  }

  r2_oe <- stats::cor(y[keep], mu_hat[keep])^2
  variance_eta <- stats::var(eta_hat[keep])
  if (!is.finite(r2_oe) || r2_oe <= 0 ||
      !is.finite(variance_eta) || variance_eta <= 0) {
    stop("The OE R2 and linear-predictor variance must be positive and finite.")
  }
  sd_y_oe <- sqrt(variance_eta / r2_oe)

  link_slope <- .sem_average_slope(
    model = model,
    data = data,
    predictor = predictor,
    type = "link",
    scale_reference = sd_x_reference,
    at = at,
    weights = average_weights,
    step = derivative_step,
    bounds = predictor_bounds,
    include_random = include_random,
    predict_function = predict_function,
    evaluation_function = evaluation_function
  )

  data.frame(
    oe_link_slope = link_slope$average_slope,
    oe_r2 = r2_oe,
    oe_variance_eta = variance_eta,
    oe_sd_y = sd_y_oe,
    standardized_oe_effect =
      link_slope$average_slope * sd_x_reference / sd_y_oe,
    stringsAsFactors = FALSE
  )
}

# Contrast predictions across a declared biological range and scale by the
# response range.
# Details: Contrasts mean predictions at declared low and high predictor values,
# scales by the destination range, and checks observed predictor support.
# Returns: A one-row data frame with contrast, support flags, and effect size.
standardize_relevant_range_path <- function(
    model,
    data,
    predictor,
    x_range_model,
    y_range_model,
    at = list(),
    average_weights = NULL,
    include_random = FALSE,
    predict_function = NULL,
    support_tolerance = sqrt(.Machine$double.eps),
    evaluation_function = NULL) {
  if (!is.numeric(x_range_model) || length(x_range_model) != 2L ||
      any(!is.finite(x_range_model)) || x_range_model[1] >= x_range_model[2]) {
    stop("`x_range_model` must contain finite increasing model-scale endpoints.")
  }
  if (!is.numeric(y_range_model) || length(y_range_model) != 2L ||
      any(!is.finite(y_range_model)) || y_range_model[1] >= y_range_model[2]) {
    stop("`y_range_model` must contain finite increasing model-scale endpoints.")
  }
  if (!predictor %in% names(data) || !is.numeric(data[[predictor]])) {
    stop("`predictor` must name a numeric column in `data`.")
  }

  evaluation_data <- .sem_apply_at(data, at = at, predictor = predictor)
  low_data <- evaluation_data
  high_data <- evaluation_data
  low_data[[predictor]] <- x_range_model[1]
  high_data[[predictor]] <- x_range_model[2]
  if (!is.null(evaluation_function)) {
    evaluation <- evaluation_function(
      evaluation_data, predictor, low_data, high_data, mode = "range"
    )
    low_data <- evaluation$low
    high_data <- evaluation$high
    if (nrow(low_data) == 0L || nrow(low_data) != nrow(high_data)) {
      stop("No common feasible rows for the declared range contrast.")
    }
    if (!is.null(average_weights)) average_weights <- average_weights[evaluation$rows]
  }

  mu_low <- .sem_predict(
    model, low_data, "response", include_random, predict_function
  )
  mu_high <- .sem_predict(
    model, high_data, "response", include_random, predict_function
  )
  individual_contrast <- mu_high - mu_low
  response_contrast <- .sem_weighted_mean(individual_contrast, average_weights)
  response_width <- diff(y_range_model)

  observed_x <- data[[predictor]]
  observed_x <- observed_x[is.finite(observed_x)]
  observed_range <- range(observed_x)
  low_supported <- x_range_model[1] >= observed_range[1] - support_tolerance
  high_supported <- x_range_model[2] <= observed_range[2] + support_tolerance

  result <- data.frame(
    rr_x_low_model = x_range_model[1],
    rr_x_high_model = x_range_model[2],
    rr_y_low_model = y_range_model[1],
    rr_y_high_model = y_range_model[2],
    rr_response_contrast = response_contrast,
    standardized_relevant_range_effect = response_contrast / response_width,
    rr_x_low_supported = low_supported,
    rr_x_high_supported = high_supported,
    rr_x_range_supported = low_supported && high_supported,
    stringsAsFactors = FALSE
  )
  if (!is.null(evaluation_function)) {
    result$rr_n_evaluation <- nrow(low_data)
    result$rr_n_excluded <- nrow(data) - nrow(low_data)
  }
  result
}

# Retrieve one node definition from the reference-scale registry.
# Details: Verifies the node registry has the required scale columns and that
# the requested biological node is represented exactly once.
# Returns: The single matching node-reference row.
.sem_node_row <- function(node_reference, node) {
  if (!is.data.frame(node_reference) ||
      !all(c("node", "model_variable", "sd_model") %in% names(node_reference))) {
    stop("`node_reference` lacks required node-scale columns.")
  }
  row <- node_reference[node_reference$node == node, , drop = FALSE]
  if (nrow(row) != 1L) {
    stop("Node `", node, "` must occur exactly once in `node_reference`.")
  }
  row
}

# Convert biological range endpoints to the fitted model scale.
# Details: Applies the declared square-root, cap, log2, identity, or z-scale
# transformation and validates endpoint and scaling metadata.
# Returns: One transformed model-scale endpoint, or NA for an undeclared value.
.sem_range_to_model_scale <- function(
    value,
    transformation,
    scaling_center = NA_real_,
    scaling_divisor = NA_real_) {
  if (is.na(value)) return(NA_real_)
  if (identical(transformation, "square root")) {
    if (value < 0) stop("Square-root range endpoints cannot be negative.")
    return(sqrt(value))
  }
  if (identical(transformation, "cap at 10")) {
    if (value < 0) stop("Richness range endpoints cannot be negative.")
    return(pmin(value, 10))
  }
  if (identical(transformation, "cap at 16; log2; z-scale")) {
    if (value < 1) stop("Tree-richness range endpoints cannot be below one.")
    value <- log2(pmin(value, 16))
  } else if (identical(transformation, "square root; z-scale")) {
    if (value < 0) stop("Square-root range endpoints cannot be negative.")
    value <- sqrt(value)
  } else if (transformation %in% c("identity", "identity (proportion)")) {
    if (identical(transformation, "identity (proportion)") &&
        (value < 0 || value > 1)) {
      stop("Proportion range endpoints must lie between zero and one.")
    }
    return(value)
  } else if (!identical(transformation, "z-scale")) {
    stop("Unsupported relevant-range transformation: `", transformation, "`.")
  }
  if (!is.finite(scaling_center) || !is.finite(scaling_divisor) ||
      scaling_divisor <= 0) {
    stop("A z-scaled node lacks valid center and divisor metadata.")
  }
  (value - scaling_center) / scaling_divisor
}

# Validate a compact biological-scale range registry, convert its endpoints to
# each fitted node scale, and report whether every required node was declared.
# Details: Requires one range declaration per node, transforms valid biological
# endpoints, checks observed support, and identifies complete required ranges.
# Returns: A list with the enriched node registry, required nodes, and ready flag.
prepare_sem_relevant_ranges <- function(
    node_reference,
    range_spec,
    required_nodes = node_reference$node) {
  node_columns <- c(
    "node", "transformation", "observed_low_biological",
    "observed_high_biological", "scaling_center", "scaling_divisor"
  )
  range_columns <- c(
    "node", "rr_low_biological", "rr_high_biological",
    "rr_rationale"
  )
  if (!all(node_columns %in% names(node_reference))) {
    stop("`node_reference` lacks relevant-range metadata columns.")
  }
  if (!all(range_columns %in% names(range_spec))) {
    stop("`range_spec` lacks required declaration columns.")
  }
  if (anyDuplicated(range_spec$node) ||
      !setequal(range_spec$node, node_reference$node)) {
    stop("`range_spec` must contain exactly one row per SEM node.")
  }
  if (!all(required_nodes %in% node_reference$node)) {
    stop("`required_nodes` contains an unknown SEM node.")
  }

  declared <- with(
    range_spec,
    !is.na(rr_low_biological) & !is.na(rr_high_biological) &
      !is.na(rr_rationale) & nzchar(rr_rationale)
  )
  any_entry <- with(
    range_spec,
    !is.na(rr_low_biological) | !is.na(rr_high_biological) |
      (!is.na(rr_rationale) & nzchar(rr_rationale))
  )
  if (any(any_entry & !declared)) {
    stop(
      "Each range needs both endpoints and a rationale; ",
      "otherwise leave all three entries as NA."
    )
  }
  if (any(
    declared & range_spec$rr_low_biological >= range_spec$rr_high_biological
  )) {
    stop("Every declared relevant range must have increasing endpoints.")
  }

  output <- dplyr::left_join(node_reference, range_spec, by = "node")
  output$rr_low_model <- mapply(
    .sem_range_to_model_scale,
    output$rr_low_biological,
    output$transformation,
    output$scaling_center,
    output$scaling_divisor,
    USE.NAMES = FALSE
  )
  output$rr_high_model <- mapply(
    .sem_range_to_model_scale,
    output$rr_high_biological,
    output$transformation,
    output$scaling_center,
    output$scaling_divisor,
    USE.NAMES = FALSE
  )
  output$rr_width_biological <-
    output$rr_high_biological - output$rr_low_biological
  output$rr_width_model <- output$rr_high_model - output$rr_low_model
  output$rr_supported <- ifelse(
    declared[match(output$node, range_spec$node)],
    output$rr_low_biological >= output$observed_low_biological &
      output$rr_high_biological <= output$observed_high_biological,
    NA
  )

  ready <- all(declared[match(required_nodes, range_spec$node)])
  if (ready) {
    required <- output[match(required_nodes, output$node), , drop = FALSE]
    if (any(!is.finite(required$rr_low_model)) ||
        any(!is.finite(required$rr_high_model)) ||
        any(required$rr_low_model >= required$rr_high_model)) {
      stop("At least one required range is invalid on its fitted model scale.")
    }
  }

  list(
    node_reference = output,
    required_nodes = required_nodes,
    ready = ready
  )
}
# Calculate all declared summaries for one direct edge. The primary and OE
# results are always calculated; the relevant-range result is calculated only
# when both nodes have complete model-scale endpoints.
# Details: Validates that the model, native dataset, and node registry agree,
# then combines response, OE, and optional relevant-range standardizations.
# Returns: One data-frame row containing metadata and all path effect summaries.
standardize_sem_path <- function(
    model,
    data,
    source_node,
    destination_node,
    node_reference,
    predictor = NULL,
    response = NULL,
    at = list(),
    average_weights = NULL,
    derivative_step = NULL,
    predictor_bounds = c(-Inf, Inf),
    include_relevant_range = TRUE,
    include_random = FALSE,
    predict_function = NULL,
    evaluation_function = NULL) {
  source_row <- .sem_node_row(node_reference, source_node)
  destination_row <- .sem_node_row(node_reference, destination_node)
  if (is.null(predictor)) predictor <- source_row$model_variable[[1]]
  if (is.null(response)) response <- destination_row$model_variable[[1]]

  model_response <- all.vars(stats::formula(model)[[2L]])
  response_expression <- stats::formula(model)[[2L]]
  count_response <- is.call(response_expression) &&
    identical(response_expression[[1L]], as.name("cbind")) &&
    length(response_expression) == 3L &&
    identical(response_expression[[2L]], as.name(response))
  if (!count_response &&
      (length(model_response) != 1L || !identical(model_response, response))) {
    stop(
      "The declared scalar response `", response,
      "` does not match the fitted model response."
    )
  }
  model_n <- tryCatch(stats::nobs(model), error = function(error) NA_real_)
  if (is.finite(model_n) && model_n != nrow(data)) {
    stop(
      "`data` has ", nrow(data), " rows but the fitted model used ", model_n,
      "; supply the exact native component sample."
    )
  }

  primary <- standardize_response_path(
    model, data, predictor,
    sd_x_reference = source_row$sd_model[[1]],
    sd_y_reference = destination_row$sd_model[[1]],
    at = at,
    average_weights = average_weights,
    derivative_step = derivative_step,
    predictor_bounds = predictor_bounds,
    include_random = include_random,
    predict_function = predict_function,
    evaluation_function = evaluation_function
  )
  oe <- standardize_oe_path(
    model, data, predictor, response,
    sd_x_reference = source_row$sd_model[[1]],
    at = at,
    average_weights = average_weights,
    derivative_step = primary$derivative_step[[1]],
    predictor_bounds = predictor_bounds,
    include_random = include_random,
    predict_function = predict_function,
    evaluation_function = evaluation_function
  )

  range_columns <- c("rr_low_model", "rr_high_model")
  range_available <- isTRUE(include_relevant_range) &&
    all(range_columns %in% names(source_row)) &&
    all(range_columns %in% names(destination_row)) &&
    all(is.finite(unlist(source_row[1, range_columns], use.names = FALSE))) &&
    all(is.finite(unlist(destination_row[1, range_columns], use.names = FALSE)))

  if (range_available) {
    relevant_range <- standardize_relevant_range_path(
      model = model,
      data = data,
      predictor = predictor,
      x_range_model = unlist(source_row[1, range_columns], use.names = FALSE),
      y_range_model = unlist(destination_row[1, range_columns], use.names = FALSE),
      at = at,
      average_weights = average_weights,
      include_random = include_random,
      predict_function = predict_function,
      evaluation_function = evaluation_function
    )
  } else {
    relevant_range <- data.frame(
      rr_x_low_model = NA_real_,
      rr_x_high_model = NA_real_,
      rr_y_low_model = NA_real_,
      rr_y_high_model = NA_real_,
      rr_response_contrast = NA_real_,
      standardized_relevant_range_effect = NA_real_,
      rr_x_low_supported = NA,
      rr_x_high_supported = NA,
      rr_x_range_supported = NA,
      stringsAsFactors = FALSE
    )
  }

  cbind(
    data.frame(
      source_node = source_node,
      destination_node = destination_node,
      predictor = predictor,
      response = response,
      random_effects = if (isTRUE(include_random)) "conditional" else "excluded",
      averaging = if (is.null(average_weights)) "equal rows" else "user supplied",
      conditioning = if (length(at) == 0L) {
        "observed covariates"
      } else {
        paste(
          paste0(
            names(at), "=",
            vapply(at, function(x) paste(x, collapse = ","), character(1))
          ),
          collapse = "; "
        )
      },
      stringsAsFactors = FALSE
    ),
    primary,
    oe,
    relevant_range
  )
}

# Configuration-driven wrapper reusable for primary, realized-parasitism,
# network, and alternative-direction SEMs. Required edge columns are edge_id,
# source_node, destination_node, model_id, and data_id. Optional columns are
# predictor, response, average_weight_variable, derivative_step, lower_model,
# upper_model, and a list-column named `at`.
# Details: Resolves the model, native dataset, conditioning, weights, bounds,
# and adapters declared for each edge, with configurable error handling.
# Returns: A data frame with one standardized result or recorded error per edge.
standardize_sem_edges <- function(
    edge_spec,
    models,
    data_sets,
    node_reference,
    include_relevant_range = TRUE,
    include_random = FALSE,
    predict_functions = NULL,
    on_error = c("stop", "record"),
    evaluation_functions = NULL) {
  on_error <- match.arg(on_error)
  required <- c(
    "edge_id", "source_node", "destination_node", "model_id", "data_id"
  )
  missing_columns <- setdiff(required, names(edge_spec))
  if (length(missing_columns) > 0L) {
    stop("`edge_spec` lacks: ", paste(missing_columns, collapse = ", "), ".")
  }
  if (anyDuplicated(edge_spec$edge_id)) stop("`edge_id` values must be unique.")
  if (is.null(names(models)) || is.null(names(data_sets))) {
    stop("`models` and `data_sets` must be named lists.")
  }

  results <- vector("list", nrow(edge_spec))
  for (i in seq_len(nrow(edge_spec))) {
    specification <- edge_spec[i, , drop = FALSE]
    model_id <- as.character(specification$model_id[[1]])
    data_id <- as.character(specification$data_id[[1]])
    if (!model_id %in% names(models)) stop("Unknown model_id: ", model_id)
    if (!data_id %in% names(data_sets)) stop("Unknown data_id: ", data_id)

    get_optional <- function(name, default = NULL) {
      if (!name %in% names(edge_spec)) return(default)
      value <- specification[[name]][[1]]
      if (length(value) == 0L || (length(value) == 1L && is.na(value))) default else value
    }

    at <- if ("at" %in% names(edge_spec)) specification$at[[1]] else list()
    if (is.null(at)) at <- list()
    data_i <- data_sets[[data_id]]
    weight_variable <- get_optional("average_weight_variable")
    average_weights <- if (is.null(weight_variable)) {
      NULL
    } else {
      if (!weight_variable %in% names(data_i)) {
        stop("Unknown averaging-weight variable: ", weight_variable)
      }
      data_i[[weight_variable]]
    }
    bounds <- c(
      get_optional("lower_model", -Inf),
      get_optional("upper_model", Inf)
    )
    predict_function <- if (!is.null(predict_functions) &&
                            model_id %in% names(predict_functions)) {
      predict_functions[[model_id]]
    } else {
      NULL
    }

    calculation <- tryCatch(
      standardize_sem_path(
        model = models[[model_id]],
        data = data_i,
        source_node = as.character(specification$source_node[[1]]),
        destination_node = as.character(specification$destination_node[[1]]),
        node_reference = node_reference,
        predictor = get_optional("predictor"),
        response = get_optional("response"),
        at = at,
        average_weights = average_weights,
        derivative_step = get_optional("derivative_step"),
        predictor_bounds = bounds,
        include_relevant_range = include_relevant_range,
        include_random = include_random,
        predict_function = predict_function,
        evaluation_function = evaluation_functions[[model_id]]
      ),
      error = function(error) error
    )

    if (inherits(calculation, "error")) {
      if (on_error == "stop") {
        stop("Edge `", specification$edge_id[[1]], "`: ", conditionMessage(calculation))
      }
      results[[i]] <- data.frame(
        edge_id = specification$edge_id[[1]],
        source_node = specification$source_node[[1]],
        destination_node = specification$destination_node[[1]],
        model_id = model_id,
        data_id = data_id,
        standardization_error = conditionMessage(calculation),
        stringsAsFactors = FALSE
      )
    } else {
      calculation$edge_id <- specification$edge_id[[1]]
      calculation$model_id <- model_id
      calculation$data_id <- data_id
      calculation$standardization_error <- NA_character_
      results[[i]] <- calculation[, c(
        "edge_id", "model_id", "data_id",
        setdiff(names(calculation), c("edge_id", "model_id", "data_id"))
      ), drop = FALSE]
    }
  }
  all_columns <- unique(unlist(lapply(results, names), use.names = FALSE))
  results <- lapply(results, function(result) {
    missing <- setdiff(all_columns, names(result))
    for (column in missing) result[[column]] <- NA
    result[, all_columns, drop = FALSE]
  })
  do.call(rbind, results)
}
# -----------------------------------------------------------------------------
# Component-model R2 summaries
# -----------------------------------------------------------------------------

# Calculate an unbiased weighted variance for component R2.
# Details: Normalizes non-negative weights and applies the unbiased weighted
# denominator; unweighted input uses the ordinary sample variance.
# Returns: One variance value.
.sem_weighted_variance <- function(x, weights = NULL) {
  if (is.null(weights)) return(stats::var(x))
  if (!is.numeric(weights) || length(weights) != length(x) ||
      any(!is.finite(weights)) || any(weights < 0) || sum(weights) <= 0) {
    stop("`weights` must be finite, non-negative, and match `x`.")
  }
  weights <- weights / sum(weights)
  denominator <- 1 - sum(weights^2)
  if (denominator <= 0) stop("At least two observations must have positive weight.")
  center <- sum(weights * x)
  sum(weights * (x - center)^2) / denominator
}

# Reconstruct the random-effect design column for an intercept or numeric slope.
# Details: Maps a random intercept to ones and a numeric random slope to its
# matching data column, stopping for unsupported term expressions.
# Returns: The numeric design vector for one random-effect term.
.sem_random_design_column <- function(term, data) {
  if (identical(term, "(Intercept)")) return(rep(1, nrow(data)))
  if (term %in% names(data) && is.numeric(data[[term]])) return(data[[term]])
  stop(
    "Cannot reconstruct random-slope term `", term,
    "`. Supply `random_variance` explicitly."
  )
}

# Average the observation-specific variance contributed by random effects.
# Details: Reconstructs each random-effect design matrix, applies its covariance
# matrix, and averages observation-specific variance contributions.
# Returns: One mean random-effect variance across all grouping terms.
.sem_random_effect_variance <- function(model, data, weights = NULL) {
  covariance_blocks <- if (inherits(model, "glmmTMB")) {
    glmmTMB::VarCorr(model)$cond
  } else if (inherits(model, "merMod")) {
    lme4::VarCorr(model)
  } else {
    return(0)
  }
  if (length(covariance_blocks) == 0L) return(0)

  individual_variance <- rep(0, nrow(data))
  for (covariance in covariance_blocks) {
    covariance <- as.matrix(covariance)
    terms <- colnames(covariance)
    if (is.null(terms)) stop("A random-effect covariance block lacks term names.")
    design <- vapply(
      terms,
      .sem_random_design_column,
      numeric(nrow(data)),
      data = data
    )
    if (is.null(dim(design))) design <- matrix(design, ncol = 1L)
    individual_variance <- individual_variance +
      rowSums((design %*% covariance) * design)
  }
  .sem_weighted_mean(individual_variance, weights)
}

# Obtain finite observation-level dispersion predictions on the response scale.
# Details: Obtains observation-level dispersion from glmmTMB, optionally keeping
# dispersion random effects, and requires positive finite values.
# Returns: One response-scale dispersion value per data row.
.sem_dispersion_prediction <- function(
    model,
    data,
    include_dispersion_random = TRUE) {
  if (inherits(model, "glmmTMB")) {
    value <- stats::predict(
      model,
      newdata = data,
      type = "disp",
      re.form = if (isTRUE(include_dispersion_random)) NULL else NA,
      allow.new.levels = TRUE
    )
  } else {
    value <- tryCatch(stats::sigma(model), error = function(error) NA_real_)
    if (!is.finite(value)) value <- 1
  }
  value <- rep(as.numeric(value), length.out = nrow(data))
  if (any(!is.finite(value)) || any(value <= 0)) {
    stop("Dispersion predictions must be positive and finite.")
  }
  value
}

# Recover binomial trial totals from weights or a two-column response.
# Details: Uses explicit totals when supplied; otherwise recovers trials from
# model weights or a two-column success/failure response.
# Returns: One positive finite trial total per observation.
.sem_binomial_trials <- function(model, data, trials = NULL) {
  if (is.character(trials)) {
    .sem_assert_scalar_character(trials, "trials")
    if (!trials %in% names(data)) stop("Unknown trial-count column: ", trials)
    trials <- data[[trials]]
  }
  if (is.null(trials) && inherits(model, "glmmTMB") &&
      !is.null(model$call$weights)) {
    trials <- tryCatch(
      eval(model$call$weights, envir = data, enclos = parent.frame()),
      error = function(error) NULL
    )
  }
  if (is.null(trials) && inherits(model, "glm") &&
      length(model$prior.weights) == nrow(data)) {
    trials <- model$prior.weights
  }
  if (is.null(trials)) trials <- rep(1, nrow(data))
  trials <- rep(as.numeric(trials), length.out = nrow(data))
  if (any(!is.finite(trials)) || any(trials <= 0)) {
    stop("Binomial trial counts must be positive and finite.")
  }
  trials
}

# Calculate family-specific conditional response variance for Nakagawa R2.
# Details: Converts fitted means and dispersion into family-appropriate
# observation-level variance, with an adapter for unsupported families.
# Returns: One positive finite conditional-response variance per row.
.sem_conditional_response_variance <- function(
    model,
    data,
    mu,
    eta,
    dispersion,
    trials = NULL,
    variance_function = NULL) {
  if (!is.null(variance_function)) {
    value <- variance_function(model, data, mu, eta, dispersion, trials)
  } else {
    family_name <- tolower(stats::family(model)$family)
    value <- switch(
      family_name,
      gaussian = dispersion^2,
      poisson = mu,
      nbinom1 = mu * (1 + dispersion),
      nbinom2 = mu * (1 + mu / dispersion),
      genpois = mu * dispersion,
      compois = mu * dispersion,
      binomial = {
        trial_count <- .sem_binomial_trials(model, data, trials)
        mu * (1 - mu) / trial_count
      },
      betabinomial = {
        trial_count <- .sem_binomial_trials(model, data, trials)
        mu * (1 - mu) * (trial_count + dispersion) /
          (trial_count * (dispersion + 1))
      },
      beta = mu * (1 - mu) / (dispersion + 1),
      gamma = mu^2 * dispersion^2,
      inverse.gaussian = mu^3 * dispersion^2,
      tweedie = {
        if (!inherits(model, "glmmTMB")) {
          stop("Supply `variance_function` for a non-glmmTMB Tweedie model.")
        }
        parameters <- glmmTMB::family_params(model)
        power <- unname(parameters[grep("Tweedie power", names(parameters))][1])
        if (!is.finite(power)) stop("Could not recover the Tweedie power.")
        dispersion * mu^power
      },
      t =,
      t_family = {
        if (!inherits(model, "glmmTMB")) {
          stop("Supply `variance_function` for a non-glmmTMB Student-t model.")
        }
        parameters <- glmmTMB::family_params(model)
        degrees_freedom <- unname(
          parameters[grep("Student-t df", names(parameters))][1]
        )
        if (!is.finite(degrees_freedom) || degrees_freedom <= 2) {
          stop("Student-t response variance is undefined when df <= 2.")
        }
        dispersion^2 * degrees_freedom / (degrees_freedom - 2)
      },
      stop(
        "Unsupported family `", stats::family(model)$family,
        "`. Supply `variance_function` explicitly."
      )
    )
  }
  value <- as.numeric(value)
  if (length(value) != nrow(data) || any(!is.finite(value)) || any(value < 0)) {
    stop("Conditional response variances must be finite, non-negative, and row-specific.")
  }
  value
}

# Calculate marginal and conditional Nakagawa R2 for one native SEM component.
# Details: Partitions fixed, random, and conditional-response variance on the
# native sample instead of relying on unsupported package approximations.
# Returns: A one-row data frame with variance components and marginal/conditional R2.
sem_component_nakagawa_r2 <- function(
    model,
    data,
    approximation = c("auto", "lognormal", "delta"),
    trials = NULL,
    average_weights = NULL,
    random_variance = NULL,
    include_dispersion_random = TRUE,
    variance_function = NULL,
    predict_function = NULL,
    mean_derivative_function = NULL) {
  approximation <- match.arg(approximation)
  model_n <- tryCatch(stats::nobs(model), error = function(error) NA_real_)
  if (is.finite(model_n) && model_n != nrow(data)) {
    stop("`data` must be the exact native component-model sample.")
  }
  if (!is.null(average_weights) && length(average_weights) != nrow(data)) {
    stop("`average_weights` must match the rows of `data`.")
  }
  if (inherits(model, "glmmTMB")) {
    zero_inflation <- length(glmmTMB::fixef(model)$zi) > 0L ||
      length(glmmTMB::VarCorr(model)$zi) > 0L
    if (zero_inflation && is.null(variance_function)) {
      stop("Zero-inflated models require an explicit `variance_function`.")
    }
  }

  eta <- .sem_predict(model, data, "link", FALSE, predict_function)
  mu <- .sem_predict(model, data, "response", FALSE, predict_function)
  fixed_variance <- .sem_weighted_variance(eta, average_weights)
  if (!is.finite(fixed_variance) || fixed_variance < 0) {
    stop("The fixed-effect linear-predictor variance is invalid.")
  }

  if (is.null(random_variance)) {
    random_variance <- .sem_random_effect_variance(model, data, average_weights)
  }
  if (!is.numeric(random_variance) || length(random_variance) != 1L ||
      !is.finite(random_variance) || random_variance < 0) {
    stop("`random_variance` must be one finite non-negative number.")
  }

  dispersion <- .sem_dispersion_prediction(
    model, data, include_dispersion_random
  )
  response_variance <- .sem_conditional_response_variance(
    model, data, mu, eta, dispersion, trials, variance_function
  )
  link_name <- stats::family(model)$link
  approximation_used <- if (approximation == "auto") {
    if (identical(link_name, "log")) "lognormal" else "delta"
  } else {
    approximation
  }
  if (approximation_used == "lognormal") {
    if (!identical(link_name, "log") || any(mu <= 0)) {
      stop("The lognormal approximation requires a log link and positive means.")
    }
    observation_variance <- log1p(response_variance / mu^2)
  } else {
    derivative <- if (is.null(mean_derivative_function)) {
      stats::family(model)$mu.eta(eta)
    } else {
      mean_derivative_function(model, data, eta, mu)
    }
    if (any(!is.finite(derivative)) || any(abs(derivative) <= .Machine$double.eps)) {
      stop("The inverse-link derivative is zero or non-finite.")
    }
    observation_variance <- response_variance / derivative^2
  }
  distribution_variance <- .sem_weighted_mean(
    observation_variance, average_weights
  )
  total_variance <- fixed_variance + random_variance + distribution_variance
  if (!is.finite(total_variance) || total_variance <= 0) {
    stop("The total latent-scale variance is not positive and finite.")
  }

  data.frame(
    family = stats::family(model)$family,
    link = link_name,
    n = nrow(data),
    approximation = approximation_used,
    dispersion_random_effects = isTRUE(include_dispersion_random),
    variance_fixed = fixed_variance,
    variance_random = random_variance,
    variance_distribution = distribution_variance,
    r2_marginal = fixed_variance / total_variance,
    r2_conditional = (fixed_variance + random_variance) / total_variance,
    stringsAsFactors = FALSE
  )
}

# Assign whole plots to reproducible, optionally stratified validation folds.
# Details: Assigns complete clusters to folds, optionally within strata, using
# a reproducible seed and approximately balanced cluster counts.
# Returns: One integer fold identifier per input row.
make_grouped_sem_folds <- function(
    data,
    group,
    folds = 5L,
    strata = NULL,
    seed = 1L) {
  .sem_assert_scalar_character(group, "group")
  if (!group %in% names(data)) stop("Unknown grouping column: ", group)
  if (!is.null(strata) && (!is.character(strata) ||
                          !all(strata %in% names(data)))) {
    stop("`strata` must be NULL or names of columns in `data`.")
  }
  if (group %in% strata) stop("`group` cannot also be a stratum column.")
  if (!is.numeric(folds) || length(folds) != 1L || !is.finite(folds) ||
      folds < 2L || folds %% 1 != 0) {
    stop("`folds` must be an integer of at least two.")
  }
  folds <- as.integer(folds)
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
    stop("`seed` must be one finite number.")
  }
  grouping_columns <- c(strata, group)
  groups <- unique(data[, grouping_columns, drop = FALSE])
  if (nrow(groups) < folds) stop("There are fewer groups than folds.")
  if (anyNA(groups)) stop("Grouping and stratum identifiers cannot be missing.")

  stratum_key <- if (is.null(strata) || length(strata) == 0L) {
    factor(rep("all", nrow(groups)))
  } else {
    do.call(
      interaction,
      c(unname(groups[, strata, drop = FALSE]), list(drop = TRUE, lex.order = TRUE))
    )
  }
  set.seed(seed)
  groups$.fold <- NA_integer_
  for (level in levels(stratum_key)) {
    indices <- which(stratum_key == level)
    indices <- sample(indices, length(indices))
    groups$.fold[indices] <- rep(seq_len(folds), length.out = length(indices))
  }

  row_key <- do.call(
    interaction,
    c(unname(data[, grouping_columns, drop = FALSE]), list(drop = TRUE))
  )
  group_key <- do.call(
    interaction,
    c(unname(groups[, grouping_columns, drop = FALSE]), list(drop = TRUE))
  )
  assigned <- groups$.fold[match(as.character(row_key), as.character(group_key))]
  if (anyNA(assigned) || !setequal(unique(assigned), seq_len(folds))) {
    stop("Fold assignment failed or left at least one fold empty.")
  }
  assigned
}

# Calculate grouped, out-of-sample response-scale predictive R2.
# Details: Refits on grouped training folds, predicts held-out rows, and compares
# model SSE with a training-mean null benchmark on the response scale.
# Returns: A sem_cv_r2 list with summary, fold diagnostics, and predictions.
sem_component_predictive_r2 <- function(
    model,
    data,
    group,
    folds = 5L,
    strata = NULL,
    seed = 1L,
    response = NULL,
    score_weights = NULL,
    fold_id = NULL,
    update_arguments = list(),
    refit_function = NULL,
    predict_function = NULL,
    validate_refits = TRUE,
    on_refit_failure = c("stop", "record")) {
  on_refit_failure <- match.arg(on_refit_failure)
  if (!is.list(update_arguments)) stop("`update_arguments` must be a list.")
  model_n <- tryCatch(stats::nobs(model), error = function(error) NA_real_)
  if (is.finite(model_n) && model_n != nrow(data)) {
    stop("`data` must be the exact native component-model sample.")
  }
  if (is.null(response)) {
    response <- all.vars(stats::formula(model)[[2L]])
    if (length(response) != 1L) {
      stop("Supply `response` explicitly for a non-scalar model response.")
    }
  }
  .sem_assert_scalar_character(response, "response")
  if (!response %in% names(data) || !is.numeric(data[[response]])) {
    stop("`response` must identify one numeric column in `data`.")
  }
  if (any(!is.finite(data[[response]]))) {
    stop("The response must be finite on every native-sample row.")
  }
  if (is.character(score_weights)) {
    .sem_assert_scalar_character(score_weights, "score_weights")
    if (!score_weights %in% names(data)) stop("Unknown scoring-weight column.")
    score_weights <- data[[score_weights]]
  }
  if (is.null(score_weights)) score_weights <- rep(1, nrow(data))
  if (!is.numeric(score_weights) || length(score_weights) != nrow(data) ||
      any(!is.finite(score_weights)) || any(score_weights < 0) ||
      sum(score_weights) <= 0) {
    stop("Scoring weights must be finite, non-negative, and match `data`.")
  }
  if (is.null(fold_id)) {
    fold_id <- make_grouped_sem_folds(data, group, folds, strata, seed)
  }
  if (length(fold_id) != nrow(data) || anyNA(fold_id)) {
    stop("`fold_id` must contain one non-missing fold identifier per row.")
  }
  cluster_key <- do.call(
    interaction,
    c(
      unname(data[, c(strata, group), drop = FALSE]),
      list(drop = TRUE, lex.order = TRUE)
    )
  )
  folds_per_cluster <- tapply(fold_id, cluster_key, function(x) length(unique(x)))
  if (any(folds_per_cluster != 1L)) {
    stop("Every complete group must be assigned to exactly one fold.")
  }
  fold_levels <- sort(unique(fold_id))
  if (length(fold_levels) < 2L) stop("At least two folds are required.")

  if (is.null(refit_function)) {
    refit_function <- function(model, training_data, update_arguments) {
      do.call(
        stats::update,
        c(list(object = model, data = training_data), update_arguments)
      )
    }
  }
  predictions <- rep(NA_real_, nrow(data))
  null_predictions <- rep(NA_real_, nrow(data))
  fold_rows <- vector("list", length(fold_levels))

  for (i in seq_along(fold_levels)) {
    test <- fold_id == fold_levels[[i]]
    training_data <- data[!test, , drop = FALSE]
    test_data <- data[test, , drop = FALSE]
    refit <- tryCatch(
      list(model = refit_function(model, training_data, update_arguments), error = NULL),
      error = function(error) list(model = NULL, error = conditionMessage(error))
    )
    failure_reason <- refit$error
    fitted_model <- refit$model
    if (is.null(failure_reason) && isTRUE(validate_refits)) {
      validation <- validate_sem_refit(fitted_model)
      if (!validation$ok) {
        failure_reason <- validation$reason
      }
    }
    if (!is.null(failure_reason)) {
      if (identical(on_refit_failure, "stop")) {
        stop("Fold ", fold_levels[[i]], " refit failed: ", failure_reason, ".")
      }
      fold_rows[[i]] <- data.frame(
        fold = fold_levels[[i]],
        n_test = sum(test),
        refit_status = "failed",
        refit_reason = failure_reason,
        sse_model = NA_real_,
        sse_null = NA_real_,
        r2_response = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }
    predictions[test] <- .sem_predict(
      fitted_model, test_data, "response", FALSE, predict_function
    )
    training_weights <- score_weights[!test]
    null_predictions[test] <- .sem_weighted_mean(
      training_data[[response]], training_weights
    )

    fold_weights <- score_weights[test]
    fold_sse <- sum(fold_weights * (test_data[[response]] - predictions[test])^2)
    fold_null_sse <- sum(
      fold_weights * (test_data[[response]] - null_predictions[test])^2
    )
    fold_rows[[i]] <- data.frame(
      fold = fold_levels[[i]],
      n_test = sum(test),
      refit_status = "available",
      refit_reason = NA_character_,
      sse_model = fold_sse,
      sse_null = fold_null_sse,
      r2_response = 1 - fold_sse / fold_null_sse,
      stringsAsFactors = FALSE
    )
  }

  fold_results <- do.call(rbind, fold_rows)
  failed_folds <- fold_results$fold[fold_results$refit_status == "failed"]
  cv_available <- length(failed_folds) == 0L
  cv_reason <- if (cv_available) {
    NA_character_
  } else {
    paste(
      paste0(
        "fold ", failed_folds, ": ",
        fold_results$refit_reason[fold_results$refit_status == "failed"]
      ),
      collapse = " | "
    )
  }

  observed <- data[[response]]
  if (cv_available) {
    sse_model <- sum(score_weights * (observed - predictions)^2)
    sse_null <- sum(score_weights * (observed - null_predictions)^2)
    if (!is.finite(sse_null) || sse_null <= 0) {
      stop("The cross-validated null-model squared error is not positive.")
    }
    positive_weight <- score_weights > 0
    correlation_r2 <- if (sum(positive_weight) >= 3L &&
                          stats::sd(observed[positive_weight]) > 0 &&
                          stats::sd(predictions[positive_weight]) > 0) {
      stats::cor(observed[positive_weight], predictions[positive_weight])^2
    } else {
      NA_real_
    }
  } else {
    sse_model <- NA_real_
    sse_null <- NA_real_
    correlation_r2 <- NA_real_
  }

  structure(
    list(
      summary = data.frame(
        response = response,
        family = stats::family(model)$family,
        link = stats::family(model)$link,
        n = nrow(data),
        folds = length(fold_levels),
        grouping = group,
        stratification = if (is.null(strata)) NA_character_ else
          paste(strata, collapse = ":"),
        weighting = if (all(score_weights == 1)) "equal rows" else "user supplied",
        cv_status = if (cv_available) "available" else "unavailable",
        cv_reason = cv_reason,
        n_folds_failed = length(failed_folds),
        sse_model = sse_model,
        sse_null = sse_null,
        r2_response_cv = if (cv_available) 1 - sse_model / sse_null else NA_real_,
        r2_correlation_cv = correlation_r2,
        stringsAsFactors = FALSE
      ),
      folds = fold_results,
      predictions = data.frame(
        row = seq_len(nrow(data)),
        fold = fold_id,
        observed = observed,
        predicted = predictions,
        null_predicted = null_predictions,
        weight = score_weights,
        stringsAsFactors = FALSE
      )
    ),
    class = "sem_cv_r2"
  )
}
# Build the complete inference, edge, and standardization objects for one SEM
# component. This is reusable when a workflow substitutes or adds a response
# component while retaining the shared node and data registries.
# Details: Converts one fitted component into shared metadata, term inference,
# directed-edge specifications, standardized effects, and reporting tables.
# Returns: A named list containing the model and all component-level SEM objects.
build_sem_component_results <- function(
    model,
    component,
    model_id,
    data_id,
    data_sets,
    node_reference,
    reference_values,
    include_relevant_range = TRUE,
    nuisance_predictors = c(
      "site", "sc_slope", "sc_elev", "sc_east", "sc_north"
    ),
    forest_age_node = "forest_age",
    coefficient_table = NULL,
    inference_function = sem_term_wald_test,
    predict_function = NULL,
    evaluation_function = NULL) {
  .sem_assert_scalar_character(component, "component")
  .sem_assert_scalar_character(model_id, "model_id")
  .sem_assert_scalar_character(data_id, "data_id")
  .sem_assert_scalar_character(forest_age_node, "forest_age_node")
  focal_dbi <- inherits(model, "gamlss") && identical(model$family[[1L]], "DBI")
  if (focal_dbi) {
    if (!isTRUE(model$converged) || is.null(coefficient_table) ||
        is.null(predict_function) || identical(inference_function, sem_term_wald_test)) {
      stop("A DBI component requires explicit prediction and bootstrap inference adapters.")
    }
  } else if (!inherits(model, "glmmTMB") ||
             !isTRUE(model$fit$convergence == 0L) ||
             !isTRUE(model$sdr$pdHess) ||
             !all(is.finite(glmmTMB::fixef(model)$cond))) {
    stop("`model` must be a valid converged glmmTMB or adapted DBI component fit.")
  }
  if (is.null(names(data_sets)) || !data_id %in% names(data_sets)) {
    stop("`data_id` must identify a named native component dataset.")
  }

  response <- all.vars(stats::formula(model))[[1]]
  component_metadata <- tibble::tibble(
    component = component,
    Response = response,
    source_model = model_id,
    data_id = data_id,
    n = if (focal_dbi) model$noObs else stats::nobs(model),
    family = if (focal_dbi) "DBI" else stats::family(model)$family,
    link = if (focal_dbi) model$mu.link else stats::family(model)$link
  )

  coefficients <- if (is.null(coefficient_table)) {
    summary(model)$coefficients$cond
  } else coefficient_table
  edge_table <- tibble::tibble(
    Response = response,
    Predictor = rownames(coefficients),
    Estimate = unname(coefficients[, "Estimate"]),
    Std.Error = unname(coefficients[, "Std. Error"]),
    Crit.Value = unname(coefficients[, "z value"]),
    P.Value = unname(coefficients[, "Pr(>|z|)"])
  ) |>
    dplyr::left_join(component_metadata, by = "Response") |>
    dplyr::filter(
      Predictor != "(Intercept)",
      !Predictor %in% nuisance_predictors,
      !startsWith(Predictor, "site")
    ) |>
    dplyr::mutate(
      term_type = dplyr::if_else(
        grepl(":", Predictor, fixed = TRUE),
        "interaction",
        "direct_path"
      ),
      interaction_type = dplyr::case_when(
        term_type == "direct_path" ~ NA_character_,
        grepl("sc_fa", Predictor, fixed = TRUE) ~ "age_moderation",
        TRUE ~ "ecological_moderation"
      ),
      conf_low_link = if ("Conf.Low" %in% colnames(coefficients)) {
        coefficients[Predictor, "Conf.Low"]
      } else Estimate - stats::qnorm(0.975) * Std.Error,
      conf_high_link = if ("Conf.High" %in% colnames(coefficients)) {
        coefficients[Predictor, "Conf.High"]
      } else Estimate + stats::qnorm(0.975) * Std.Error
    )
  edge_table <- edge_table |>
    dplyr::bind_cols(dplyr::bind_rows(Map(
      function(term, term_type) {
        inference_function(model, term, term_type, reference_values)
      },
      edge_table$Predictor,
      edge_table$term_type
    )))

  source_lookup <- node_reference |>
    dplyr::select(source_node = node, Predictor = model_variable)
  destination_lookup <- node_reference |>
    dplyr::select(destination_node = node, Response = model_variable)
  forest_age_variable <- node_reference |>
    dplyr::filter(node == forest_age_node) |>
    dplyr::pull(model_variable)
  if (length(forest_age_variable) != 1L) {
    stop("`forest_age_node` must identify exactly one node-reference row.")
  }

  age_moderated <- edge_table |>
    dplyr::filter(interaction_type == "age_moderation") |>
    dplyr::transmute(
      Response,
      Predictor = sub("(^sc_fa:|:sc_fa$)", "", Predictor),
      age_moderated = TRUE
    )

  edge_spec <- edge_table |>
    dplyr::filter(
      term_type == "direct_path",
      Predictor != forest_age_variable
    ) |>
    dplyr::left_join(source_lookup, by = "Predictor") |>
    dplyr::left_join(destination_lookup, by = "Response") |>
    dplyr::left_join(age_moderated, by = c("Response", "Predictor")) |>
    dplyr::mutate(
      age_moderated = dplyr::coalesce(age_moderated, FALSE),
      edge_id = paste(source_node, destination_node, sep = "__"),
      model_id = source_model,
      predictor = Predictor,
      response = Response,
      at = lapply(
        age_moderated,
        function(is_moderated) {
          if (is_moderated) list(sc_fa = 0) else list()
        }
      ),
      average_weight_variable = NA_character_
    ) |>
    dplyr::select(
      edge_id, source_node, destination_node,
      model_id, data_id, predictor, response,
      age_moderated, at, average_weight_variable
    )

  age_effect_spec <- edge_table |>
    dplyr::filter(
      term_type == "direct_path",
      Predictor == forest_age_variable
    ) |>
    dplyr::left_join(destination_lookup, by = "Response") |>
    dplyr::mutate(
      edge_id = paste(forest_age_node, destination_node, sep = "__"),
      source_node = forest_age_node,
      model_id = source_model,
      predictor = Predictor,
      response = Response,
      age_moderated = NA,
      at = rep(list(list()), dplyr::n()),
      average_weight_variable = NA_character_
    ) |>
    dplyr::select(
      edge_id, source_node, destination_node,
      model_id, data_id, predictor, response,
      age_moderated, at, average_weight_variable
    )

  standardization_spec <- dplyr::bind_rows(edge_spec, age_effect_spec)
  standardized_effects <- standardize_sem_edges(
    edge_spec = standardization_spec,
    models = stats::setNames(list(model), model_id),
    data_sets = data_sets,
    node_reference = node_reference,
    include_relevant_range = include_relevant_range,
    include_random = FALSE,
    on_error = "stop",
    predict_functions = stats::setNames(list(predict_function), model_id),
    evaluation_functions = stats::setNames(list(evaluation_function), model_id)
  )
  term_standardizations <- standardized_effects |>
    dplyr::left_join(
      standardization_spec |>
        dplyr::select(edge_id, age_moderated, average_weight_variable),
      by = "edge_id"
    )

  term_results <- edge_table |>
    dplyr::mutate(
      term_role = dplyr::case_when(
        term_type == "interaction" & interaction_type == "age_moderation" ~
          "age_moderation",
        term_type == "interaction" ~ "ecological_interaction",
        Predictor == forest_age_variable ~ "forest_age_main_effect",
        TRUE ~ "main_ecological_path"
      )
    ) |>
    dplyr::transmute(
      component,
      model_id = source_model,
      data_id,
      response = Response,
      predictor = Predictor,
      term_role,
      interaction_type,
      native_n = n,
      family,
      link,
      coefficient_estimate = Estimate,
      coefficient_std_error = Std.Error,
      coefficient_conf_low = conf_low_link,
      coefficient_conf_high = conf_high_link,
      coefficient_test_statistic = Crit.Value,
      coefficient_p_value = P.Value,
      inference_estimand,
      inference_method,
      reference_conditioning,
      estimate = inference_estimate,
      std_error = inference_standard_error,
      conf_low = inference_conf_low,
      conf_high = inference_conf_high,
      df = NA_real_,
      test_statistic = inference_test_statistic,
      p_value = inference_p_value
    ) |>
    dplyr::left_join(
      term_standardizations |>
        dplyr::select(
          model_id, data_id, response, predictor,
          edge_id, source_node, destination_node, age_moderated,
          average_weight_variable, random_effects, averaging, conditioning,
          n_evaluation, derivative_step, average_partial_effect,
          sd_x_reference, sd_y_reference, standardized_response_effect,
          oe_link_slope, oe_r2, oe_variance_eta, oe_sd_y,
          standardized_oe_effect,
          rr_x_low_model, rr_x_high_model, rr_y_low_model, rr_y_high_model,
          rr_response_contrast, standardized_relevant_range_effect,
          rr_x_low_supported, rr_x_high_supported, rr_x_range_supported,
          standardization_error,
          dplyr::any_of(c("rr_n_evaluation", "rr_n_excluded"))
        ),
      by = c("model_id", "data_id", "response", "predictor")
    ) |>
    dplyr::mutate(
      standardization_status = dplyr::case_when(
        term_role %in% c("main_ecological_path", "forest_age_main_effect") &
          is.na(standardization_error) ~ "available",
        term_role %in% c("main_ecological_path", "forest_age_main_effect") ~
          paste0("failed: ", standardization_error),
        TRUE ~ "not calculated: interaction"
      ),
      .after = age_moderated
    )

  standardized_rows <- term_results$term_role %in%
    c("main_ecological_path", "forest_age_main_effect")
  stopifnot(
    !anyNA(standardization_spec$source_node),
    !anyNA(standardization_spec$destination_node),
    !anyDuplicated(standardization_spec$edge_id),
    all(is.na(standardized_effects$standardization_error)),
    all(is.finite(standardized_effects$standardized_response_effect)),
    all(is.finite(standardized_effects$standardized_oe_effect)),
    all(is.finite(term_results$estimate)),
    all(is.finite(term_results$std_error)),
    all(is.finite(term_results$standardized_response_effect[standardized_rows])),
    all(is.na(term_results$standardized_response_effect[!standardized_rows]))
  )
  if (include_relevant_range) {
    stopifnot(
      all(is.finite(standardized_effects$standardized_relevant_range_effect)),
      all(standardized_effects$rr_x_range_supported)
    )
  }

  list(
    model = stats::setNames(list(model), model_id),
    component_metadata = component_metadata,
    edge_table = edge_table,
    edge_spec = edge_spec,
    standardization_spec = standardization_spec,
    standardized_effects = standardized_effects,
    term_results = term_results
  )
}
# -----------------------------------------------------------------------------
# Network-SEM focal inference, feasibility, R2, and export
# -----------------------------------------------------------------------------
# These helpers support the bounded focal-richness component and the
# part-whole richness constraints specific to Network_SEM.R.

# Convert a requested bootstrap success threshold into a replicate count.
# Details: Values at or below one are treated as proportions; larger values
# are treated as counts. Impossible or non-positive thresholds stop early.
# Returns: The minimum required number of successful bootstrap refits.
.sem_required_successes <- function(total, minimum_success) {
  if (!is.numeric(minimum_success) || length(minimum_success) != 1L ||
      !is.finite(minimum_success) || minimum_success <= 0) {
    stop("`minimum_success` must be a positive count or proportion.")
  }
  required <- if (minimum_success <= 1) {
    ceiling(total * minimum_success)
  } else {
    ceiling(minimum_success)
  }
  if (required > total) stop("Required successful replicates exceed `times`.")
  required
}

# Run the hierarchical parametric bootstrap for focal parasitoid richness.
# Details: Simulates bounded DBI counts with new plot and year effects, refits
# the unchanged location and dispersion model, retries failed fits, and keeps
# the response, seeds, warnings, errors, and fit diagnostics for every attempt.
# Percentile intervals apply to location-logit coefficients, not the full SEM.
# Returns: Coefficient summaries, replicate audits, draws, settings, and fits.
bootstrap_focal_dbi <- function(
    model, data, times = 2000L, seed = 123L, confidence = 0.95,
    minimum_success = 0.95, progress_every = 10L,
    cycle_limits = c(300L, 1000L, 1000L), tolerance = 1e-6) {
  if (!requireNamespace("gamlss", quietly = TRUE) ||
      !requireNamespace("gamlss.dist", quietly = TRUE)) {
    stop("The focal bootstrap requires gamlss and gamlss.dist.")
  }
  valid_integer <- function(x, lower = 0) {
    is.numeric(x) && length(x) == 1L && is.finite(x) &&
      x >= lower && x <= .Machine$integer.max && x == floor(x)
  }
  if (!valid_integer(times, 2) || !valid_integer(seed) ||
      !valid_integer(progress_every) || seed + 2 * times - 1 > .Machine$integer.max) {
    stop("Invalid replicate count, seed or progress interval.")
  }
  if (length(cycle_limits) != 3L ||
      !all(vapply(cycle_limits, valid_integer, logical(1), lower = 1)) ||
      length(confidence) != 1L || !is.finite(confidence) ||
      confidence <= 0 || confidence >= 1 ||
      length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0) {
    stop("Supply three positive cycle limits and valid confidence/tolerance values.")
  }
  required_successes <- max(2L, .sem_required_successes(times, minimum_success))
  times <- as.integer(times)
  if (!inherits(model, "gamlss") || !isTRUE(model$converged) ||
      !identical(model$family[[1L]], "DBI") ||
      !identical(model$mu.link, "logit") || !identical(model$sigma.link, "log")) {
    stop("Supply the converged focal DBI model with logit/log links.")
  }
  mu_formula <- stats::formula(model, what = "mu")
  sigma_formula <- stats::formula(model, what = "sigma")
  random_terms <- c("random(plot_id)", "random(year)")
  if (!identical(colnames(model$mu.s), random_terms) ||
      length(attr(stats::terms(sigma_formula), "term.labels")) != 0L ||
      attr(stats::terms(sigma_formula), "intercept") != 1L ||
      !identical(mu_formula[[2L]],
                 quote(cbind(focal_rich, focal_max - focal_rich))) ||
      any(model$weights != 1) || any(model$mu.offset != 0) ||
      any(model$sigma.offset != 0)) {
    stop("Expected unweighted focal counts, constant dispersion and plot/year random intercepts.")
  }
  needed <- unique(c(all.vars(mu_formula), "par_rich"))
  if (!is.data.frame(data) || !all(needed %in% names(data))) {
    stop("The focal model's analysis columns are required.")
  }
  data <- as.data.frame(data)[, needed, drop = FALSE]
  if (anyNA(data) || nrow(data) != model$noObs ||
      !identical(as.numeric(data$focal_rich), as.numeric(model$y)) ||
      !identical(as.numeric(data$focal_max), as.numeric(model$bd)) ||
      any(data$focal_max != pmin(4L, data$par_rich)) ||
      any(data$focal_max < 1 | data$focal_max %% 1 != 0) ||
      any(data$focal_rich < 0 | data$focal_rich > data$focal_max |
            data$focal_rich %% 1 != 0)) {
    stop("Supply the same complete, ordered and bounded observations as the fitted model.")
  }
  beta <- stats::coef(model, what = "mu")
  beta <- beta[!is.na(beta)]
  fixed_formula <- stats::update.formula(
    mu_formula, ~ . - random(plot_id) - random(year)
  )
  fixed_contrasts <- attr(model$mu.x, "contrasts")
  fixed_contrasts <- fixed_contrasts[
    intersect(names(fixed_contrasts), all.vars(fixed_formula))
  ]
  x <- stats::model.matrix(fixed_formula, data = data,
                           contrasts.arg = fixed_contrasts)
  if (!all(names(beta) %in% colnames(x)) ||
      !isTRUE(all.equal(unname(x[, names(beta), drop = FALSE]),
                        unname(model$mu.x[, names(beta), drop = FALSE])))) {
    stop("Predictor values or row order differ from the fitted model.")
  }
  smoothers <- lapply(1:2, function(j) gamlss::getSmo(model, "mu", which = j))
  group_index <- lapply(c("plot_id", "year"), function(nm) {
    as.integer(factor(data[[nm]]))
  })
  # Check grouping against the factor labels saved by the random() smoothers.
  for (j in 1:2) {
    saved <- match(smoothers[[j]]$factor, unique(smoothers[[j]]$factor))
    supplied <- match(group_index[[j]], unique(group_index[[j]]))
    if (!identical(saved, supplied)) stop("Random-effect grouping differs from the fit.")
  }
  random_sd <- vapply(smoothers, function(s) s$sigb, numeric(1))
  sigma_hat <- stats::fitted(model, what = "sigma")
  if (any(!is.finite(c(beta, random_sd, sigma_hat))) ||
      any(random_sd < 0) || any(sigma_hat <= 0) ||
      diff(range(sigma_hat)) > 1e-8) stop("Invalid generating parameters.")
  sigma_hat <- sigma_hat[[1L]]
  # Match the pilot: remove the stored smooth contributions from the fitted
  # predictor. GAMLSS backfitting can leave a small difference from X %*% beta;
  # retain that discrepancy for review instead of silently changing the pilot.
  eta_fixed <- as.numeric(stats::qlogis(stats::fitted(model, what = "mu")) -
                            rowSums(model$mu.s))
  if (any(!is.finite(eta_fixed))) stop("Invalid fitted location predictor.")
  predictor_difference <- max(abs(
    eta_fixed - as.numeric(x[, names(beta), drop = FALSE] %*% beta)
  ))

  # Isolate simulation seeds from refit RNG use, and restore the caller's RNG.
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = ".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  simulation_seeds <- as.integer(seed + 2 * (seq_len(times) - 1L))
  refit_seeds <- simulation_seeds + 1L
  n <- nrow(data)
  simulated_response <- matrix(NA_integer_, n, times,
    dimnames = list(rownames(data), as.character(seq_len(times))))
  draws <- matrix(NA_real_, times, length(beta),
    dimnames = list(as.character(seq_len(times)), names(beta)))
  audit <- data.frame(
    replicate = seq_len(times), simulation_seed = simulation_seeds,
    refit_seed = refit_seeds, accepted = FALSE, attempts = 0L,
    selected_attempt = NA_integer_, cycles = NA_integer_,
    global_deviance = NA_real_, sigma = NA_real_, plot_SD = NA_real_,
    year_SD = NA_real_, elapsed_seconds = NA_real_,
    warning_text = "", error_text = "", stringsAsFactors = FALSE
  )
  attempt_records <- list()
  retry_starts <- c(sigma_hat, sigma_hat, 1.4)
  started <- Sys.time()
  # Preserve formula lookup without requiring gamlss to be attached by the caller.
  environment(mu_formula) <- list2env(
    list(random = gamlss::random), parent = environment(mu_formula)
  )

  for (b in seq_len(times)) {
    tick <- proc.time()[["elapsed"]]
    set.seed(simulation_seeds[[b]])
    simulation_warnings <- character()
    y <- tryCatch(withCallingHandlers({
      eta <- eta_fixed
      for (j in 1:2) {
        u <- stats::rnorm(max(group_index[[j]]), sd = random_sd[[j]])
        eta <- eta + u[group_index[[j]]]
      }
      mu <- stats::plogis(eta)
      vapply(seq_len(n), function(i) {
        support <- seq.int(0L, data$focal_max[[i]])
        p <- gamlss.dist::dDBI(support, mu = mu[[i]], sigma = sigma_hat,
                              bd = data$focal_max[[i]])
        if (any(!is.finite(p)) || any(p < 0) || abs(sum(p) - 1) > 1e-6) {
          stop("Invalid finite-support simulation probabilities.")
        }
        as.integer(sample(support, 1L, prob = p / sum(p)))
      }, integer(1))
    }, warning = function(w) {
      simulation_warnings <<- c(simulation_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) e)
    audit$warning_text[[b]] <- paste(unique(simulation_warnings), collapse = " | ")
    if (inherits(y, "error")) {
      audit$error_text[[b]] <- paste("Simulation:", conditionMessage(y))
    } else {
      simulated_response[, b] <- y
      boot_data <- data
      boot_data$focal_rich <- y
      replicate_warnings <- simulation_warnings
      for (attempt in seq_along(cycle_limits)) {
        messages <- character()
        set.seed(refit_seeds[[b]])
        a <- data.frame(
          replicate = b, attempt = attempt, refit_seed = refit_seeds[[b]],
          cycle_limit = cycle_limits[[attempt]], sigma_start = retry_starts[[attempt]],
          converged = FALSE, accepted = FALSE, cycles = NA_integer_,
          global_deviance = NA_real_, sigma = NA_real_, plot_SD = NA_real_,
          year_SD = NA_real_, warning_text = "", error_text = ""
        )
        result <- tryCatch(withCallingHandlers({
          fit <- gamlss::gamlss(
            formula = mu_formula, sigma.formula = ~ 1,
            family = gamlss.dist::DBI(mu.link = "logit", sigma.link = "log"),
            sigma.start = retry_starts[[attempt]], data = boot_data,
            control = gamlss::gamlss.control(
              n.cyc = cycle_limits[[attempt]], c.crit = tolerance, trace = FALSE
            )
          )
          a$converged <- isTRUE(fit$converged)
          a$cycles <- fit$iter
          a$global_deviance <- fit$G.deviance
          a$sigma <- stats::fitted(fit, what = "sigma")[[1L]]
          a$plot_SD <- gamlss::getSmo(fit, "mu", which = 1)$sigb
          a$year_SD <- gamlss::getSmo(fit, "mu", which = 2)$sigb
          bb <- stats::coef(fit, what = "mu")[names(beta)]
          if (!a$converged) stop("Refit did not converge.")
          if (any(!is.finite(c(bb, a$global_deviance, a$sigma, a$plot_SD, a$year_SD))) ||
              a$sigma <= 0 || a$plot_SD < 0 || a$year_SD < 0) {
            stop("Invalid refitted parameter estimates.")
          }
          bb
        }, warning = function(w) {
          messages <<- c(messages, conditionMessage(w))
          invokeRestart("muffleWarning")
        }), error = function(e) e)
        a$warning_text <- paste(unique(messages), collapse = " | ")
        a$accepted <- !inherits(result, "error")
        if (!a$accepted) a$error_text <- conditionMessage(result)
        attempt_records[[length(attempt_records) + 1L]] <- a
        audit$attempts[[b]] <- attempt
        replicate_warnings <- unique(c(replicate_warnings, messages))
        audit$warning_text[[b]] <- paste(replicate_warnings, collapse = " | ")
        audit$error_text[[b]] <- a$error_text
        if (a$accepted) {
          draws[b, ] <- result
          audit$accepted[[b]] <- TRUE
          audit$selected_attempt[[b]] <- attempt
          for (nm in c("cycles", "global_deviance", "sigma", "plot_SD", "year_SD")) {
            audit[b, nm] <- a[[nm]]
          }
          break
        }
      }
    }
    audit$elapsed_seconds[[b]] <- proc.time()[["elapsed"]] - tick
    if (progress_every > 0 && (b %% progress_every == 0 || b == times)) {
      message("Focal DBI bootstrap: ", b, "/", times, "; accepted ",
              sum(audit$accepted), "; retried ", sum(audit$attempts > 1L), ".")
    }
  }
  n_success <- sum(audit$accepted)
  sufficient <- n_success >= required_successes
  coefficients <- data.frame(
    term = names(beta), estimate = unname(beta), bootstrap_mean = NA_real_,
    estimated_bias = NA_real_, bootstrap_SE = NA_real_,
    percentile_lower = NA_real_, percentile_upper = NA_real_,
    confidence = confidence, n_success = n_success, n_requested = times,
    interval_status = if (sufficient) "available_model_based" else "insufficient_refits",
    row.names = NULL
  )
  if (n_success >= 2L) {
    successful <- draws[audit$accepted, , drop = FALSE]
    coefficients$bootstrap_mean <- colMeans(successful)
    coefficients$estimated_bias <- colMeans(successful) - unname(beta)
    coefficients$bootstrap_SE <- apply(successful, 2, stats::sd)
    if (sufficient) {
      limits <- apply(successful, 2, stats::quantile,
                      probs = c((1 - confidence) / 2, (1 + confidence) / 2), names = FALSE)
      coefficients$percentile_lower <- limits[1L, ]
      coefficients$percentile_upper <- limits[2L, ]
    }
  }
  if (!sufficient) warning("Too few successful focal refits; percentile intervals withheld.")
  attempts <- if (length(attempt_records)) do.call(rbind, attempt_records) else data.frame(
    replicate = integer(), attempt = integer(), warning_text = character(), error_text = character()
  )
  list(
    coefficients = coefficients, audit = audit, attempts = attempts, draws = draws,
    simulated_response = simulated_response, reference_model = model, data = data,
    generating_parameters = list(beta = beta, fixed_predictor = eta_fixed,
                                max_fixed_predictor_difference = predictor_difference,
                                sigma = sigma_hat,
                                plot_SD = random_sd[[1L]], year_SD = random_sd[[2L]]),
    settings = list(
      times = times, seed = seed, confidence = confidence,
      required_successes = required_successes, cycle_limits = cycle_limits,
      sigma_starts = retry_starts, tolerance = tolerance,
      rng_kind = RNGkind(), started = started, finished = Sys.time(),
      method = "Hierarchical parametric bootstrap; percentile intervals on location-logit scale",
      limitations = paste(
        "Conditional on observed predictors and the DBI/Gaussian random-effect model;",
        "GAMLSS estimation approximations remain. No complete-SEM uncertainty.",
        "Accepted fits with warnings require review. Intervals exclude failed refits."
      )
    ),
    summary = data.frame(
      requested = times, accepted = n_success, failed = times - n_success,
      retried = sum(audit$attempts > 1L),
      accepted_after_retry = sum(audit$accepted & audit$attempts > 1L),
      replicates_with_warnings = sum(nzchar(audit$warning_text)),
      accepted_attempts_with_warnings = sum(attempts$accepted & nzchar(attempts$warning_text)),
      success_rate = n_success / times, intervals_available = sufficient,
      elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
    ),
    session_info = utils::sessionInfo()
  )
}

# Save a complete focal-richness bootstrap and its readable audit tables.
# Details: Adds checksums for the supplied provenance files, preserves the full
# result as RDS, and exports coefficients, summaries, attempts, and draws as CSV.
# Returns: The output directory invisibly after all files have been written.
export_focal_dbi_bootstrap <- function(result, output_dir, source_paths = character()) {
  .sem_assert_scalar_character(output_dir, "output_dir")
  if (length(source_paths) && any(!file.exists(source_paths))) {
    stop("A bootstrap provenance file is missing.")
  }
  result$source_checksums <- data.frame(
    path = normalizePath(source_paths, mustWork = TRUE),
    md5 = unname(tools::md5sum(source_paths)), row.names = NULL
  )
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (!dir.exists(output_dir)) stop("Could not create bootstrap output directory.")
  tables <- list(
    coefficients = result$coefficients, summary = result$summary,
    replicate_audit = result$audit, attempt_audit = result$attempts,
    coefficient_draws = data.frame(
      replicate = result$audit$replicate, accepted = result$audit$accepted,
      result$draws, check.names = FALSE
    )
  )
  saveRDS(result, file.path(output_dir, "focal_dbi_bootstrap.rds"))
  for (nm in names(tables)) {
    utils::write.csv(tables[[nm]],
      file.path(output_dir, paste0("focal_dbi_", nm, ".csv")),
      row.names = FALSE, na = "NA"
    )
  }
  invisible(output_dir)
}

# Calculate exact conditional moments for the bounded focal-richness model.
# Details: Reconstructs fixed-effect location predictors and evaluates the full
# DBI probability mass function at each observation-specific richness bound.
# Returns: Link predictors, expected counts, variances, and upper bounds.
sem_focal_dbi_moments <- function(model, newdata, eta_shift = 0) {
  if (!inherits(model, "gamlss") || !identical(model$family[[1L]], "DBI") ||
      model$mu.link != "logit" || model$sigma.link != "log") {
    stop("Expected the focal DBI model with logit/log links.")
  }
  sigma <- as.numeric(stats::fitted(model, what = "sigma"))
  if (any(!is.finite(sigma)) || any(sigma <= 0) || diff(range(sigma)) > 1e-8) {
    stop("The focal prediction adapter requires constant positive dispersion.")
  }
  fixed_formula <- stats::update.formula(
    stats::formula(model, what = "mu"), ~ . - random(plot_id) - random(year)
  )
  terms <- stats::delete.response(stats::terms(fixed_formula))
  contrasts <- attr(model$mu.x, "contrasts")
  contrasts <- contrasts[intersect(names(contrasts), all.vars(terms))]
  x <- stats::model.matrix(terms, newdata, contrasts.arg = contrasts)
  beta <- stats::coef(model, what = "mu")
  beta <- beta[!is.na(beta)]
  if (!all(names(beta) %in% colnames(x)) || nrow(x) != nrow(newdata)) {
    stop("Cannot reconstruct the focal fixed-effect prediction matrix.")
  }
  eta <- as.numeric(x[, names(beta), drop = FALSE] %*% beta) + eta_shift
  bd <- pmin(4, newdata$par_rich)
  if (any(!is.finite(bd)) || any(bd < 1 | bd %% 1 != 0)) {
    stop("Focal predictions require integer total richness and valid bounds.")
  }
  moments <- vapply(seq_along(eta), function(i) {
    support <- seq.int(0L, bd[[i]])
    log_p <- gamlss.dist::dDBI(
      support, mu = stats::plogis(eta[[i]]), sigma = sigma[[1L]],
      bd = bd[[i]], log = TRUE
    )
    p <- exp(log_p - max(log_p))
    p <- p / sum(p)
    mean <- sum(support * p)
    c(mean = mean, variance = sum((support - mean)^2 * p))
  }, numeric(2))
  if (any(!is.finite(moments))) stop("Non-finite focal DBI moments.")
  data.frame(eta = eta, mean = moments[1L, ], variance = moments[2L, ], bound = bd)
}

# Provide focal DBI predictions in the format used by SEM standardization.
# Details: Delegates to the exact-moment calculation and deliberately excludes
# realized random effects so predictions match all other pathway summaries.
# Returns: A vector on either the link scale or expected-count scale.
sem_predict_focal_dbi <- function(model, newdata, type, include_random = FALSE) {
  if (isTRUE(include_random)) stop("The focal SEM adapter excludes realized random effects.")
  moments <- sem_focal_dbi_moments(model, newdata)
  if (type == "link") moments$eta else moments$mean
}

# Convert saved focal-bootstrap draws into SEM coefficient inference.
# Details: Verifies the bootstrap against the fitted model and ordered data,
# then calculates coefficient or reference-condition slope percentile intervals.
# It does not turn percentile intervals into Wald tests or bootstrap P values.
# Returns: A coefficient table, inference callback, and interval-status label.
sem_focal_dbi_inference_adapter <- function(model, data, bootstrap) {
  beta <- stats::coef(model, what = "mu")
  beta <- beta[!is.na(beta)]
  reference_beta <- bootstrap$generating_parameters$beta
  same_data <- setequal(names(data), names(bootstrap$data)) &&
    isTRUE(all.equal(data[, names(bootstrap$data), drop = FALSE],
                     bootstrap$data, check.attributes = FALSE))
  if (!same_data || !identical(names(beta), names(reference_beta)) ||
      !isTRUE(all.equal(unname(beta), unname(reference_beta), tolerance = 1e-8)) ||
      !identical(deparse(stats::formula(model, what = "mu")),
                 deparse(stats::formula(bootstrap$reference_model, what = "mu"))) ||
      !isTRUE(all.equal(as.numeric(stats::fitted(model, what = "sigma")),
                       as.numeric(stats::fitted(bootstrap$reference_model, what = "sigma"))))) {
    stop("The focal bootstrap does not match this model and its ordered data.")
  }
  accepted <- bootstrap$audit$accepted
  draws <- bootstrap$draws[accepted, names(beta), drop = FALSE]
  if (nrow(draws) < 2L || any(!is.finite(draws))) {
    stop("At least two finite accepted focal bootstrap draws are required.")
  }
  sufficient <- nrow(draws) >= bootstrap$settings$required_successes
  confidence <- bootstrap$settings$confidence
  summarize <- function(gradient) {
    values <- as.numeric(draws %*% gradient)
    limits <- if (sufficient) stats::quantile(
      values, c((1 - confidence) / 2, (1 + confidence) / 2), names = FALSE
    ) else c(NA_real_, NA_real_)
    c(estimate = sum(beta * gradient), se = stats::sd(values),
      lower = limits[1L], upper = limits[2L])
  }
  raw <- t(vapply(seq_along(beta), function(i) {
    gradient <- numeric(length(beta))
    gradient[[i]] <- 1
    summarize(gradient)
  }, numeric(4)))
  coefficient_table <- cbind(
    Estimate = raw[, "estimate"], `Std. Error` = raw[, "se"],
    `z value` = NA_real_, `Pr(>|z|)` = NA_real_,
    Conf.Low = raw[, "lower"], Conf.High = raw[, "upper"]
  )
  rownames(coefficient_table) <- names(beta)
  inference <- function(model, term, term_type, reference_values) {
    if (!term %in% names(beta)) stop("Unknown focal coefficient: ", term)
    gradient <- stats::setNames(numeric(length(beta)), names(beta))
    gradient[[term]] <- 1
    conditioning <- character()
    if (term_type == "direct_path") {
      for (interaction in names(beta)[grepl(":", names(beta), fixed = TRUE)]) {
        parts <- strsplit(interaction, ":", fixed = TRUE)[[1L]]
        if (!term %in% parts) next
        moderators <- setdiff(parts, term)
        if (!all(moderators %in% names(reference_values))) stop("Missing moderator reference.")
        gradient[[interaction]] <- prod(reference_values[moderators])
        conditioning <- c(conditioning, paste0(
          moderators, "=", format(reference_values[moderators], digits = 8)
        ))
      }
    }
    result <- summarize(gradient)
    tibble::tibble(
      inference_estimand = if (term_type == "direct_path") {
        "reference_location_logit_slope"
      } else "location_logit_interaction_coefficient",
      inference_method = paste(
        "Hierarchical parametric bootstrap; percentile intervals; P values not estimated;",
        nrow(draws), "accepted of", bootstrap$settings$times
      ),
      reference_conditioning = if (length(conditioning)) {
        paste(unique(conditioning), collapse = "; ")
      } else if (term_type == "direct_path") "no interacting moderators" else NA_character_,
      inference_estimate = unname(result[["estimate"]]),
      inference_standard_error = unname(result[["se"]]),
      inference_conf_low = unname(result[["lower"]]),
      inference_conf_high = unname(result[["upper"]]),
      inference_test_statistic = NA_real_, inference_p_value = NA_real_
    )
  }
  list(coefficients = coefficient_table, inference = inference,
       interval_status = if (sufficient) "available_model_based" else "insufficient_refits")
}

# Build the feasibility rules used for network-path contrasts.
# Details: Uses integer contrasts for total-to-focal richness and restricts
# part-whole richness perturbations to biologically possible combinations.
# Returns: A callback that supplies feasible low/high data and evaluation rows.
sem_network_evaluation <- function(response) {
  force(response)
  function(data, predictor, low, high, mode) {
    rows <- seq_len(nrow(data))
    step <- NULL
    if (mode == "slope" && response == "focal_rich" && predictor == "par_rich") {
      low$par_rich <- pmax(min(data$par_rich), data$par_rich - 1)
      high$par_rich <- pmin(max(data$par_rich), data$par_rich + 1)
      step <- 1
    }
    # Observed focal richness is a predictor only in the downstream equations.
    if (response %in% c("cell_prate", "niche_overlap_hl") &&
        predictor %in% c("par_rich", "focal_rich")) {
      lower <- if (predictor == "par_rich") pmax(2, data$focal_rich) else 0
      upper <- if (predictor == "focal_rich") pmin(4, data$par_rich) else Inf
      if (mode == "slope") {
        low[[predictor]] <- pmax(low[[predictor]], lower)
        high[[predictor]] <- pmin(high[[predictor]], upper)
      } else {
        rows <- which(low[[predictor]] >= lower & high[[predictor]] <= upper)
        low <- low[rows, , drop = FALSE]
        high <- high[rows, , drop = FALSE]
      }
    }
    list(low = low, high = high, rows = rows, step = step)
  }
}

# Calculate family-aware marginal and conditional R2 for network components.
# Details: Uses exact DBI moments for focal richness, the complete zero-inflated
# beta mixture variance for niche overlap, and the shared R2 method otherwise.
# Returns: One-row component R2 metadata and variance decomposition.
sem_network_component_r2 <- function(model, data) {
  if (inherits(model, "gamlss")) {
    moments <- sem_focal_dbi_moments(model, data)
    h <- 1e-5
    derivative <- (sem_focal_dbi_moments(model, data, h)$mean -
                     sem_focal_dbi_moments(model, data, -h)$mean) / (2 * h)
    if (any(!is.finite(derivative)) || any(derivative <= 0)) {
      stop("Invalid exact-mean derivative for DBI R2.")
    }
    random_variance <- sum(vapply(1:2, function(j) {
      gamlss::getSmo(model, "mu", which = j)$sigb^2
    }, numeric(1)))
    fixed_variance <- stats::var(moments$eta)
    distribution_variance <- mean(moments$variance / derivative^2)
    total <- fixed_variance + random_variance + distribution_variance
    result <- data.frame(
      family = "DBI", link = model$mu.link, n = nrow(data), approximation = "delta",
      dispersion_random_effects = FALSE,
      variance_fixed = fixed_variance, variance_random = random_variance,
      variance_distribution = distribution_variance,
      r2_marginal = fixed_variance / total,
      r2_conditional = (fixed_variance + random_variance) / total
    )
    result$r2_definition <- paste(
      "Approximate Nakagawa decomposition on DBI location-logit scale;",
      "exact count moments and mean derivative; conditional on observed bounds"
    )
  } else if (length(glmmTMB::fixef(model)$zi)) {
    if (stats::family(model)$family != "beta" ||
        !identical(names(glmmTMB::fixef(model)$zi), "(Intercept)") ||
        length(glmmTMB::VarCorr(model)$zi)) {
      stop("The network mixture R2 adapter requires intercept-only zero-inflated beta.")
    }
    pi_zero <- stats::plogis(unname(glmmTMB::fixef(model)$zi[[1L]]))
    variance_function <- function(model, data, mu, eta, dispersion, trials) {
      positive_mean <- stats::plogis(eta)
      (1 - pi_zero) * positive_mean * (1 - positive_mean) / (dispersion + 1) +
        pi_zero * (1 - pi_zero) * positive_mean^2
    }
    derivative_function <- function(model, data, eta, mu) {
      positive_mean <- stats::plogis(eta)
      (1 - pi_zero) * positive_mean * (1 - positive_mean)
    }
    result <- sem_component_nakagawa_r2(
      model, data, approximation = "delta", variance_function = variance_function,
      mean_derivative_function = derivative_function
    )
    result$r2_definition <- paste(
      "Approximate Nakagawa decomposition on positive-overlap logit scale;",
      "full mixture variance and overall-mean derivative; constant zero probability"
    )
  } else {
    result <- sem_component_nakagawa_r2(model, data)
    result$r2_definition <- "Established family-aware Nakagawa link-scale decomposition"
  }
  if (any(!is.finite(unlist(result[c("r2_marginal", "r2_conditional")])))) {
    stop("Non-finite network component R2.")
  }
  result
}

# Summarize the network SEM's main ecological paths.
# Details: Extracts the directed paths from component term tables and attaches
# matching age interactions plus any additional moderation terms.
# Returns: One reporting row per main ecological path.
network_sem_path_table <- function(component_results) {
  terms <- dplyr::bind_rows(lapply(component_results, `[[`, "term_results"))
  paths <- terms[terms$term_role == "main_ecological_path", , drop = FALSE]
  age_rows <- terms[terms$term_role == "age_moderation", , drop = FALSE]
  age_rows$main_predictor <- sub("(^sc_fa:|:sc_fa$)", "", age_rows$predictor)
  age_fields <- c("predictor", "estimate", "std_error", "conf_low", "conf_high", "p_value")
  for (field in age_fields) {
    column <- if (field == "predictor") "age_interaction_term" else paste0("age_interaction_", field)
    paths[[column]] <- age_rows[[field]][match(
      paste(paths$model_id, paths$predictor), paste(age_rows$model_id, age_rows$main_predictor)
    )]
  }
  paths$moderation_terms <- vapply(seq_len(nrow(paths)), function(i) {
    candidates <- terms$predictor[terms$model_id == paths$model_id[[i]] &
                                  terms$term_role %in% c("age_moderation", "ecological_interaction")]
    paste(candidates[vapply(strsplit(candidates, ":", fixed = TRUE),
                            function(x) paths$predictor[[i]] %in% x, logical(1))], collapse = ";")
  }, character(1))
  paths
}

# Save a complete network-SEM handoff and its publication-facing tables.
# Details: Validates required objects, flattens list columns for CSV, and exports
# graph and missing-path results when the handoff contains a reduced graph.
# Returns: The output directory invisibly after all files have been written.
export_network_sem <- function(
    handoff,
    output_dir,
    file_stem = "network_sem") {
  .sem_assert_scalar_character(output_dir, "output_dir")
  .sem_assert_scalar_character(file_stem, "file_stem")
  required <- c("primary_sem_models", "primary_sem_term_results", "path_table",
                "component_r2_table", "node_reference_table")
  if (!all(required %in% names(handoff))) stop("Incomplete network SEM handoff.")
  columns <- c(
    "component", "model_id", "data_id", "response", "predictor", "term_role", "native_n",
    "family", "link", "inference_estimand", "inference_method", "reference_conditioning",
    "estimate", "std_error", "conf_low", "conf_high", "test_statistic", "p_value",
    "edge_id", "source_node", "destination_node", "age_moderated",
    "standardized_response_effect", "standardized_oe_effect", "standardized_relevant_range_effect",
    "response_estimand", "oe_estimand", "effect_method", "rr_n_evaluation", "rr_n_excluded",
    "uncertainty_status", "standardization_status"
  )
  tables <- list(
    terms = handoff$primary_sem_term_results[, columns, drop = FALSE],
    paths = handoff$path_table, component_r2 = handoff$component_r2_table,
    node_reference = handoff$node_reference_table
  )
  if (!is.null(handoff$graph)) {
    stopifnot(all(c("basis_set", "missing_path_results", "global_fit") %in% names(handoff)))
    tables <- c(tables, list(
      directed_edges = handoff$graph$directed_edges,
      nondirectional_associations = handoff$graph$nondirectional_associations,
      basis_set = handoff$basis_set, missing_paths = handoff$missing_path_results,
      global_fit = handoff$global_fit
    ))
  }
  # Flatten conditioning sets and named reference values for CSV; retain the
  # original list-columns and full fitted-model audit in the RDS handoff.
  for (name in names(tables)) {
    list_columns <- vapply(tables[[name]], is.list, logical(1))
    tables[[name]][list_columns] <- lapply(tables[[name]][list_columns], function(column) {
      vapply(column, function(value) {
        if (!is.null(names(value))) value <- paste0(names(value), "=", unlist(value))
        paste(value, collapse = ";")
      }, character(1))
    })
  }
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (!dir.exists(output_dir)) stop("Could not create the network output directory.")
  saveRDS(handoff, file.path(output_dir, paste0(file_stem, "_handoff.rds")))
  for (name in names(tables)) utils::write.csv(
    tables[[name]],
    file.path(output_dir, paste0(file_stem, "_", name, ".csv")),
    row.names = FALSE, na = "NA"
  )
  invisible(output_dir)
}

# -----------------------------------------------------------------------------
# Parsimonious-SEM reduction, testing, and export
# -----------------------------------------------------------------------------
# These helpers apply the recorded reductions, rebuild the pathway graph,
# test its missing paths, and prepare the publication-facing results.

# Compare one manually reduced component with its full parent. The function
# reports nested-model evidence, model validity, reference-conditioned
# coefficient stability, and optional grouped predictive R2 on shared folds.
# Details: Requires matched native samples and families, then compares AIC/LRT,
# retained-coefficient stability, and optionally shared-fold predictive R2.
# Returns: A list of fit tables, LRT, stability results, and CV details.
compare_component_reduction <- function(
    full_model,
    reduced_model,
    data,
    reference_values,
    focal_terms = NULL,
    run_cv = FALSE,
    folds = 5L,
    seed = 1201L,
    group = "plot_id",
    strata = "site",
    update_arguments = list()) {
  models <- list(full = full_model, reduced = reduced_model)

  if (any(vapply(models, stats::nobs, numeric(1)) != nrow(data)) ||
      !identical(stats::family(full_model)[c("family", "link")],
                 stats::family(reduced_model)[c("family", "link")])) {
    stop("Models must use the same native sample, response family, and link.")
  }
  model_data_calls <- lapply(models, function(model) model$call$data)
  if (!identical(model_data_calls$full, model_data_calls$reduced)) {
    stop(
      "Models must be fitted using the same data-object expression. ",
      "Refit the parent on `component_data` before deriving the reduced model."
    )
  }

  model_comparison <- dplyr::bind_rows(lapply(names(models), function(model_id) {
    fit <- models[[model_id]]
    model_loglik <- stats::logLik(fit)
    tibble::tibble(
      model = model_id,
      converged = isTRUE(fit$fit$convergence == 0L),
      pd_hessian = isTRUE(fit$sdr$pdHess),
      n = stats::nobs(fit),
      parameters = attr(model_loglik, "df"),
      logLik = as.numeric(model_loglik),
      AIC = stats::AIC(fit)
    )
  })) |>
    dplyr::mutate(delta_AIC = AIC - min(AIC))

  likelihood_ratio_test <- stats::anova(
    reduced_model,
    full_model,
    test = "Chisq"
  )

  full_coefficients <- names(glmmTMB::fixef(full_model)$cond)
  reduced_coefficients <- names(glmmTMB::fixef(reduced_model)$cond)
  common_terms <- setdiff(
    intersect(full_coefficients, reduced_coefficients),
    "(Intercept)"
  )
  if (is.null(focal_terms)) focal_terms <- common_terms
  if (!all(focal_terms %in% common_terms)) {
    stop("Every `focal_terms` entry must remain in both models.")
  }

  coefficient_stability <- dplyr::bind_rows(lapply(focal_terms, function(term) {
    term_type <- if (grepl(":", term, fixed = TRUE)) {
      "interaction"
    } else {
      "direct_path"
    }
    full_result <- sem_term_wald_test(
      full_model, term, term_type, reference_values
    )
    reduced_result <- sem_term_wald_test(
      reduced_model, term, term_type, reference_values
    )
    tibble::tibble(
      term = term,
      full_estimate = full_result$inference_estimate,
      reduced_estimate = reduced_result$inference_estimate,
      estimate_change = reduced_estimate - full_estimate,
      change_in_full_se = estimate_change /
        full_result$inference_standard_error,
      full_p_value = full_result$inference_p_value,
      reduced_p_value = reduced_result$inference_p_value
    )
  }))

  cv_results <- NULL
  if (isTRUE(run_cv)) {
    fold_id <- make_grouped_sem_folds(
      data,
      group = group,
      strata = strata,
      folds = folds,
      seed = seed
    )
    cv_results <- lapply(models, function(model) {
      sem_component_predictive_r2(
        model = model,
        data = data,
        group = group,
        strata = strata,
        fold_id = fold_id,
        update_arguments = update_arguments,
        on_refit_failure = "record"
      )
    })
    cv_summary <- dplyr::bind_rows(
      lapply(names(cv_results), function(model_id) {
        dplyr::mutate(cv_results[[model_id]]$summary, model = model_id, .before = 1)
      })
    )
    if (all(cv_summary$cv_status == "available")) {
      cv_summary <- cv_summary |>
        dplyr::mutate(
          delta_r2_response_cv = r2_response_cv -
            r2_response_cv[model == "full"]
        )
    } else {
      cv_summary$delta_r2_response_cv <- NA_real_
    }
  } else {
    cv_summary <- tibble::tibble()
  }

  list(
    model_comparison = model_comparison,
    likelihood_ratio_test = likelihood_ratio_test,
    coefficient_stability = coefficient_stability,
    cv_summary = cv_summary,
    cv_details = cv_results
  )
}

# Replay a documented sequence of manual term-removal decisions. Each candidate
# is compared with the last accepted model; rejected candidates do not alter the
# parent model, and a "review" row is evaluated but left pending.
# Details: Refits the maximal parent on native data and evaluates removals in
# order; accepted steps update the parent, while rejected steps do not.
# Returns: Final and candidate models, comparisons, term results, and decision log.
run_component_reduction_plan <- function(
    component,
    maximal_model,
    data,
    plan,
    reference_values,
    control) {
  required_columns <- c(
    "step", "removed_term", "decision", "run_cv", "rationale"
  )
  missing_columns <- setdiff(required_columns, names(plan))
  if (length(missing_columns) > 0L) {
    stop(
      "The reduction plan is missing: ",
      paste(missing_columns, collapse = ", "),
      "."
    )
  }

  plan <- plan |>
    dplyr::select(dplyr::all_of(required_columns)) |>
    dplyr::arrange(step)

  if (nrow(plan) == 0L ||
      anyDuplicated(plan$step) ||
      any(is.na(plan$step)) ||
      !all(plan$step == seq_len(nrow(plan)))) {
    stop("`step` must contain the unique consecutive integers 1, 2, ... .")
  }
  if (any(!plan$decision %in% c("accept", "reject", "review"))) {
    stop("`decision` must be \"accept\", \"reject\", or \"review\".")
  }
  review_steps <- which(plan$decision == "review")
  if (length(review_steps) > 1L ||
      (length(review_steps) == 1L && review_steps != nrow(plan))) {
    stop("At most one `review` row is allowed, and it must be the last row.")
  }
  completed_steps <- plan$decision %in% c("accept", "reject")
  if (any(completed_steps & !nzchar(trimws(plan$rationale)))) {
    stop("Every accepted or rejected step requires a rationale.")
  }

  current_model <- update(maximal_model, data = data, control = control)
  current_label <- "maximal"
  model_history <- list(maximal = current_model)
  comparison_history <- list()
  log_rows <- vector("list", nrow(plan))
  pending_candidate <- NULL

  for (row_id in seq_len(nrow(plan))) {
    plan_row <- plan[row_id, , drop = FALSE]
    removed_term <- plan_row$removed_term[[1]]
    parent_label <- current_label
    candidate_label <- paste0("step_", plan_row$step[[1]])

    if (!removed_term %in% names(glmmTMB::fixef(current_model)$cond)) {
      stop(
        "Term `", removed_term, "` is not present in the conditional model ",
        "used as the parent for step ", plan_row$step[[1]], "."
      )
    }

    reduction_formula <- stats::as.formula(
      paste(". ~ . -", removed_term),
      env = environment(stats::formula(current_model))
    )
    candidate_model <- update(
      current_model,
      reduction_formula,
      data = data,
      control = control
    )
    comparison <- compare_component_reduction(
      full_model = current_model,
      reduced_model = candidate_model,
      data = data,
      reference_values = reference_values,
      run_cv = plan_row$run_cv[[1]],
      update_arguments = list(control = control)
    )

    parent_summary <- comparison$model_comparison |>
      dplyr::filter(model == "full")
    candidate_summary <- comparison$model_comparison |>
      dplyr::filter(model == "reduced")
    lrt_row <- comparison$likelihood_ratio_test[
      nrow(comparison$likelihood_ratio_test), , drop = FALSE
    ]
    maximum_change <- if (nrow(comparison$coefficient_stability) > 0L) {
      max(abs(comparison$coefficient_stability$change_in_full_se), na.rm = TRUE)
    } else {
      NA_real_
    }
    if (!is.finite(maximum_change)) maximum_change <- NA_real_

    parent_cv_r2 <- NA_real_
    candidate_cv_r2 <- NA_real_
    cv_status <- "not run"
    if (nrow(comparison$cv_summary) > 0L) {
      parent_cv <- comparison$cv_summary |>
        dplyr::filter(model == "full")
      candidate_cv <- comparison$cv_summary |>
        dplyr::filter(model == "reduced")
      parent_cv_r2 <- parent_cv$r2_response_cv[[1]]
      candidate_cv_r2 <- candidate_cv$r2_response_cv[[1]]
      cv_status <- paste(unique(comparison$cv_summary$cv_status), collapse = "; ")
    }

    log_rows[[row_id]] <- tibble::tibble(
      component = component,
      step = plan_row$step[[1]],
      parent_model = parent_label,
      candidate_model = candidate_label,
      removed_term = removed_term,
      decision = plan_row$decision[[1]],
      converged = candidate_summary$converged[[1]],
      pd_hessian = candidate_summary$pd_hessian[[1]],
      parent_AIC = parent_summary$AIC[[1]],
      candidate_AIC = candidate_summary$AIC[[1]],
      delta_AIC = candidate_AIC - parent_AIC,
      LRT_chisq = lrt_row[["Chisq"]],
      LRT_df = lrt_row[["Chi Df"]],
      LRT_p_value = lrt_row[["Pr(>Chisq)"]],
      max_abs_change_in_parent_se = maximum_change,
      cv_status = cv_status,
      parent_cv_r2 = parent_cv_r2,
      candidate_cv_r2 = candidate_cv_r2,
      delta_cv_r2 = candidate_cv_r2 - parent_cv_r2,
      rationale = plan_row$rationale[[1]]
    )
    comparison_history[[candidate_label]] <- comparison
    model_history[[candidate_label]] <- candidate_model

    if (plan_row$decision[[1]] == "accept") {
      current_model <- candidate_model
      current_label <- candidate_label
    } else if (plan_row$decision[[1]] == "review") {
      pending_candidate <- candidate_model
      break
    }
  }

  current_terms <- setdiff(
    names(glmmTMB::fixef(current_model)$cond),
    "(Intercept)"
  )
  current_term_results <- dplyr::bind_rows(lapply(current_terms, function(term) {
    term_type <- if (grepl(":", term, fixed = TRUE)) {
      "interaction"
    } else {
      "direct_path"
    }
    dplyr::bind_cols(
      tibble::tibble(term = term, term_type = term_type),
      sem_term_wald_test(
        model = current_model,
        term = term,
        term_type = term_type,
        reference_values = reference_values
      )
    )
  })) |>
    dplyr::transmute(
      term,
      estimate = inference_estimate,
      p_value = inference_p_value
    )

  list(
    current_model = current_model,
    current_term_results = current_term_results,
    pending_candidate = pending_candidate,
    model_history = model_history,
    comparison_history = comparison_history,
    removal_log = dplyr::bind_rows(log_rows)
  )
}

# Refit a finalized component, add any basis-set conditioning terms missing
# from that response model, and then test its declared focal term block against
# the conditioned null. Results remain provisional until the complete
# parsimonious SEM and final basis set exist.
# Details: Adds declared conditioning terms to the finalized response model,
# then compares that null with a focal-path add-back on the same native rows.
# Returns: A list with the test summary, specifications, refits, and comparisons.
audit_component_missing_paths <- function(
    component,
    final_model,
    data,
    path_spec,
    reference_values,
    control) {
  required_columns <- c("path_id", "focal_term", "added_terms")
  missing_columns <- setdiff(required_columns, names(path_spec))
  if (length(missing_columns) > 0L) {
    stop(
      "The missing-path specification is missing: ",
      paste(missing_columns, collapse = ", "),
      "."
    )
  }
  if (nrow(path_spec) == 0L || anyDuplicated(path_spec$path_id)) {
    stop("`path_spec` must contain at least one uniquely named path.")
  }
  if (!is.list(path_spec$added_terms)) {
    stop("`added_terms` must be a list-column of character vectors.")
  }
  if (!"conditioning_terms" %in% names(path_spec)) {
    path_spec$conditioning_terms <- rep(
      list(character()),
      nrow(path_spec)
    )
  }
  if (!is.list(path_spec$conditioning_terms)) {
    stop("`conditioning_terms` must be a list-column of character vectors.")
  }

  final_model <- update(final_model, data = data, control = control)
  final_coefficients <- names(glmmTMB::fixef(final_model)$cond)
  comparison_history <- list()
  conditioned_null_models <- list()
  addback_models <- list()
  summary_rows <- vector("list", nrow(path_spec))

  for (row_id in seq_len(nrow(path_spec))) {
    path_id <- path_spec$path_id[[row_id]]
    focal_term <- path_spec$focal_term[[row_id]]
    added_terms <- unlist(path_spec$added_terms[[row_id]], use.names = FALSE)
    conditioning_terms <- unlist(
      path_spec$conditioning_terms[[row_id]],
      use.names = FALSE
    )

    if (!is.character(path_id) || length(path_id) != 1L ||
        is.na(path_id) || !nzchar(path_id) ||
        !is.character(focal_term) || length(focal_term) != 1L ||
        is.na(focal_term) || !nzchar(focal_term) ||
        length(added_terms) == 0L || anyNA(added_terms) ||
        any(!nzchar(added_terms)) ||
        !is.character(conditioning_terms) || anyNA(conditioning_terms) ||
        any(!nzchar(conditioning_terms))) {
      stop("Each missing-path specification must use non-empty term names.")
    }
    if (!focal_term %in% added_terms) {
      stop("The focal term for `", path_id, "` must occur in `added_terms`.")
    }
    if (anyDuplicated(c(added_terms, conditioning_terms))) {
      stop(
        "Add-back and conditioning terms must be unique for `",
        path_id, "`."
      )
    }
    new_conditioning_terms <- setdiff(
      conditioning_terms,
      final_coefficients
    )
    conditioned_null_model <- final_model
    if (length(new_conditioning_terms) > 0L) {
      conditioning_formula <- stats::as.formula(
        paste(
          ". ~ . +",
          paste(new_conditioning_terms, collapse = " + ")
        ),
        env = environment(stats::formula(final_model))
      )
      conditioned_null_model <- update(
        final_model,
        conditioning_formula,
        data = data,
        control = control
      )
    }
    null_coefficients <- names(
      glmmTMB::fixef(conditioned_null_model)$cond
    )
    already_present <- intersect(added_terms, null_coefficients)
    if (length(already_present) > 0L) {
      stop(
        "The conditioned null model already contains the add-back term(s) for `",
        path_id, "`: ", paste(already_present, collapse = ", "), "."
      )
    }

    addback_formula <- stats::as.formula(
      paste(". ~ . +", paste(added_terms, collapse = " + ")),
      env = environment(stats::formula(final_model))
    )
    addback_model <- update(
      conditioned_null_model,
      addback_formula,
      data = data,
      control = control
    )
    comparison <- compare_component_reduction(
      full_model = addback_model,
      reduced_model = conditioned_null_model,
      data = data,
      reference_values = reference_values,
      run_cv = FALSE,
      update_arguments = list(control = control)
    )
    focal_test <- sem_term_wald_test(
      model = addback_model,
      term = focal_term,
      term_type = "direct_path",
      reference_values = reference_values
    )

    null_summary <- comparison$model_comparison |>
      dplyr::filter(model == "reduced")
    addback_summary <- comparison$model_comparison |>
      dplyr::filter(model == "full")
    lrt_row <- comparison$likelihood_ratio_test[
      nrow(comparison$likelihood_ratio_test), , drop = FALSE
    ]

    summary_rows[[row_id]] <- tibble::tibble(
      component = component,
      path_id = path_id,
      focal_term = focal_term,
      added_terms = paste(added_terms, collapse = " + "),
      conditioning_terms = paste(conditioning_terms, collapse = " + "),
      n = stats::nobs(final_model),
      final_converged = isTRUE(final_model$fit$convergence == 0L),
      final_pd_hessian = isTRUE(final_model$sdr$pdHess),
      conditioned_null_converged = null_summary$converged[[1]],
      conditioned_null_pd_hessian = null_summary$pd_hessian[[1]],
      addback_converged = addback_summary$converged[[1]],
      addback_pd_hessian = addback_summary$pd_hessian[[1]],
      final_AIC = stats::AIC(final_model),
      conditioned_null_AIC = null_summary$AIC[[1]],
      addback_AIC = addback_summary$AIC[[1]],
      delta_AIC_addback = addback_AIC - conditioned_null_AIC,
      LRT_chisq = lrt_row[["Chisq"]],
      LRT_df = lrt_row[["Chi Df"]],
      LRT_p_value = lrt_row[["Pr(>Chisq)"]],
      reference_estimate = focal_test$inference_estimate[[1]],
      reference_std_error = focal_test$inference_standard_error[[1]],
      reference_p_value = focal_test$inference_p_value[[1]],
      provisional = TRUE
    )
    comparison_history[[path_id]] <- comparison
    conditioned_null_models[[path_id]] <- conditioned_null_model
    addback_models[[path_id]] <- addback_model
  }

  list(
    summary = dplyr::bind_rows(summary_rows),
    final_model = final_model,
    path_spec = path_spec,
    comparison_history = comparison_history,
    conditioned_null_models = conditioned_null_models,
    addback_models = addback_models
  )
}

# Assemble a reduced SEM graph from its finalized component models, including
# both directed ecological edges and explicitly declared nondirectional
# associations. The maximal directed-edge registry remains the naming and
# metadata template for refreshing edge presence and age-moderation status.
# Details: Reconciles retained coefficients with the maximal edge registry,
# verifies exactly the expected removals, and carries declared associations.
# Returns: A list of directed edges, associations, removed edges, and model terms.
assemble_reduced_sem_graph <- function(
    models,
    maximal_edge_spec,
    expected_removed_edges,
    nondirectional_associations = NULL,
    age_variable = "sc_fa") {
  required_columns <- c(
    "edge_id", "model_id", "predictor", "source_node", "destination_node"
  )
  missing_columns <- setdiff(required_columns, names(maximal_edge_spec))
  if (length(missing_columns) > 0L) {
    stop(
      "`maximal_edge_spec` lacks: ",
      paste(missing_columns, collapse = ", "),
      "."
    )
  }
  if (is.null(names(models)) || any(!nzchar(names(models))) ||
      anyDuplicated(names(models))) {
    stop("`models` must be a uniquely named list.")
  }
  required_model_ids <- unique(maximal_edge_spec$model_id)
  if (!setequal(names(models), required_model_ids)) {
    stop(
      "The reduced-model registry must contain exactly these model IDs: ",
      paste(required_model_ids, collapse = ", "),
      "."
    )
  }
  if (anyDuplicated(maximal_edge_spec$edge_id)) {
    stop("The maximal edge registry contains duplicated edge IDs.")
  }
  if (!is.character(expected_removed_edges) ||
      anyNA(expected_removed_edges) ||
      any(!nzchar(expected_removed_edges)) ||
      anyDuplicated(expected_removed_edges)) {
    stop("`expected_removed_edges` must contain unique non-empty IDs.")
  }
  unknown_expected_edges <- setdiff(
    expected_removed_edges,
    maximal_edge_spec$edge_id
  )
  if (length(unknown_expected_edges) > 0L) {
    stop(
      "Unknown expected removed edge(s): ",
      paste(unknown_expected_edges, collapse = ", "),
      "."
    )
  }
  .sem_assert_scalar_character(age_variable, "age_variable")

  valid_models <- vapply(models, function(model) {
    if (inherits(model, "gamlss") && identical(model$family[[1L]], "DBI")) {
      coefficients <- stats::coef(model, what = "mu")
      fixed <- coefficients[!startsWith(names(coefficients), "random(")]
      return(isTRUE(model$converged) && length(fixed) > 0L && all(is.finite(fixed)))
    }
    inherits(model, "glmmTMB") &&
      isTRUE(model$fit$convergence == 0L) &&
      isTRUE(model$sdr$pdHess) &&
      all(is.finite(glmmTMB::fixef(model)$cond))
  }, logical(1))
  if (any(!valid_models)) {
    stop(
      "Invalid reduced component model(s): ",
      paste(names(valid_models)[!valid_models], collapse = ", "),
      "."
    )
  }

  model_terms <- lapply(
    models,
    function(model) if (inherits(model, "gamlss")) {
      names(stats::coef(model, what = "mu"))
    } else names(glmmTMB::fixef(model)$cond)
  )
  retained <- mapply(
    function(model_id, predictor) {
      predictor %in% model_terms[[model_id]]
    },
    maximal_edge_spec$model_id,
    maximal_edge_spec$predictor,
    USE.NAMES = FALSE
  )
  edge_spec <- maximal_edge_spec[retained, , drop = FALSE]
  edge_spec$age_moderated <- mapply(
    function(model_id, predictor) {
      possible_interactions <- c(
        paste0(predictor, ":", age_variable),
        paste0(age_variable, ":", predictor)
      )
      any(possible_interactions %in% model_terms[[model_id]])
    },
    edge_spec$model_id,
    edge_spec$predictor,
    USE.NAMES = FALSE
  )
  edge_spec$at <- lapply(
    edge_spec$age_moderated,
    function(is_moderated) {
      if (is_moderated) stats::setNames(list(0), age_variable) else list()
    }
  )
  removed_edges <- setdiff(
    maximal_edge_spec$edge_id,
    edge_spec$edge_id
  )
  if (!setequal(removed_edges, expected_removed_edges)) {
    stop(
      "The fitted reduced models do not reproduce the declared edge set. ",
      "Unexpected removals: ",
      paste(setdiff(removed_edges, expected_removed_edges), collapse = ", "),
      "; expected removals still present: ",
      paste(setdiff(expected_removed_edges, removed_edges), collapse = ", "),
      "."
    )
  }
  if (anyDuplicated(edge_spec$edge_id)) {
    stop("The reduced edge registry contains duplicated edge IDs.")
  }

  if (is.null(nondirectional_associations)) {
    nondirectional_associations <- tibble::tibble(
      association_id = character(),
      node_1 = character(),
      node_2 = character(),
      association_type = character(),
      native_n = numeric(),
      estimate = numeric(),
      conditioning = character()
    )
  } else {
    required_association_columns <- c(
      "association_id", "node_1", "node_2", "association_type",
      "native_n", "estimate", "conditioning"
    )
    missing_association_columns <- setdiff(
      required_association_columns,
      names(nondirectional_associations)
    )
    if (length(missing_association_columns) > 0L) {
      stop(
        "`nondirectional_associations` lacks: ",
        paste(missing_association_columns, collapse = ", "),
        "."
      )
    }
    association_text_columns <- c(
      "association_id", "node_1", "node_2",
      "association_type", "conditioning"
    )
    valid_associations <- nrow(nondirectional_associations) > 0L &&
      !anyDuplicated(nondirectional_associations$association_id) &&
      all(vapply(
        nondirectional_associations[association_text_columns],
        function(column) {
          is.character(column) && !anyNA(column) && all(nzchar(column))
        },
        logical(1)
      )) &&
      all(nondirectional_associations$node_1 !=
            nondirectional_associations$node_2) &&
      all(is.finite(nondirectional_associations$native_n)) &&
      all(nondirectional_associations$native_n > 0) &&
      all(is.finite(nondirectional_associations$estimate)) &&
      all(abs(nondirectional_associations$estimate) <= 1) &&
      all(nondirectional_associations$node_1 %in%
            c(maximal_edge_spec$source_node,
              maximal_edge_spec$destination_node)) &&
      all(nondirectional_associations$node_2 %in%
            c(maximal_edge_spec$source_node,
              maximal_edge_spec$destination_node)) &&
      !any(nondirectional_associations$association_id %in%
             maximal_edge_spec$edge_id)
    if (!isTRUE(valid_associations)) {
      stop("The nondirectional-association registry is invalid.")
    }
    nondirectional_associations <- tibble::as_tibble(
      nondirectional_associations
    )
  }

  list(
    directed_edges = tibble::as_tibble(edge_spec),
    nondirectional_associations = nondirectional_associations,
    edge_spec = tibble::as_tibble(edge_spec),
    removed_edges = removed_edges,
    model_terms = model_terms
  )
}

# Attach the declared ecological conditioning sets and response-model metadata
# to every missing edge in a finalized reduced graph. Conditioning sets remain
# an explicit scientific declaration rather than being inferred from row-level
# data with incompatible observational grains.
# Details: Joins every removed edge to its scientifically declared conditioning
# set and response-model metadata, confirming the path remains absent.
# Returns: A tibble with one directed-separation claim per removed edge.
derive_sem_basis_set <- function(
    maximal_edge_spec,
    removed_edges,
    conditioning_spec,
    models) {
  required_edge_columns <- c(
    "edge_id", "source_node", "destination_node", "model_id", "data_id",
    "predictor", "response"
  )
  missing_edge_columns <- setdiff(
    required_edge_columns,
    names(maximal_edge_spec)
  )
  if (length(missing_edge_columns) > 0L) {
    stop(
      "`maximal_edge_spec` lacks: ",
      paste(missing_edge_columns, collapse = ", "),
      "."
    )
  }
  required_conditioning_columns <- c(
    "edge_id", "conditioning_nodes", "conditioning_terms"
  )
  missing_conditioning_columns <- setdiff(
    required_conditioning_columns,
    names(conditioning_spec)
  )
  if (length(missing_conditioning_columns) > 0L) {
    stop(
      "`conditioning_spec` lacks: ",
      paste(missing_conditioning_columns, collapse = ", "),
      "."
    )
  }
  if (!is.list(conditioning_spec$conditioning_nodes) ||
      !is.list(conditioning_spec$conditioning_terms)) {
    stop(
      "`conditioning_nodes` and `conditioning_terms` must be list-columns."
    )
  }
  if (anyDuplicated(conditioning_spec$edge_id) ||
      !setequal(conditioning_spec$edge_id, removed_edges)) {
    stop("`conditioning_spec` must declare each removed edge exactly once.")
  }
  valid_conditioning <- vapply(
    conditioning_spec$conditioning_nodes,
    function(nodes) {
      is.character(nodes) && length(nodes) > 0L &&
        !anyNA(nodes) && all(nzchar(nodes)) && !anyDuplicated(nodes)
    },
    logical(1)
  )
  if (any(!valid_conditioning)) {
    stop("Every basis claim needs unique non-empty conditioning-node names.")
  }
  valid_conditioning_terms <- vapply(
    conditioning_spec$conditioning_terms,
    function(terms) {
      is.character(terms) && !anyNA(terms) &&
        all(nzchar(terms)) && !anyDuplicated(terms)
    },
    logical(1)
  )
  if (any(!valid_conditioning_terms)) {
    stop("Every basis claim needs unique conditioning-term names.")
  }

  basis_rows <- match(removed_edges, maximal_edge_spec$edge_id)
  if (anyNA(basis_rows)) stop("At least one removed edge is unknown.")
  basis_set <- maximal_edge_spec[basis_rows, , drop = FALSE]
  conditioning_rows <- match(
    basis_set$edge_id,
    conditioning_spec$edge_id
  )
  basis_set$conditioning_nodes <-
    conditioning_spec$conditioning_nodes[conditioning_rows]
  basis_set$conditioning_terms <-
    conditioning_spec$conditioning_terms[conditioning_rows]

  unknown_models <- setdiff(basis_set$model_id, names(models))
  if (length(unknown_models) > 0L) {
    stop(
      "The basis set refers to unknown model ID(s): ",
      paste(unknown_models, collapse = ", "),
      "."
    )
  }
  paths_still_present <- mapply(
    function(model_id, predictor) {
      predictor %in% names(glmmTMB::fixef(models[[model_id]])$cond)
    },
    basis_set$model_id,
    basis_set$predictor,
    USE.NAMES = FALSE
  )
  if (any(paths_still_present)) {
    stop(
      "A declared missing path remains in its finalized model: ",
      paste(basis_set$edge_id[paths_still_present], collapse = ", "),
      "."
    )
  }

  basis_set$native_n <- vapply(
    basis_set$model_id,
    function(model_id) stats::nobs(models[[model_id]]),
    numeric(1)
  )
  basis_set$family <- vapply(
    basis_set$model_id,
    function(model_id) stats::family(models[[model_id]])$family,
    character(1)
  )
  basis_set$link <- vapply(
    basis_set$model_id,
    function(model_id) stats::family(models[[model_id]])$link,
    character(1)
  )

  tibble::as_tibble(basis_set) |>
    dplyr::transmute(
      claim_id = edge_id,
      source_node,
      destination_node,
      model_id,
      data_id,
      predictor,
      response,
      conditioning_nodes,
      conditioning_terms,
      native_n,
      family,
      link
    )
}

# Refit and test every claim in a finalized multiscale basis set using the
# declared downstream response model and its native data. The function keeps
# the component-level audit objects while returning one validated results table.
# Details: Groups claims by response model, runs native-sample add-back tests,
# and requires valid 1-df comparisons with converged, positive-Hessian fits.
# Returns: A list containing claim specifications, full audits, and result rows.
audit_sem_basis_set <- function(
    basis_set,
    models,
    data_sets,
    model_data_ids,
    component_metadata,
    reference_values,
    control) {
  required_columns <- c(
    "claim_id", "model_id", "data_id", "predictor", "native_n",
    "family", "link", "conditioning_nodes", "conditioning_terms"
  )
  missing_columns <- setdiff(required_columns, names(basis_set))
  if (length(missing_columns) > 0L) {
    stop(
      "`basis_set` lacks: ",
      paste(missing_columns, collapse = ", "),
      "."
    )
  }
  if (nrow(basis_set) == 0L || anyDuplicated(basis_set$claim_id)) {
    stop("`basis_set` must contain at least one uniquely named claim.")
  }
  if (is.null(names(models)) || is.null(names(data_sets)) ||
      is.null(names(model_data_ids))) {
    stop("Model, dataset, and model-to-data registries must be named.")
  }
  unknown_models <- setdiff(basis_set$model_id, names(models))
  unknown_data <- setdiff(basis_set$data_id, names(data_sets))
  if (length(unknown_models) > 0L || length(unknown_data) > 0L) {
    stop("The basis set refers to an unknown model or native dataset.")
  }
  mapped_data <- unname(model_data_ids[basis_set$model_id])
  if (anyNA(mapped_data) || !identical(mapped_data, basis_set$data_id)) {
    stop("At least one basis claim does not use its model's native dataset.")
  }

  metadata_columns <- c("component", "source_model")
  if (!all(metadata_columns %in% names(component_metadata))) {
    stop("`component_metadata` lacks component-to-model identifiers.")
  }
  component_rows <- match(
    unique(basis_set$model_id),
    component_metadata$source_model
  )
  if (anyNA(component_rows)) {
    stop("At least one basis-set model lacks component metadata.")
  }
  component_by_model <- stats::setNames(
    component_metadata$component[component_rows],
    component_metadata$source_model[component_rows]
  )

  claims_by_model <- split(basis_set, basis_set$model_id)
  path_specs <- lapply(claims_by_model, function(claims) {
    tibble::tibble(
      path_id = claims$claim_id,
      focal_term = claims$predictor,
      added_terms = as.list(claims$predictor),
      conditioning_terms = claims$conditioning_terms
    )
  })
  audits <- lapply(names(path_specs), function(model_id) {
    audit_component_missing_paths(
      component = component_by_model[[model_id]],
      final_model = models[[model_id]],
      data = data_sets[[model_data_ids[[model_id]]]],
      path_spec = path_specs[[model_id]],
      reference_values = reference_values,
      control = control
    )
  })
  names(audits) <- names(path_specs)

  results <- lapply(audits, `[[`, "summary") |>
    dplyr::bind_rows() |>
    dplyr::mutate(provisional = FALSE) |>
    dplyr::left_join(
      basis_set |>
        dplyr::select(
          claim_id, data_id, conditioning_nodes,
          native_n, family, link
        ),
      by = c("path_id" = "claim_id")
    )

  valid_results <- nrow(results) == nrow(basis_set) &&
    setequal(results$path_id, basis_set$claim_id) &&
    all(results$n == results$native_n) &&
    all(results$final_converged) &&
    all(results$final_pd_hessian) &&
    all(results$conditioned_null_converged) &&
    all(results$conditioned_null_pd_hessian) &&
    all(results$addback_converged) &&
    all(results$addback_pd_hessian) &&
    all(results$LRT_df == 1) &&
    all(is.finite(results$LRT_p_value)) &&
    all(results$LRT_p_value > 0) &&
    all(results$LRT_p_value <= 1)
  if (!isTRUE(valid_results)) {
    stop("At least one final basis-set claim failed validation.")
  }

  list(
    path_specs = path_specs,
    audits = audits,
    results = results
  )
}

# Fisher's C for a validated set of directed-separation tests. This function
# deliberately does not calculate a global AIC or AICc for multiscale models.
# Details: Applies -2 times the sum of log P values with 2k degrees of freedom;
# it deliberately avoids a global AIC across incompatible observational grains.
# Returns: A one-row tibble with Fisher's C, df, P value, and interpretation.
sem_fisher_c <- function(
    missing_path_results,
    p_value_column = "LRT_p_value") {
  .sem_assert_scalar_character(p_value_column, "p_value_column")
  if (!p_value_column %in% names(missing_path_results)) {
    stop("Unknown P-value column: `", p_value_column, "`.")
  }
  p_values <- missing_path_results[[p_value_column]]
  if (length(p_values) == 0L || any(!is.finite(p_values)) ||
      any(p_values <= 0) || any(p_values > 1)) {
    stop("Fisher's C requires finite P values in the interval (0, 1].")
  }
  validation_columns <- c(
    "final_converged", "final_pd_hessian",
    "addback_converged", "addback_pd_hessian"
  )
  if (all(validation_columns %in% names(missing_path_results))) {
    fit_status <- unlist(
      missing_path_results[validation_columns],
      use.names = FALSE
    )
    if (anyNA(fit_status) || any(!fit_status)) {
      stop("Fisher's C cannot use a basis claim with an invalid fitted model.")
    }
  }

  fisher_c <- -2 * sum(log(p_values))
  degrees_freedom <- 2L * length(p_values)
  p_value <- stats::pchisq(
    fisher_c,
    df = degrees_freedom,
    lower.tail = FALSE
  )

  tibble::tibble(
    n_independence_claims = length(p_values),
    fisher_C = fisher_c,
    df = degrees_freedom,
    p_value = p_value,
    interpretation = if (p_value > 0.05) {
      "basis-set claims are consistent with the parsimonious graph"
    } else {
      "at least one basis-set claim is inconsistent with the parsimonious graph"
    },
    scope = paste(
      "Fisher's C from native-response add-back LRTs;",
      "no global AIC or AICc because components use different grains"
    )
  )
}

# Assemble one publication-facing row per retained directed path. Link-scale
# inference comes from the finalized component models; the three standardized
# point estimates reuse the node-anchored definitions above.
# Details: Combines reference-conditioned Wald inference, age-interaction
# information, and all three standardized effect definitions for retained edges.
# Returns: One publication-facing data-frame row per directed path.
build_sem_path_table <- function(
    edge_spec,
    models,
    data_sets,
    node_reference,
    reference_values,
    include_relevant_range = TRUE) {
  required <- c(
    "edge_id", "source_node", "destination_node", "model_id", "data_id",
    "predictor", "response", "age_moderated"
  )
  if (!all(required %in% names(edge_spec)) ||
      anyDuplicated(edge_spec$edge_id)) {
    stop("`edge_spec` lacks required columns or has duplicated edge IDs.")
  }

  standardized <- standardize_sem_edges(
    edge_spec = edge_spec,
    models = models,
    data_sets = data_sets,
    node_reference = node_reference,
    include_relevant_range = include_relevant_range,
    include_random = FALSE,
    on_error = "stop"
  )

  inference <- dplyr::bind_rows(lapply(seq_len(nrow(edge_spec)), function(i) {
    edge <- edge_spec[i, , drop = FALSE]
    model <- models[[edge$model_id[[1]]]]
    coefficients <- glmmTMB::fixef(model)$cond
    covariance <- as.matrix(stats::vcov(model)$cond)
    predictor <- edge$predictor[[1]]
    path_test <- sem_term_wald_test(
      model, predictor, "direct_path", reference_values
    )
    interactions <- names(coefficients)[vapply(
      strsplit(names(coefficients), ":", fixed = TRUE),
      function(parts) length(parts) > 1L && predictor %in% parts,
      logical(1)
    )]
    age_interaction <- interactions[vapply(
      strsplit(interactions, ":", fixed = TRUE),
      function(parts) "sc_fa" %in% parts,
      logical(1)
    )]
    if (length(age_interaction) > 1L ||
        !identical(length(age_interaction) == 1L,
                   isTRUE(edge$age_moderated[[1]]))) {
      stop("Age-moderation metadata disagrees with edge `", edge$edge_id, "`.")
    }
    age_term <- if (length(age_interaction) == 1L) age_interaction else NA_character_
    age_estimate <- if (is.na(age_term)) NA_real_ else coefficients[[age_term]]
    age_se <- if (is.na(age_term)) NA_real_ else sqrt(covariance[age_term, age_term])

    tibble::tibble(
      edge_id = edge$edge_id[[1]],
      native_n = stats::nobs(model),
      family = stats::family(model)$family,
      link = stats::family(model)$link,
      moderation_terms = paste(interactions, collapse = ";"),
      inference_estimand = path_test$inference_estimand[[1]],
      reference_conditioning = path_test$reference_conditioning[[1]],
      estimate = path_test$inference_estimate[[1]],
      std_error = path_test$inference_standard_error[[1]],
      conf_low = path_test$inference_conf_low[[1]],
      conf_high = path_test$inference_conf_high[[1]],
      test_statistic = path_test$inference_test_statistic[[1]],
      p_value = path_test$inference_p_value[[1]],
      age_interaction_term = age_term,
      age_interaction_estimate = age_estimate,
      age_interaction_std_error = age_se,
      age_interaction_conf_low = age_estimate - stats::qnorm(0.975) * age_se,
      age_interaction_conf_high = age_estimate + stats::qnorm(0.975) * age_se,
      age_interaction_p_value = if (is.na(age_term)) {
        NA_real_
      } else {
        2 * stats::pnorm(-abs(age_estimate / age_se))
      }
    )
  }))

  standardization_columns <- c(
    "edge_id", "random_effects", "averaging", "conditioning",
    "average_partial_effect", "sd_x_reference", "sd_y_reference",
    "standardized_response_effect", "oe_r2", "standardized_oe_effect",
    "rr_x_low_model", "rr_x_high_model", "rr_y_low_model",
    "rr_y_high_model", "rr_response_contrast",
    "standardized_relevant_range_effect", "rr_x_range_supported",
    "standardization_error"
  )

  edge_spec |>
    dplyr::select(dplyr::all_of(required)) |>
    dplyr::left_join(inference, by = "edge_id") |>
    dplyr::left_join(
      standardized |>
        dplyr::select(dplyr::all_of(standardization_columns)),
      by = "edge_id"
    )
}

# Save the finalized parsimonious SEM as one lossless R handoff accompanied by
# flat CSV summaries for inspection outside R. The RDS is authoritative because
# it preserves fitted models, nested graph objects, and list-columns.
# Details: Validates consistency among models, graph, paths, logs, and R2; then
# flattens list-columns for CSV while preserving lossless objects in the RDS.
# Returns: A list containing the in-memory handoff and all written file paths.
export_parsimonious_sem_handoff <- function(
    models,
    graph,
    basis_set,
    missing_path_results,
    global_fit,
    path_table,
    reduction_logs,
    component_r2,
    output_dir,
    file_stem = "parsimonious_sem") {
  .sem_assert_scalar_character(output_dir, "output_dir")
  .sem_assert_scalar_character(file_stem, "file_stem")
  if (is.null(names(models)) || any(!nzchar(names(models))) ||
      anyDuplicated(names(models))) {
    stop("`models` must be a uniquely named list.")
  }
  required_graph_objects <- c(
    "directed_edges", "nondirectional_associations",
    "removed_edges", "model_terms"
  )
  if (!is.list(graph) ||
      !all(required_graph_objects %in% names(graph))) {
    stop("`graph` is not a complete assembled SEM graph.")
  }
  if (!is.list(reduction_logs) || is.null(names(reduction_logs)) ||
      any(!nzchar(names(reduction_logs))) || anyDuplicated(names(reduction_logs))) {
    stop("`reduction_logs` must be a uniquely named list.")
  }
  if (!is.data.frame(path_table) ||
      !setequal(path_table$edge_id, graph$directed_edges$edge_id) ||
      anyDuplicated(path_table$edge_id)) {
    stop("`path_table` must contain each retained directed edge exactly once.")
  }
  if (!is.data.frame(component_r2) ||
      !setequal(component_r2$model_id, names(models))) {
    stop("Component R2 must cover every finalized model.")
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  handoff <- list(
    sem_models = models,
    graph = graph,
    basis_set = basis_set,
    missing_path_results = missing_path_results,
    global_fit = global_fit,
    path_table = path_table,
    reduction_logs = reduction_logs,
    component_r2 = component_r2
  )

  collapse_character_list <- function(values) {
    vapply(
      values,
      function(value) paste(unlist(value, use.names = FALSE), collapse = ";"),
      character(1)
    )
  }
  format_named_list <- function(values) {
    vapply(values, function(value) {
      if (length(value) == 0L) return("")
      value_names <- names(value)
      if (is.null(value_names)) {
        paste(unlist(value, use.names = FALSE), collapse = ";")
      } else {
        paste(
          paste0(value_names, "=", unlist(value, use.names = FALSE)),
          collapse = ";"
        )
      }
    }, character(1))
  }

  directed_edges_export <- graph$directed_edges
  if ("at" %in% names(directed_edges_export)) {
    directed_edges_export$at <- format_named_list(directed_edges_export$at)
  }
  basis_set_export <- basis_set
  if ("conditioning_nodes" %in% names(basis_set_export)) {
    basis_set_export$conditioning_nodes <- collapse_character_list(
      basis_set_export$conditioning_nodes
    )
  }
  if ("conditioning_terms" %in% names(basis_set_export)) {
    basis_set_export$conditioning_terms <- collapse_character_list(
      basis_set_export$conditioning_terms
    )
  }
  missing_paths_export <- missing_path_results
  if ("conditioning_nodes" %in% names(missing_paths_export)) {
    missing_paths_export$conditioning_nodes <- collapse_character_list(
      missing_paths_export$conditioning_nodes
    )
  }
  reduction_log_export <- dplyr::bind_rows(
    reduction_logs,
    .id = "component_registry"
  )

  output_files <- c(
    handoff = file.path(output_dir, paste0(file_stem, "_handoff.rds")),
    directed_edges = file.path(
      output_dir, paste0(file_stem, "_directed_edges.csv")
    ),
    nondirectional_associations = file.path(
      output_dir, paste0(file_stem, "_nondirectional_associations.csv")
    ),
    basis_set = file.path(output_dir, paste0(file_stem, "_basis_set.csv")),
    missing_paths = file.path(
      output_dir, paste0(file_stem, "_missing_paths.csv")
    ),
    global_fit = file.path(
      output_dir, paste0(file_stem, "_global_fit.csv")
    ),
    path_table = file.path(
      output_dir, paste0(file_stem, "_paths.csv")
    ),
    reduction_log = file.path(
      output_dir, paste0(file_stem, "_reduction_log.csv")
    ),
    component_r2 = file.path(
      output_dir, paste0(file_stem, "_component_r2.csv")
    )
  )

  utils::write.csv(
    directed_edges_export,
    output_files[["directed_edges"]],
    row.names = FALSE
  )
  utils::write.csv(
    graph$nondirectional_associations,
    output_files[["nondirectional_associations"]],
    row.names = FALSE
  )
  utils::write.csv(
    basis_set_export,
    output_files[["basis_set"]],
    row.names = FALSE
  )
  utils::write.csv(
    missing_paths_export,
    output_files[["missing_paths"]],
    row.names = FALSE
  )
  utils::write.csv(
    global_fit,
    output_files[["global_fit"]],
    row.names = FALSE
  )
  utils::write.csv(
    path_table,
    output_files[["path_table"]],
    row.names = FALSE
  )
  utils::write.csv(
    reduction_log_export,
    output_files[["reduction_log"]],
    row.names = FALSE
  )
  utils::write.csv(
    component_r2,
    output_files[["component_r2"]],
    row.names = FALSE
  )
  saveRDS(handoff, output_files[["handoff"]])

  list(
    handoff = handoff,
    files = output_files
  )
}

# -----------------------------------------------------------------------------
# Model-refit validation
# -----------------------------------------------------------------------------

# Confirm that a refitted component has finite coefficients and a valid fit.
# Details: For glmmTMB fits, requires optimizer convergence, a positive-definite
# Hessian, and finite fixed effects; other classes must have finite coefficients.
# Returns: A list with an `ok` flag and a concise failure reason.
validate_sem_refit <- function(model) {
  coefficients <- tryCatch(
    if (inherits(model, "glmmTMB")) {
      unlist(glmmTMB::fixef(model), use.names = TRUE)
    } else {
      stats::coef(model)
    },
    error = function(error) numeric(0)
  )
  finite_coefficients <- length(coefficients) > 0L &&
    all(is.finite(coefficients))

  if (inherits(model, "glmmTMB")) {
    convergence <- isTRUE(model$fit$convergence == 0L)
    positive_hessian <- isTRUE(model$sdr$pdHess)
    ok <- convergence && positive_hessian && finite_coefficients
    reason <- paste(
      c(
        if (!convergence) "optimizer convergence code was nonzero",
        if (!positive_hessian) "Hessian was not positive definite",
        if (!finite_coefficients) "fixed coefficients were non-finite"
      ),
      collapse = "; "
    )
  } else {
    ok <- finite_coefficients
    reason <- if (ok) "" else "coefficients were absent or non-finite"
  }
  list(ok = ok, reason = reason)
}
