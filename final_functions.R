###-------------------------------------------------------------------------------------
# Helper functions for data analysis of manuscript "Forest growth strengthens diversity effects on multitrophic interactions"
# Author: Massimo Martini

###-------------------------------------------------------------------------------------

#calculating forest age effect from different models
calculate_forest_age_effect <- function(model, variable, forest_age_var = "sc_fa") {
  
  # Get data from model
  data <- model$frame
  
  # Get forest age range
  fa_range <- range(data[[forest_age_var]], na.rm = TRUE)
  
  # Create prediction data (mean values for numeric variables only)
  pred_data <- data.frame(lapply(data, function(x) {
    if (is.numeric(x)) {
      rep(mean(x, na.rm = TRUE), 4)
    } else {
      rep(x[1], 4)  # Use first value for non-numeric
    }
  }))
  
  # Set forest age to min/max and variable to +/- 0.5 SD
  var_mean <- mean(data[[variable]], na.rm = TRUE)
  var_sd <- sd(data[[variable]], na.rm = TRUE)
  pred_data[[forest_age_var]] <- rep(fa_range, each = 2)
  pred_data[[variable]] <- rep(c(var_mean - 0.5 * var_sd, var_mean + 0.5 * var_sd), 2)
  
  # Get predictions (suppress random effects warning)
  pred <- predict(model, newdata = pred_data, type = "link", allow.new.levels = TRUE)
  
  # Calculate effects
  effect_min <- pred[2] - pred[1]  # at min forest age
  effect_max <- pred[4] - pred[3]  # at max forest age
  
  # Results with verification info
  list(
    variable = variable,
    forest_age_range = fa_range,
    effect_at_min_age = effect_min,
    effect_at_max_age = effect_max,
    effect_difference = effect_max - effect_min,
    
    # Verification info
    prediction_setup = data.frame(
      scenario = c("min_fa_low_var", "min_fa_high_var", "max_fa_low_var", "max_fa_high_var"),
      forest_age = pred_data[[forest_age_var]],
      variable_value = pred_data[[variable]],
      prediction = pred
    ),
    variable_contrast = paste("Mean ± 0.5 SD:", 
                              round(var_mean - 0.5 * var_sd, 3), "to", 
                              round(var_mean + 0.5 * var_sd, 3))
  )
}


###-------------------------------------------------------------------------------------
#Sensitivity analysis
compare_models <- function(model1, model2) {
  # Extract summary coefficient tables
  coef1 <- summary(model1)$coefficients$cond
  coef2 <- summary(model2)$coefficients$cond
  
  # Align parameter names (assumes models have the same set/order of fixed effects)
  params <- intersect(rownames(coef1), rownames(coef2)) # to cover cases where parameters differ
  coef1 <- coef1[params, , drop = FALSE]
  coef2 <- coef2[params, , drop = FALSE]
  
  # Build the comparison table
  comp <- data.frame(
    Parameter     = params,
    Estimate_1    = coef1[,"Estimate"],
    SE_1          = coef1[,"Std. Error"],
    P_1           = coef1[,"Pr(>|z|)"],
    Estimate_2    = coef2[,"Estimate"],
    SE_2          = coef2[,"Std. Error"],
    P_2           = coef2[,"Pr(>|z|)"]
  )
  
  # Add significance and direction change columns
  comp$Signif_1   <- comp$P_1 < 0.05
  comp$Signif_2   <- comp$P_2 < 0.05
  comp$SignifChange <- comp$Signif_1 != comp$Signif_2
  comp$DirChange    <- sign(comp$Estimate_1) != sign(comp$Estimate_2)
  
  return(comp)
}


#---------------------------------------------------------------------------------------------
#Remove path analysis sub-models from the global environment once they are already safely inside the various lists 
rm_path_objects <- function(...) {
  objs <- list(...)
  all_names <- unlist(lapply(objs, function(x) {
    nms <- names(x)
    # Keep non-empty names, skip if NA or blank
    nms[!is.na(nms) & trimws(nms) != ""]
  }))
  # Only keep names that are valid identifiers and exist in the global workspace
  valid_names <- all_names[make.names(all_names) == all_names & all_names %in% ls(envir = .GlobalEnv)]
  if (length(valid_names) > 0) rm(list = valid_names, envir = .GlobalEnv)
}


#---------------------------------------------------------------------------------------------
#calculating model dispersion parameters and creating a table for printing
get_disp_row <- function(mod, name, nsim = 500) {
  sim <- simulateResiduals(mod, n = nsim, plot = FALSE, refit = FALSE)
  td  <- suppressWarnings(testDispersion(sim, plot = FALSE))
  tibble(
    Model      = name,
    Family     = tryCatch(family(mod)$family, error = function(e) NA_character_),
    Dispersion = if ("statistic" %in% names(td)) as.numeric(td$statistic) else NA_real_,
    P_value    = if ("p.value"   %in% names(td)) as.numeric(td$p.value)    else NA_real_
  )
}


#---------------------------------------------------------------------------------------------
#setting a pretty names map
pretty_map <- c(
  "Tree_rich."   = "sc_logtr",
  "Stand_vol."  = "sc_sv",
  "Tree_FD"       = "sc_fd",
  "Stand_age"    = "sc_fa",
  "Host_abund."    = "sc_ha",
  "Host_rich."     = "sc_hr",
  "Par_abund."   = "sc_pa",
  "Par_rich."    = "sc_pr",
  "Par_rich."    = "sc_pr10",
  "Slope"         = "sc_slope",
  "Elevation"     = "sc_elev",
  "Eastness"      = "sc_east",
  "Northness"     = "sc_north",
  "Host_abund." = "sc_cells",
  "Focal_rich"  = "sc_frich"
)
#printing an anova type I table
make_type1_table <- function(mods, anova_obj, pretty_map,
                             digits = 3,
                             baseline_label = "(baseline: intercept + RE)") {
  stopifnot(is.list(mods), length(mods) >= 1)
  
  # ---- helpers --------------------------------------------------------------
  get_terms <- function(m) attr(terms(stats::formula(m)), "term.labels")
  
  # reverse lookup code -> pretty
  code_to_pretty <- function(code) {
    hits <- names(pretty_map)[match(code, unname(pretty_map))]
    ifelse(is.na(hits), code, hits)
  }
  pretty_one_term <- function(term) {
    parts <- strsplit(term, ":", fixed = TRUE)[[1]]
    # use ASCII 'x' to avoid encoding issues in CSV/Excel
    paste(vapply(parts, code_to_pretty, character(1)), collapse = " x ")
  }
  
  # Normalize/standardize ANOVA column names to a canonical set
  std_names <- function(df) {
    nn_raw <- names(df)
    nn_std <- gsub("[^A-Za-z]", "", nn_raw)  # drop spaces, symbols
    out_names <- character(length(nn_std))
    for (i in seq_along(nn_std)) {
      s <- nn_std[i]
      if (s == "Chisq") out_names[i] <- "Chisq"
      else if (s %in% c("PrChisq","Prchisq","PrChisQ","PrGTChisq")) out_names[i] <- "Pr(>Chisq)"
      else if (s %in% c("ChiDf","ChiDF","ChDf")) out_names[i] <- "Chi Df"
      else if (s %in% c("AIC")) out_names[i] <- "AIC"
      else if (s %in% c("Df","df")) out_names[i] <- "Df"
      else if (s %in% c("logLik","loglik")) out_names[i] <- "logLik"
      else if (s %in% c("deviance","Deviance")) out_names[i] <- "deviance"
      else if (s %in% c("BIC")) out_names[i] <- "BIC"
      else out_names[i] <- nn_raw[i]  # keep original if unknown
    }
    names(df) <- out_names
    df
  }
  
  # ---- align objects & names ------------------------------------------------
  rn <- names(mods)
  if (is.null(rn)) rn <- paste0("m", seq_along(mods) - 1L)
  
  a <- as.data.frame(anova_obj, stringsAsFactors = FALSE)
  a <- std_names(a)
  stopifnot(nrow(a) == length(rn))
  rownames(a) <- rn
  
  # ---- compute Added terms --------------------------------------------------
  term_list <- lapply(mods, get_terms)
  added <- character(length(mods))
  
  # --- baseline row: show random factor structure ---
  re_terms <- lme4::findbars(formula(mods[[1]]))
  if (length(re_terms) > 0) {
    re_names <- vapply(re_terms, function(x) deparse(x[[3]]), character(1))
    re_names <- gsub(":", "/", re_names)
    re_names <- gsub("_id", "", re_names, fixed = TRUE)
    re_names <- unique(re_names)
    for (r in re_names) {
      if (grepl("/", r)) {
        parent <- sub("/.*", "", r)
        re_names <- setdiff(re_names, parent)
      }
    }
    re_names <- c(re_names[grepl("/", re_names)], re_names[!grepl("/", re_names)])
    added[1] <- paste0("Random factors: ", paste(re_names, collapse = " + "))
  } else {
    added[1] <- "(no random factors)"
  }
  
  # --- detect terms added in later models ---
  for (i in 2:length(mods)) {
    new_terms <- setdiff(term_list[[i]], term_list[[i - 1]])
    added[i] <- if (length(new_terms) == 0) "—" else
      paste(vapply(new_terms, pretty_one_term, character(1)), collapse = " + ")
  }
  
  # ---- choose/ensure column order (keep only the essentials) ----------------
  desired <- c("AIC", "Chi Df", "Chisq", "Pr(>Chisq)")
  present <- intersect(desired, names(a))
  extra <- setdiff(names(a), present)  # we won't include extras
  a2 <- a[, present, drop = FALSE]
  
  # --- add significance stars column right after p-value --------------------
  if ("Pr(>Chisq)" %in% names(a2)) {
    sig <- cut(a2[["Pr(>Chisq)"]],
               breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
               labels = c("***", "**", "*", ".", ""))
    idx <- which(names(a2) == "Pr(>Chisq)")
    a2 <- cbind(a2[seq_len(idx)], Signif = sig, a2[-seq_len(idx), drop = FALSE])
  }
  
  # ---- build output ---------------------------------------------------------
  out <- cbind(
    Model = rn,
    Added = added,
    a2,
    deparse.level = 0
  )
  
  # numeric formatting
  num_cols <- setdiff(names(out), c("Model","Added","Signif"))
  
  for (cl in num_cols) {
    if (cl == "Pr(>Chisq)" && is.numeric(out[[cl]])) {
      pvals <- out[[cl]]
      out[[cl]] <- ifelse(
        pvals < 0.001,
        "<0.001",
        sprintf(paste0("%.", digits, "f"), pvals)
      )
    } else if (is.numeric(out[[cl]])) {
      out[[cl]] <- round(out[[cl]], digits)
    }
  }
  
  out
}


# ------------------------------------------------------------------------------
# Clean model-summary exporter for glmmTMB models
#
# Usage:
#   source("functions.R")
#   model_summaries <- export_glmmTMB_markdown(
#     models = models,
#     file = "model_summaries.md",
#     digits = 3
#   )
# ------------------------------------------------------------------------------

sig_stars <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*",
                              ifelse(p < 0.1, ".", "")))))
}

fmt_value <- function(x, digits = 3) {
  if (length(x) == 0 || is.null(x)) return(NA_character_)
  if (is.character(x)) return(x)
  ifelse(is.na(x), "NA", formatC(as.numeric(x), format = "f", digits = digits))
}

fmt_p <- function(x, digits = 3) {
  ifelse(is.na(x), "NA",
         ifelse(x < 0.001, "<0.001", formatC(x, format = "f", digits = digits)))
}

fmt_sig <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", format(signif(as.numeric(x), digits), scientific = FALSE, trim = TRUE))
}

fmt_docx_numeric_cells <- function(x, digits = 3, p_value = FALSE) {
  x_chr <- as.character(x)
  x_num <- suppressWarnings(as.numeric(x_chr))
  is_num <- !is.na(x_num) & nzchar(trimws(x_chr))

  out <- x_chr
  out[is.na(out)] <- ""
  if (any(is_num)) {
    out[is_num] <- if (p_value) fmt_p(x_num[is_num], digits) else fmt_sig(x_num[is_num], digits)
  }
  out
}

capture_warnings <- function(expr,
                             object_name = "script_warnings",
                             envir = .GlobalEnv,
                             append = TRUE,
                             quiet = TRUE) {
  captured <- list()

  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      captured[[length(captured) + 1L]] <<- list(
        message = conditionMessage(w),
        call = paste(deparse(conditionCall(w)), collapse = " ")
      )

      if (quiet) invokeRestart("muffleWarning")
    }
  )

  warning_table <- data.frame(
    Index = seq_along(captured),
    Message = vapply(captured, `[[`, character(1), "message"),
    Call = vapply(captured, `[[`, character(1), "call"),
    stringsAsFactors = FALSE
  )

  if (append && exists(object_name, envir = envir, inherits = FALSE)) {
    old <- get(object_name, envir = envir, inherits = FALSE)
    if (is.data.frame(old) && nrow(old) > 0) {
      warning_table <- rbind(old, warning_table)
      warning_table$Index <- seq_len(nrow(warning_table))
    }
  }

  assign(object_name, warning_table, envir = envir)
  invisible(value)
}

pad_markdown_table <- function(df, digits = 3, p_cols = character(0)) {
  if (is.null(df) || nrow(df) == 0) return(character(0))
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)

  for (nm in names(df)) {
    if (is.numeric(df[[nm]])) {
      df[[nm]] <- if (nm %in% p_cols) fmt_p(df[[nm]], digits) else fmt_value(df[[nm]], digits)
    } else {
      df[[nm]] <- ifelse(is.na(df[[nm]]), "NA", as.character(df[[nm]]))
    }
  }

  widths <- vapply(names(df), function(nm) {
    max(nchar(nm, type = "width"), nchar(df[[nm]], type = "width"), na.rm = TRUE)
  }, numeric(1))

  pad_row <- function(vals) {
    vals <- as.character(vals)
    paste(mapply(function(val, width) sprintf(paste0("%-", width, "s"), val),
                 vals, widths),
          collapse = " | ")
  }

  c(
    "```",
    pad_row(names(df)),
    paste(mapply(function(width) paste(rep("-", width), collapse = ""),
                 widths),
          collapse = "-|-"),
    apply(df, 1, pad_row),
    "```"
  )
}

default_pretty_names <- function() {
  c(
    sc_logtr = "Tree richness",
    sc_sv = "Stand volume",
    sc_fd = "Tree FD",
    sc_fa = "Stand age",
    sc_elev = "Elevation",
    sc_east = "Eastness",
    sc_north = "Northness",
    sc_slope = "Slope",
    sc_temp = "Annual temperature",
    sc_humid = "Annual humidity",
    sc_cells = "Host abundance",
    sc_hr = "Host richness",
    sc_pr = "Parasitoid richness",
    sc_pr10 = "Parasitoid richness",
    sc_netsize = "Network size",
    sc_links = "Number of links",
    sc_linkdense = "Linkage density",
    sc_h2 = "H2",
    sc_mdprime = "Mean d-prime",
    sc_niche = "Niche overlap",
    sc_intev = "Interaction evenness",
    sc_robust = "Robustness",
    sc_frich = "Focal parasitoid richness"
  )
}

pretty_term <- function(term, pretty_names = default_pretty_names()) {
  if (term %in% c("(Intercept)", "(Intercept).1")) return(term)

  clean_piece <- function(x) {
    if (x %in% names(pretty_names)) return(pretty_names[[x]])

    y <- x
    for (code in names(pretty_names)) {
      y <- gsub(paste0("\\b", code, "\\b"), pretty_names[[code]], y)
    }
    y
  }

  pieces <- strsplit(term, ":", fixed = TRUE)[[1]]
  pieces <- vapply(pieces, clean_piece, character(1))
  paste(pieces, collapse = " x ")
}

safe_family <- function(model) {
  fam <- try(stats::family(model), silent = TRUE)
  if (inherits(fam, "try-error")) {
    return(list(family = NA_character_, link = NA_character_))
  }
  list(
    family = if (!is.null(fam$family)) fam$family else NA_character_,
    link = if (!is.null(fam$link)) fam$link else NA_character_
  )
}

extract_fixed_table <- function(model, component = c("cond", "zi", "disp"),
                                pretty_names = default_pretty_names()) {
  component <- match.arg(component)
  coefs <- try(coef(summary(model)), silent = TRUE)
  if (inherits(coefs, "try-error") || is.null(coefs[[component]]) || nrow(coefs[[component]]) == 0) {
    return(NULL)
  }

  tab <- as.data.frame(coefs[[component]], check.names = FALSE)
  stat_col <- intersect(c("z value", "t value"), names(tab))
  p_col <- intersect(c("Pr(>|z|)", "Pr(>|t|)"), names(tab))

  out <- data.frame(
    Term = vapply(rownames(tab), pretty_term, character(1), pretty_names = pretty_names),
    Estimate = tab[["Estimate"]],
    `Std. Error` = tab[["Std. Error"]],
    check.names = FALSE
  )

  if (length(stat_col) > 0) out[[stat_col[1]]] <- tab[[stat_col[1]]]
  if (length(p_col) > 0) {
    out[["P value"]] <- tab[[p_col[1]]]
    out[["Sig."]] <- sig_stars(tab[[p_col[1]]])
  }

  out
}

extract_random_effects <- function(model) {
  if (requireNamespace("broom.mixed", quietly = TRUE)) {
    ran <- try(broom.mixed::tidy(model, effects = "ran_pars"), silent = TRUE)
    if (!inherits(ran, "try-error") && !is.null(ran) && nrow(ran) > 0) {
      ran <- subset(ran, component %in% c("cond", NA) & grepl("^sd__", term))
      if (nrow(ran) > 0) {
        return(data.frame(
          Component = "conditional",
          Group = ran$group,
          Term = sub("^sd__", "", ran$term),
          `Std. Dev.` = ran$estimate,
          Variance = ran$estimate^2,
          check.names = FALSE
        ))
      }
    }
  }

  vc <- try(as.data.frame(VarCorr(model)), silent = TRUE)
  if (!inherits(vc, "try-error") && !is.null(vc) && nrow(vc) > 0) {
    if (all(c("component", "grp", "var1", "var2", "sdcor") %in% names(vc))) {
      vc <- subset(vc, component == "cond" & is.na(var2))
      if (nrow(vc) > 0) {
        return(data.frame(
          Component = "conditional",
          Group = vc$grp,
          Term = ifelse(is.na(vc$var1), "(Intercept)", vc$var1),
          `Std. Dev.` = vc$sdcor,
          Variance = vc$sdcor^2,
          check.names = FALSE
        ))
      }
    }
  }

  NULL
}

extract_dispersion_random_from_summary <- function(model) {
  txt <- capture.output(suppressWarnings(summary(model)))
  start <- grep("^Dispersion model:", txt)
  if (length(start) == 0) return(NULL)

  start <- start[1] + 1L
  if (start > length(txt)) return(NULL)

  stop_at <- length(txt)
  for (i in seq(from = start, to = length(txt))) {
    line <- trimws(txt[[i]])
    if (i > start && (
      line == "" ||
      grepl("^Number of obs:", line) ||
      grepl("^Conditional model:", line) ||
      grepl("^Zero-inflation model:", line)
    )) {
      stop_at <- i - 1L
      break
    }
  }

  block <- txt[start:stop_at]
  block <- block[nzchar(trimws(block))]
  if (length(block) < 2) return(NULL)

  header_idx <- grep("Groups\\s+Name\\s+Variance\\s+Std\\.Dev\\.", block)
  if (length(header_idx) == 0) return(NULL)

  rows <- block[(header_idx[1] + 1L):length(block)]
  rows <- rows[nzchar(trimws(rows))]
  if (length(rows) == 0) return(NULL)

  parsed <- lapply(rows, function(line) {
    fields <- strsplit(trimws(line), "\\s+")[[1]]
    if (length(fields) < 4) return(NULL)

    numeric_at <- which(!is.na(suppressWarnings(as.numeric(fields))))
    if (length(numeric_at) < 2) return(NULL)

    var_pos <- numeric_at[1]
    sd_pos <- numeric_at[2]
    if (var_pos < 3) return(NULL)

    group <- fields[1]
    term <- paste(fields[2:(var_pos - 1L)], collapse = " ")
    variance <- suppressWarnings(as.numeric(fields[var_pos]))
    std_dev <- suppressWarnings(as.numeric(fields[sd_pos]))
    corr <- if (sd_pos < length(fields)) paste(fields[(sd_pos + 1L):length(fields)], collapse = " ") else NA_character_

    data.frame(
      Group = group,
      Term = term,
      Variance = variance,
      `Std. Dev.` = std_dev,
      Correlation = corr,
      check.names = FALSE
    )
  })

  parsed <- parsed[!vapply(parsed, is.null, logical(1))]
  if (length(parsed) == 0) return(NULL)

  out <- do.call(rbind, parsed)
  rownames(out) <- NULL
  out
}

extract_dispersion_parameters <- function(model, pretty_names = default_pretty_names()) {
  out <- list()

  disp_fixed <- extract_fixed_table(model, "disp", pretty_names)
  if (!is.null(disp_fixed)) out$dispersion_fixed <- disp_fixed

  vc <- try(as.data.frame(VarCorr(model)), silent = TRUE)
  if (!inherits(vc, "try-error") && !is.null(vc) && nrow(vc) > 0 &&
      "component" %in% names(vc)) {
    disp_vc <- subset(vc, component == "disp")
    if (nrow(disp_vc) > 0) {
      keep <- intersect(c("grp", "var1", "var2", "vcov", "sdcor"), names(disp_vc))
      disp_vc <- disp_vc[, keep, drop = FALSE]
      names(disp_vc) <- c("Group", "Term 1", "Term 2", "Variance/Covariance", "Std.Dev./Corr")[seq_along(keep)]
      out$dispersion_random <- disp_vc
    }
  }

  if (is.null(out$dispersion_random)) {
    disp_summary <- extract_dispersion_random_from_summary(model)
    if (!is.null(disp_summary)) out$dispersion_random <- disp_summary
  }

  family_name <- safe_family(model)$family
  if (!is.na(family_name) && grepl("^genpois", family_name, ignore.case = TRUE)) {
    sigma_val <- try(sigma(model), silent = TRUE)
    if (!inherits(sigma_val, "try-error") && is.finite(sigma_val)) {
      out$family_dispersion <- data.frame(
        Parameter = "Generalized Poisson dispersion",
        Estimate = as.numeric(sigma_val),
        check.names = FALSE
      )
    }
  }

  if (length(out) == 0) NULL else out
}

fit_null_model <- function(model) {
  tryCatch(
    suppressWarnings(update(model, . ~ 1)),
    error = function(e) structure(list(error = conditionMessage(e)), class = "summary_export_error")
  )
}

likelihood_r2 <- function(model, null_model = NULL) {
  if (is.null(null_model)) null_model <- fit_null_model(model)
  if (inherits(null_model, "summary_export_error")) {
    return(list(
      values = c(McFadden = NA_real_, CoxSnell = NA_real_, Nagelkerke = NA_real_),
      note = paste("Likelihood R2 failed:", null_model$error)
    ))
  }

  ll_full <- try(as.numeric(stats::logLik(model)), silent = TRUE)
  ll_null <- try(as.numeric(stats::logLik(null_model)), silent = TRUE)
  n <- try(stats::nobs(model), silent = TRUE)

  if (inherits(ll_full, "try-error") || inherits(ll_null, "try-error") ||
      inherits(n, "try-error") || !is.finite(ll_full) || !is.finite(ll_null) ||
      !is.finite(n) || n <= 0) {
    return(list(
      values = c(McFadden = NA_real_, CoxSnell = NA_real_, Nagelkerke = NA_real_),
      note = "Likelihood R2 failed because log-likelihood or n was not finite."
    ))
  }

  cox_snell <- 1 - exp((2 / n) * (ll_null - ll_full))
  max_cox_snell <- 1 - exp((2 / n) * ll_null)

  list(
    values = c(
      McFadden = 1 - (ll_full / ll_null),
      CoxSnell = cox_snell,
      Nagelkerke = if (abs(max_cox_snell) < .Machine$double.eps) NA_real_ else cox_snell / max_cox_snell
    ),
    note = NA_character_
  )
}

nakagawa_r2 <- function(model, tolerances = c(1e-08, 1e-10, 1e-12, 0)) {
  if (!requireNamespace("performance", quietly = TRUE)) {
    return(list(
      values = c(Marginal = NA_real_, Conditional = NA_real_),
      note = "Nakagawa R2 not calculated because package 'performance' is unavailable."
    ))
  }

  messages <- character(0)
  for (tol in tolerances) {
    res <- tryCatch(
      suppressWarnings(performance::r2_nakagawa(model, tolerance = tol)),
      error = function(e) {
        messages <<- c(messages, conditionMessage(e))
        NULL
      }
    )

    vals <- suppressWarnings(try(unlist(res), silent = TRUE))
    if (!inherits(vals, "try-error") && length(vals) > 0) {
      m_idx <- grep("marginal", names(vals), ignore.case = TRUE)[1]
      c_idx <- grep("conditional", names(vals), ignore.case = TRUE)[1]
      m <- suppressWarnings(as.numeric(vals[m_idx]))
      c <- suppressWarnings(as.numeric(vals[c_idx]))
      if (is.finite(m) || is.finite(c)) {
        note <- if (identical(tol, tolerances[1])) {
          NA_character_
        } else {
          paste0("Nakagawa R2 required a less sensitive singularity tolerance (",
                 format(tol, scientific = TRUE), ").")
        }
        return(list(values = c(Marginal = m, Conditional = c), note = note))
      }
    }
  }

  if (requireNamespace("MuMIn", quietly = TRUE)) {
    mm <- try(suppressWarnings(MuMIn::r.squaredGLMM(model)), silent = TRUE)
    if (!inherits(mm, "try-error") && !is.null(mm) &&
        all(c("R2m", "R2c") %in% colnames(mm))) {
      return(list(
        values = c(Marginal = as.numeric(mm[1, "R2m"]), Conditional = as.numeric(mm[1, "R2c"])),
        note = "Nakagawa R2 was calculated with MuMIn::r.squaredGLMM() after performance::r2_nakagawa() failed."
      ))
    }
  }

  list(
    values = c(Marginal = NA_real_, Conditional = NA_real_),
    note = paste(c("Nakagawa R2 failed for this model.", unique(messages)), collapse = " ")
  )
}

model_fit_table <- function(model, lr2, nr2) {
  fam <- safe_family(model)

  data.frame(
    Family = fam$family,
    Link = fam$link,
    Observations = tryCatch(stats::nobs(model), error = function(e) NA_integer_),
    logLik = tryCatch(as.numeric(stats::logLik(model)), error = function(e) NA_real_),
    AIC = tryCatch(stats::AIC(model), error = function(e) NA_real_),
    R2_Nakagawa_marginal = nr2$values[["Marginal"]],
    R2_Nakagawa_conditional = nr2$values[["Conditional"]],
    R2_McFadden = lr2$values[["McFadden"]],
    R2_CoxSnell = lr2$values[["CoxSnell"]],
    R2_Nagelkerke = lr2$values[["Nagelkerke"]],
    check.names = FALSE
  )
}

export_glmmTMB_markdown <- function(models,
                                    file = "model_summaries.md",
                                    digits = 3,
                                    pretty_names = default_pretty_names(),
                                    null_models = NULL,
                                    title = "Model summaries") {
  stopifnot(is.list(models))
  if (length(models) == 0) stop("`models` is empty.")

  model_names <- names(models)
  if (is.null(model_names)) model_names <- paste0("Model_", seq_along(models))
  model_names[!nzchar(model_names)] <- paste0("Model_", which(!nzchar(model_names)))

  if (!is.null(null_models) && is.null(names(null_models))) {
    stop("If `null_models` is supplied, it must be a named list.")
  }

  con <- file(file, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)

  writeLines(c(
    paste0("# ", title),
    "",
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "R2 reporting:",
    "",
    "- Nakagawa R2 is reported as marginal / conditional when available.",
    "- Likelihood-based pseudo-R2 is reported as McFadden / Cox-Snell / Nagelkerke.",
    "- If Nakagawa R2 required a less sensitive singularity tolerance, the model section reports that explicitly.",
    "- R2 values from different families or dispersion structures should be compared cautiously.",
    "",
    "---",
    ""
  ), con)

  results <- vector("list", length(models))
  names(results) <- model_names

  for (i in seq_along(models)) {
    model <- models[[i]]
    model_name <- model_names[[i]]
    null_model <- if (!is.null(null_models) && model_name %in% names(null_models)) null_models[[model_name]] else NULL

    lr2 <- likelihood_r2(model, null_model)
    nr2 <- nakagawa_r2(model)
    fit <- model_fit_table(model, lr2, nr2)

    writeLines(c(paste0("## ", model_name), ""), con)
    writeLines(pad_markdown_table(fit, digits), con)
    writeLines("", con)

    notes <- c(lr2$note, nr2$note)
    notes <- notes[!is.na(notes) & nzchar(notes)]
    if (length(notes) > 0) {
      writeLines(c("**R2 notes**", "", paste0("- ", notes), ""), con)
    }

    cond <- extract_fixed_table(model, "cond", pretty_names)
    if (!is.null(cond)) {
      writeLines(c("### Fixed effects: conditional", ""), con)
      writeLines(pad_markdown_table(cond, digits, p_cols = "P value"), con)
      writeLines(c("", "Significance codes: *** < 0.001, ** < 0.01, * < 0.05, . < 0.1", ""), con)
    }

    zi <- extract_fixed_table(model, "zi", pretty_names)
    if (!is.null(zi)) {
      writeLines(c("### Fixed effects: zero-inflation", ""), con)
      writeLines(pad_markdown_table(zi, digits, p_cols = "P value"), con)
      writeLines(c("", "Significance codes: *** < 0.001, ** < 0.01, * < 0.05, . < 0.1", ""), con)
    }

    disp <- extract_dispersion_parameters(model, pretty_names)
    if (!is.null(disp)) {
      if (!is.null(disp$family_dispersion)) {
        writeLines(c("### Family dispersion", ""), con)
        writeLines(pad_markdown_table(disp$family_dispersion, digits), con)
        writeLines("", con)
      }
      if (!is.null(disp$dispersion_fixed)) {
        writeLines(c("### Fixed effects: dispersion model", ""), con)
        writeLines(pad_markdown_table(disp$dispersion_fixed, digits, p_cols = "P value"), con)
        writeLines("", con)
      }
      if (!is.null(disp$dispersion_random)) {
        writeLines(c("### Random effects: dispersion/autocorrelation", ""), con)
        writeLines(pad_markdown_table(disp$dispersion_random, digits), con)
        writeLines("", con)
      }
    }

    ran <- extract_random_effects(model)
    if (!is.null(ran)) {
      writeLines(c("### Random effects: conditional", ""), con)
      writeLines(pad_markdown_table(ran, digits), con)
      writeLines("", con)
    }

    writeLines(c("---", ""), con)

    results[[i]] <- list(
      fit = fit,
      r2_notes = notes,
      fixed_conditional = cond,
      fixed_zero_inflation = zi,
      dispersion = disp,
      random_conditional = ran
    )
  }

  message("Wrote model summaries to: ", normalizePath(file, winslash = "/", mustWork = FALSE))
  invisible(results)
}

table_for_docx <- function(model, lr2, nr2, pretty_names = default_pretty_names(), digits = 3) {
  fit <- model_fit_table(model, lr2, nr2)

  fit_rows <- data.frame(
    Section = "Model fit",
    Term = names(fit),
    Estimate = as.character(unlist(fit[1, ], use.names = FALSE)),
    `Std. Error` = "",
    Statistic = "",
    `P value` = "",
    Sig. = "",
    check.names = FALSE
  )

  cond <- extract_fixed_table(model, "cond", pretty_names)
  if (!is.null(cond)) {
    names(cond)[names(cond) %in% c("z value", "t value")] <- "Statistic"
    cond$Section <- "Fixed effects"
  }

  zi <- extract_fixed_table(model, "zi", pretty_names)
  if (!is.null(zi)) {
    names(zi)[names(zi) %in% c("z value", "t value")] <- "Statistic"
    zi$Section <- "Zero-inflation"
  }

  disp <- extract_dispersion_parameters(model, pretty_names)
  disp_rows <- NULL
  if (!is.null(disp)) {
    if (!is.null(disp$family_dispersion)) {
      disp_rows <- rbind(
        disp_rows,
        data.frame(
          Section = "Family dispersion",
          Term = disp$family_dispersion$Parameter,
          Estimate = disp$family_dispersion$Estimate,
          `Std. Error` = "",
          Statistic = "",
          `P value` = "",
          Sig. = "",
          check.names = FALSE
        )
      )
    }

    if (!is.null(disp$dispersion_fixed)) {
      tmp <- disp$dispersion_fixed
      names(tmp)[names(tmp) %in% c("z value", "t value")] <- "Statistic"
      disp_rows <- rbind(
        disp_rows,
        data.frame(
          Section = "Dispersion fixed effects",
          Term = tmp$Term,
          Estimate = tmp$Estimate,
          `Std. Error` = tmp$`Std. Error`,
          Statistic = if ("Statistic" %in% names(tmp)) tmp$Statistic else "",
          `P value` = if ("P value" %in% names(tmp)) tmp$`P value` else "",
          Sig. = if ("Sig." %in% names(tmp)) tmp$Sig. else "",
          check.names = FALSE
        )
      )
    }

    if (!is.null(disp$dispersion_random)) {
      tmp <- disp$dispersion_random
      term <- if ("Term" %in% names(tmp)) tmp$Term else if ("Term 1" %in% names(tmp)) tmp[["Term 1"]] else ""
      estimate <- if ("Std. Dev." %in% names(tmp)) tmp[["Std. Dev."]] else if ("Std.Dev./Corr" %in% names(tmp)) tmp[["Std.Dev./Corr"]] else NA
      extra <- if ("Correlation" %in% names(tmp)) tmp$Correlation else if ("Term 2" %in% names(tmp)) tmp[["Term 2"]] else ""
      disp_rows <- rbind(
        disp_rows,
        data.frame(
          Section = "Dispersion/autocorrelation",
          Term = paste(tmp$Group, term, extra),
          Estimate = estimate,
          `Std. Error` = "",
          Statistic = "",
          `P value` = "",
          Sig. = "",
          check.names = FALSE
        )
      )
    }
  }

  ran <- extract_random_effects(model)
  ran_rows <- NULL
  if (!is.null(ran)) {
    ran_rows <- data.frame(
      Section = "Random effects",
      Term = paste(ran$Group, ran$Term),
      Estimate = ran$`Std. Dev.`,
      `Std. Error` = "",
      Statistic = "",
      `P value` = "",
      Sig. = "",
      check.names = FALSE
    )
  }

  pieces <- list(fit_rows, cond, zi, disp_rows, ran_rows)
  pieces <- pieces[!vapply(pieces, is.null, logical(1))]

  all_names <- unique(unlist(lapply(pieces, names)))
  pieces <- lapply(pieces, function(x) {
    missing <- setdiff(all_names, names(x))
    for (nm in missing) x[[nm]] <- ""
    x[, all_names, drop = FALSE]
  })

  out <- do.call(rbind, pieces)
  keep <- c("Section", "Term", "Estimate", "Std. Error", "Statistic", "P value", "Sig.")
  out <- out[, intersect(keep, names(out)), drop = FALSE]

  for (nm in names(out)) {
    if (nm == "P value") {
      out[[nm]] <- fmt_docx_numeric_cells(out[[nm]], digits, p_value = TRUE)
    } else if (nm %in% c("Estimate", "Std. Error", "Statistic")) {
      out[[nm]] <- fmt_docx_numeric_cells(out[[nm]], digits, p_value = FALSE)
    } else {
      out[[nm]] <- ifelse(is.na(out[[nm]]), "", as.character(out[[nm]]))
    }
  }

  out
}

export_glmmTMB_docx <- function(models,
                                file = "model_summaries.docx",
                                digits = 3,
                                pretty_names = default_pretty_names(),
                                null_models = NULL,
                                title = "Model summaries") {
  stopifnot(is.list(models))
  if (length(models) == 0) stop("`models` is empty.")
  if (!requireNamespace("officer", quietly = TRUE)) {
    stop("Package 'officer' is required to write .docx files. Install it with install.packages('officer').")
  }
  if (!requireNamespace("flextable", quietly = TRUE)) {
    stop("Package 'flextable' is required to write formatted Word tables. Install it with install.packages('flextable').")
  }
  if (!is.null(null_models) && is.null(names(null_models))) {
    stop("If `null_models` is supplied, it must be a named list.")
  }

  model_names <- names(models)
  if (is.null(model_names)) model_names <- paste0("Model_", seq_along(models))
  model_names[!nzchar(model_names)] <- paste0("Model_", which(!nzchar(model_names)))

  doc <- officer::read_docx()
  doc <- officer::body_add_par(doc, title, style = "heading 1")
  doc <- officer::body_add_par(doc, paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), style = "Normal")
  doc <- officer::body_add_par(doc, "R2 values are reported for transparency and should be compared cautiously across different families or dispersion structures.", style = "Normal")

  results <- vector("list", length(models))
  names(results) <- model_names

  for (i in seq_along(models)) {
    model <- models[[i]]
    model_name <- model_names[[i]]
    null_model <- if (!is.null(null_models) && model_name %in% names(null_models)) null_models[[model_name]] else NULL

    lr2 <- likelihood_r2(model, null_model)
    nr2 <- nakagawa_r2(model)
    tab <- table_for_docx(model, lr2, nr2, pretty_names, digits)

    doc <- officer::body_add_par(doc, model_name, style = "heading 2")

    notes <- c(lr2$note, nr2$note)
    notes <- notes[!is.na(notes) & nzchar(notes)]
    if (length(notes) > 0) {
      doc <- officer::body_add_par(doc, paste("R2 notes:", paste(notes, collapse = " ")), style = "Normal")
    }

    ft <- flextable::flextable(tab)
    ft <- flextable::theme_booktabs(ft)
    ft <- flextable::fontsize(ft, size = 9, part = "all")
    ft <- flextable::bold(ft, part = "header")
    ft <- flextable::align(ft, align = "left", part = "all")
    ft <- flextable::autofit(ft)

    doc <- flextable::body_add_flextable(doc, ft)
    doc <- officer::body_add_par(doc, "", style = "Normal")

    results[[i]] <- tab
  }

  print(doc, target = file)
  message("Wrote Word model summaries to: ", normalizePath(file, winslash = "/", mustWork = FALSE))
  invisible(results)
}

format_sem_coef_table <- function(sem_obj, digits = 3, standardize = "scale") {
  if (!requireNamespace("piecewiseSEM", quietly = TRUE)) {
    stop("Package 'piecewiseSEM' is required to extract path-analysis coefficients.")
  }

  tab <- piecewiseSEM::coefs(sem_obj, standardize = standardize)
  tab <- as.data.frame(tab, stringsAsFactors = FALSE, check.names = FALSE)
  names(tab)[names(tab) == "" | is.na(names(tab))] <- "Sig."

  p_cols <- grepl("^P$|^P[._ ]?Value$|^p[._ ]?value$|Pr\\(", names(tab), ignore.case = TRUE)
  text_cols <- grepl("response|predictor|path|direction|sig|signif", names(tab), ignore.case = TRUE)

  for (j in seq_along(tab)) {
    nm <- names(tab)[[j]]
    if (p_cols[[j]]) {
      tab[[nm]] <- fmt_docx_numeric_cells(tab[[nm]], digits, p_value = TRUE)
    } else if (!text_cols[[j]]) {
      tab[[nm]] <- fmt_docx_numeric_cells(tab[[nm]], digits, p_value = FALSE)
    } else {
      tab[[nm]] <- ifelse(is.na(tab[[nm]]), "", as.character(tab[[nm]]))
    }
  }

  tab
}

export_piecewiseSEM_markdown <- function(sems,
                                         file = "sem_path_results.md",
                                         digits = 3,
                                         standardize = "scale",
                                         title = "Path analysis results") {
  stopifnot(is.list(sems))
  if (length(sems) == 0) stop("`sems` is empty.")

  sem_names <- names(sems)
  if (is.null(sem_names)) sem_names <- paste0("SEM_", seq_along(sems))
  sem_names[!nzchar(sem_names)] <- paste0("SEM_", which(!nzchar(sem_names)))

  con <- file(file, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)

  writeLines(c(
    paste0("# ", title),
    "",
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    paste0("Standardization: ", standardize),
    "",
    "Numeric values are shown to three significant figures by default; p-values below 0.001 are shown as <0.001.",
    "",
    "---",
    ""
  ), con)

  results <- vector("list", length(sems))
  names(results) <- sem_names

  for (i in seq_along(sems)) {
    sem_name <- sem_names[[i]]
    tab <- format_sem_coef_table(sems[[i]], digits = digits, standardize = standardize)

    writeLines(c(paste0("## ", sem_name), ""), con)
    writeLines(pad_markdown_table(tab, digits = digits), con)
    writeLines(c("", "---", ""), con)

    results[[i]] <- tab
  }

  message("Wrote SEM path results to: ", normalizePath(file, winslash = "/", mustWork = FALSE))
  invisible(results)
}

export_piecewiseSEM_docx <- function(sems,
                                     file = "sem_path_results.docx",
                                     digits = 3,
                                     standardize = "scale",
                                     title = "Path analysis results") {
  stopifnot(is.list(sems))
  if (length(sems) == 0) stop("`sems` is empty.")
  if (!requireNamespace("officer", quietly = TRUE)) {
    stop("Package 'officer' is required to write .docx files. Install it with install.packages('officer').")
  }
  if (!requireNamespace("flextable", quietly = TRUE)) {
    stop("Package 'flextable' is required to write formatted Word tables. Install it with install.packages('flextable').")
  }

  sem_names <- names(sems)
  if (is.null(sem_names)) sem_names <- paste0("SEM_", seq_along(sems))
  sem_names[!nzchar(sem_names)] <- paste0("SEM_", which(!nzchar(sem_names)))

  doc <- officer::read_docx()
  doc <- officer::body_add_par(doc, title, style = "heading 1")
  doc <- officer::body_add_par(doc, paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), style = "Normal")
  doc <- officer::body_add_par(doc, paste0("Standardization: ", standardize), style = "Normal")
  doc <- officer::body_add_par(doc, "Numeric values are shown to three significant figures by default; p-values below 0.001 are shown as <0.001.", style = "Normal")

  results <- vector("list", length(sems))
  names(results) <- sem_names

  for (i in seq_along(sems)) {
    sem_name <- sem_names[[i]]
    tab <- format_sem_coef_table(sems[[i]], digits = digits, standardize = standardize)

    doc <- officer::body_add_par(doc, sem_name, style = "heading 2")

    ft <- flextable::flextable(tab)
    ft <- flextable::theme_booktabs(ft)
    ft <- flextable::fontsize(ft, size = 9, part = "all")
    ft <- flextable::bold(ft, part = "header")
    ft <- flextable::align(ft, align = "left", part = "all")
    ft <- flextable::autofit(ft)

    doc <- flextable::body_add_flextable(doc, ft)
    doc <- officer::body_add_par(doc, "", style = "Normal")

    results[[i]] <- tab
  }

  print(doc, target = file)
  message("Wrote Word SEM path results to: ", normalizePath(file, winslash = "/", mustWork = FALSE))
  invisible(results)
}

format_export_table <- function(tab, digits = 3) {
  tab <- as.data.frame(tab, stringsAsFactors = FALSE, check.names = FALSE)
  p_cols <- grepl("^P$|^P[._ ]?Value$|^p[._ ]?value$|Pr\\(", names(tab), ignore.case = TRUE)

  for (j in seq_along(tab)) {
    nm <- names(tab)[[j]]
    if (p_cols[[j]]) {
      tab[[nm]] <- fmt_docx_numeric_cells(tab[[nm]], digits, p_value = TRUE)
    } else if (is.numeric(tab[[nm]])) {
      tab[[nm]] <- fmt_sig(tab[[nm]], digits)
    } else {
      tab[[nm]] <- ifelse(is.na(tab[[nm]]), "", as.character(tab[[nm]]))
    }
  }

  tab
}

get_export_table_n <- function(tab, table_name, table_index, n = NULL) {
  if (!is.null(n)) {
    if (!is.null(names(n)) && table_name %in% names(n)) return(n[[table_name]])
    if (length(n) >= table_index) return(n[[table_index]])
  }

  attr_n <- attr(tab, "n", exact = TRUE)
  if (!is.null(attr_n)) return(attr_n)

  n_cols <- intersect(c("n", "N", "nobs", "Nobs", "Observations"), names(tab))
  if (length(n_cols) > 0) {
    vals <- unique(tab[[n_cols[[1]]]])
    vals <- vals[!is.na(vals)]
    if (length(vals) == 1) return(vals[[1]])
  }

  NULL
}

export_table_list_markdown <- function(tables,
                                       file = "tables.md",
                                       digits = 3,
                                       title = "Model comparison tables",
                                       n = NULL) {
  stopifnot(is.list(tables))
  if (length(tables) == 0) stop("`tables` is empty.")

  table_names <- names(tables)
  if (is.null(table_names)) table_names <- paste0("Table_", seq_along(tables))
  table_names[!nzchar(table_names)] <- paste0("Table_", which(!nzchar(table_names)))

  con <- file(file, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)

  writeLines(c(
    paste0("# ", title),
    "",
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "Numeric values are shown to three significant figures by default; p-values below 0.001 are shown as <0.001.",
    "",
    "---",
    ""
  ), con)

  results <- vector("list", length(tables))
  names(results) <- table_names

  for (i in seq_along(tables)) {
    tab <- format_export_table(tables[[i]], digits = digits)
    table_n <- get_export_table_n(tables[[i]], table_names[[i]], i, n = n)
    writeLines(c(paste0("## ", table_names[[i]]), ""), con)
    if (!is.null(table_n)) writeLines(c(paste0("n = ", table_n), ""), con)
    writeLines(pad_markdown_table(tab, digits = digits), con)
    writeLines(c("", "---", ""), con)
    results[[i]] <- tab
  }

  message("Wrote tables to: ", normalizePath(file, winslash = "/", mustWork = FALSE))
  invisible(results)
}

export_table_list_docx <- function(tables,
                                   file = "tables.docx",
                                   digits = 3,
                                   title = "Model comparison tables",
                                   n = NULL) {
  stopifnot(is.list(tables))
  if (length(tables) == 0) stop("`tables` is empty.")
  if (!requireNamespace("officer", quietly = TRUE)) {
    stop("Package 'officer' is required to write .docx files. Install it with install.packages('officer').")
  }
  if (!requireNamespace("flextable", quietly = TRUE)) {
    stop("Package 'flextable' is required to write formatted Word tables. Install it with install.packages('flextable').")
  }

  table_names <- names(tables)
  if (is.null(table_names)) table_names <- paste0("Table_", seq_along(tables))
  table_names[!nzchar(table_names)] <- paste0("Table_", which(!nzchar(table_names)))

  doc <- officer::read_docx()
  doc <- officer::body_add_par(doc, title, style = "heading 1")
  doc <- officer::body_add_par(doc, paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), style = "Normal")
  doc <- officer::body_add_par(doc, "Numeric values are shown to three significant figures by default; p-values below 0.001 are shown as <0.001.", style = "Normal")

  results <- vector("list", length(tables))
  names(results) <- table_names

  for (i in seq_along(tables)) {
    tab <- format_export_table(tables[[i]], digits = digits)
    table_n <- get_export_table_n(tables[[i]], table_names[[i]], i, n = n)
    doc <- officer::body_add_par(doc, table_names[[i]], style = "heading 2")
    if (!is.null(table_n)) {
      doc <- officer::body_add_par(doc, paste0("n = ", table_n), style = "Normal")
    }

    ft <- flextable::flextable(tab)
    ft <- flextable::theme_booktabs(ft)
    ft <- flextable::fontsize(ft, size = 9, part = "all")
    ft <- flextable::bold(ft, part = "header")
    ft <- flextable::align(ft, align = "left", part = "all")
    ft <- flextable::autofit(ft)

    doc <- flextable::body_add_flextable(doc, ft)
    doc <- officer::body_add_par(doc, "", style = "Normal")
    results[[i]] <- tab
  }

  print(doc, target = file)
  message("Wrote Word tables to: ", normalizePath(file, winslash = "/", mustWork = FALSE))
  invisible(results)
}

get_glmmTMB_predictors <- function(model, component = "cond") {
  f <- try(stats::formula(model, component = component), silent = TRUE)
  if (inherits(f, "try-error")) f <- stats::formula(model)

  fixed_formula <- if (requireNamespace("lme4", quietly = TRUE)) {
    lme4::nobars(f)
  } else {
    f
  }

  term_labels <- attr(stats::terms(fixed_formula), "term.labels")
  term_labels <- term_labels[!grepl("\\|", term_labels)]
  if (length(term_labels) == 0) return(character(0))

  vars <- unique(unlist(lapply(term_labels, function(term) {
    all.vars(stats::as.formula(paste("~", term)))
  })))

  response_vars <- all.vars(fixed_formula[[2]])
  setdiff(vars, response_vars)
}
