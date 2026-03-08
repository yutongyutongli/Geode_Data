#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------
# Utility functions
# -----------------------------
cohen_d <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  n1 <- length(x)
  n2 <- length(y)
  if (n1 < 2 || n2 < 2) return(NA_real_)
  s1 <- stats::sd(x)
  s2 <- stats::sd(y)
  sp <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
  if (!is.finite(sp) || sp == 0) return(NA_real_)
  (mean(x) - mean(y)) / sp
}

calc_basic_effects <- function(dt, outcome, group_col = "binary_group", g1 = "ASD", g2 = "TD") {
  sub <- dt[!is.na(get(outcome)) & get(group_col) %in% c(g1, g2)]
  x <- sub[get(group_col) == g1, get(outcome)]
  y <- sub[get(group_col) == g2, get(outcome)]

  test <- tryCatch(stats::t.test(x, y), error = function(e) NULL)

  data.table(
    outcome = outcome,
    group1 = g1,
    n1 = length(x),
    mean1 = mean(x),
    sd1 = stats::sd(x),
    group2 = g2,
    n2 = length(y),
    mean2 = mean(y),
    sd2 = stats::sd(y),
    cohen_d = cohen_d(x, y),
    t_p = if (!is.null(test)) test$p.value else NA_real_
  )
}

confusion_metrics <- function(y_true, y_pred) {
  y_true <- as.integer(y_true)
  y_pred <- as.integer(y_pred)

  tp <- sum(y_true == 1 & y_pred == 1)
  tn <- sum(y_true == 0 & y_pred == 0)
  fp <- sum(y_true == 0 & y_pred == 1)
  fn <- sum(y_true == 1 & y_pred == 0)

  acc <- (tp + tn) / (tp + tn + fp + fn)
  prec <- if ((tp + fp) == 0) NA_real_ else tp / (tp + fp)
  rec <- if ((tp + fn) == 0) NA_real_ else tp / (tp + fn)
  f1 <- if (is.na(prec) || is.na(rec) || (prec + rec) == 0) NA_real_ else 2 * prec * rec / (prec + rec)
  spec <- if ((tn + fp) == 0) NA_real_ else tn / (tn + fp)
  bal_acc <- mean(c(rec, spec), na.rm = TRUE)

  # TD as positive class
  prec_td <- if ((tn + fn) == 0) NA_real_ else tn / (tn + fn)
  rec_td <- if ((tn + fp) == 0) NA_real_ else tn / (tn + fp)
  f1_td <- if (is.na(prec_td) || is.na(rec_td) || (prec_td + rec_td) == 0) NA_real_ else 2 * prec_td * rec_td / (prec_td + rec_td)
  macro_f1 <- mean(c(f1, f1_td), na.rm = TRUE)

  data.table(
    accuracy = acc,
    balanced_accuracy = bal_acc,
    f1 = f1,
    macro_f1 = macro_f1,
    precision = prec,
    recall = rec,
    specificity = spec,
    tp = tp,
    tn = tn,
    fp = fp,
    fn = fn
  )
}

stratified_folds <- function(y, k = 5, seed = 123) {
  set.seed(seed)
  idx1 <- which(y == 1)
  idx0 <- which(y == 0)
  idx1 <- sample(idx1)
  idx0 <- sample(idx0)

  folds <- vector("list", k)
  for (i in seq_len(k)) {
    folds[[i]] <- c(idx1[seq(i, length(idx1), by = k)], idx0[seq(i, length(idx0), by = k)])
  }
  folds
}

cv_logit_perf <- function(dt, repeats = 10, k = 5) {
  # binary_group: ASD=1, TD=0
  dat <- dt[!is.na(perf) & binary_group %in% c("ASD", "TD"), .(perf, binary_group)]
  dat[, y := as.integer(binary_group == "ASD")]

  all_pred <- integer(0)
  all_true <- integer(0)

  for (r in seq_len(repeats)) {
    folds <- stratified_folds(dat$y, k = k, seed = 1000 + r)
    for (i in seq_len(k)) {
      test_idx <- folds[[i]]
      train_idx <- setdiff(seq_len(nrow(dat)), test_idx)

      train <- dat[train_idx]
      test <- dat[test_idx]

      mdl <- stats::glm(y ~ perf, data = train, family = stats::binomial())
      p_hat <- stats::predict(mdl, newdata = test, type = "response")
      y_hat <- as.integer(p_hat >= 0.5)

      all_pred <- c(all_pred, y_hat)
      all_true <- c(all_true, test$y)
    }
  }

  confusion_metrics(all_true, all_pred)
}

cv_logit_multifeature <- function(dt, features, repeats = 10, k = 5) {
  keep <- c("binary_group", features)
  dat <- dt[, ..keep]
  dat <- dat[complete.cases(dat)]
  dat <- dat[binary_group %in% c("ASD", "TD")]
  dat[, y := as.integer(binary_group == "ASD")]

  if (nrow(dat) < 30) stop("Too few complete cases for multifeature classifier.")

  all_pred <- integer(0)
  all_true <- integer(0)

  fml <- stats::as.formula(paste("y ~", paste(features, collapse = " + ")))

  for (r in seq_len(repeats)) {
    folds <- stratified_folds(dat$y, k = k, seed = 2000 + r)
    for (i in seq_len(k)) {
      test_idx <- folds[[i]]
      train_idx <- setdiff(seq_len(nrow(dat)), test_idx)
      train <- dat[train_idx]
      test <- dat[test_idx]

      mdl <- stats::glm(fml, data = train, family = stats::binomial())
      p_hat <- stats::predict(mdl, newdata = test, type = "response")
      y_hat <- as.integer(p_hat >= 0.5)

      all_pred <- c(all_pred, y_hat)
      all_true <- c(all_true, test$y)
    }
  }

  out <- confusion_metrics(all_true, all_pred)
  out[, n_complete := nrow(dat)]
  out
}

cor_table <- function(dt, x, y, group_col) {
  out <- dt[!is.na(get(x)) & !is.na(get(y)), {
    if (.N < 3) {
      list(n = .N, r = NA_real_, p = NA_real_)
    } else {
      ct <- stats::cor.test(get(x), get(y), method = "pearson")
      list(n = .N, r = unname(ct$estimate), p = ct$p.value)
    }
  }, by = group_col]
  setnames(out, group_col, "group")
  out[, `:=`(x = x, y = y)]
  out[, grouping := group_col]
  out[, .(x, y, grouping, group, n, r, p)]
}

# -----------------------------
# Load and preprocess
# -----------------------------
setwd("/Users/doq/RCode")
load("preprocessed_data.RData")

# mirror analysisCode.R logic
alldata[, binary_group := fifelse(grepl("ASD", group, ignore.case = TRUE), "ASD", "TD")]

alldata[, group_verbal := fcase(
  is.na(vinabcstd), "flag",
  binary_group == "ASD" & vinabcstd >= 75, "ASD",
  binary_group == "TD" & vinabcstd >= 75, "TD",
  binary_group == "ASD" & vinabcstd < 75, "ASD (low verbal)",
  binary_group == "TD" & vinabcstd < 75, "flag_low_td",
  default = "Uncategorized"
)]

alldata[, srs_normal := fifelse(
  (binary_group == "ASD" & srs2total > 20) | (binary_group == "TD" & srs2total < 120), 1, 0
)]

alldata[, above_crit := fifelse((nTrials == 200 & perf > 0.616), 1, 0)]

subject_summary <- unique(alldata[, .(
  username, subID, game_version, source, group, binary_group, group_verbal, srs_normal,
  vinabcstd, srs2total, bistotal,
  aasp_low_reg_raw, aasp_sen_seek_raw, aasp_sen_ses_raw, aasp_sen_avoid_raw,
  perf, trainperf, nTrials, above_crit
)])

# Main manuscript-like subset
main_dt <- subject_summary[
  game_version == "ft" & srs_normal == 1 & above_crit == 1 & group_verbal != "flag" & group_verbal != "flag_low_td"
]

# Classification subset for ASD-vs-TD claim (not conditioning on above_crit)
class_dt <- subject_summary[
  game_version == "ft" & srs_normal == 1 & binary_group %in% c("ASD", "TD")
]

# Derive additional game features (subject level) to match manuscript metrics
trial_dt <- alldata[
  game_version == "ft" & srs_normal == 1 & above_crit == 1 &
    group_verbal != "flag" & group_verbal != "flag_low_td" &
    trial %in% 11:200,
  .(subID, choice, dflash, winstay, loseswitch, iti, rt)
]

fit_sub_glm <- function(df_sub) {
  out <- tryCatch({
    mdl <- stats::glm(choice ~ dflash + winstay + loseswitch, data = df_sub, family = stats::binomial())
    cf <- stats::coef(mdl)
    data.table(
      Sidebias = unname(cf["(Intercept)"]),
      Slope = unname(cf["dflash"]),
      WS = unname(cf["winstay"]),
      LS = unname(cf["loseswitch"])
    )
  }, error = function(e) {
    data.table(Sidebias = NA_real_, Slope = NA_real_, WS = NA_real_, LS = NA_real_)
  })
  out
}

glm_features <- trial_dt[, fit_sub_glm(.SD), by = subID]
timing_features <- trial_dt[, .(mean_iti = mean(iti, na.rm = TRUE), mean_rt = mean(rt, na.rm = TRUE)), by = subID]

main_dt <- merge(main_dt, glm_features, by = "subID", all.x = TRUE)
main_dt <- merge(main_dt, timing_features, by = "subID", all.x = TRUE)
main_dt[, `:=`(
  abs_Sidebias = abs(Sidebias),
  abs_WS = abs(WS),
  abs_LS = abs(LS)
)]

# -----------------------------
# 1) Groupwise effect sizes
# -----------------------------
outcomes <- c("perf", "trainperf", "Slope", "abs_Sidebias", "abs_WS", "abs_LS", "mean_iti", "mean_rt")
effects <- rbindlist(lapply(outcomes, function(v) calc_basic_effects(main_dt, v)))

# Add three-level pairwise effect sizes for perf (TD vs ASD, TD vs ASD low verbal, ASD vs ASD low verbal)
pair_groups <- list(
  c("TD", "ASD"),
  c("TD", "ASD (low verbal)"),
  c("ASD", "ASD (low verbal)")
)

pairwise_perf <- rbindlist(lapply(pair_groups, function(gp) {
  sub <- main_dt[group_verbal %in% gp & !is.na(perf)]
  x <- sub[group_verbal == gp[1], perf]
  y <- sub[group_verbal == gp[2], perf]
  test <- tryCatch(stats::t.test(x, y), error = function(e) NULL)
  data.table(
    outcome = "perf",
    group1 = gp[1],
    n1 = length(x),
    mean1 = mean(x),
    sd1 = stats::sd(x),
    group2 = gp[2],
    n2 = length(y),
    mean2 = mean(y),
    sd2 = stats::sd(y),
    cohen_d = cohen_d(x, y),
    t_p = if (!is.null(test)) test$p.value else NA_real_
  )
}))

# -----------------------------
# 2) Classifier metrics (accuracy/F1)
# -----------------------------
class_metrics <- cv_logit_perf(class_dt, repeats = 20, k = 5)
majority_baseline <- max(
  class_dt[binary_group == "ASD", .N],
  class_dt[binary_group == "TD", .N]
) / nrow(class_dt)
class_metrics[, `:=`(
  n = nrow(class_dt),
  asd_n = class_dt[binary_group == "ASD", .N],
  td_n = class_dt[binary_group == "TD", .N],
  majority_class_accuracy = majority_baseline,
  model = "5-fold CV logistic regression: binary_group ~ perf (20 repeats)"
)]

multi_features <- c("perf", "trainperf", "Slope", "abs_Sidebias", "abs_WS", "abs_LS", "mean_iti")
class_metrics_multi <- cv_logit_multifeature(main_dt[, c("binary_group", multi_features), with = FALSE], multi_features, repeats = 20, k = 5)
class_metrics_multi[, `:=`(
  n = nrow(main_dt),
  asd_n = main_dt[binary_group == "ASD", .N],
  td_n = main_dt[binary_group == "TD", .N],
  majority_class_accuracy = max(main_dt[binary_group == "ASD", .N], main_dt[binary_group == "TD", .N]) / nrow(main_dt),
  model = paste0("5-fold CV logistic regression: binary_group ~ ", paste(multi_features, collapse = " + "), " (20 repeats)")
)]

# -----------------------------
# 3) Correlations: pooled vs within-group
# -----------------------------
corr_specs <- list(
  c("srs2total", "perf"),
  c("vinabcstd", "perf")
)

corr_pooled <- rbindlist(lapply(corr_specs, function(sp) {
  x <- sp[1]
  y <- sp[2]
  sub <- main_dt[!is.na(get(x)) & !is.na(get(y))]
  ct <- stats::cor.test(sub[[x]], sub[[y]], method = "pearson")
  data.table(x = x, y = y, grouping = "pooled", group = "Pooled", n = nrow(sub), r = unname(ct$estimate), p = ct$p.value)
}))

corr_binary <- rbindlist(lapply(corr_specs, function(sp) cor_table(main_dt, sp[1], sp[2], "binary_group")))
corr_verbal <- rbindlist(lapply(corr_specs, function(sp) cor_table(main_dt, sp[1], sp[2], "group_verbal")))
correlations <- rbind(corr_pooled, corr_binary, corr_verbal, fill = TRUE)

# -----------------------------
# Save outputs
# -----------------------------
out_dir <- "results/stats_additional"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fwrite(effects, file.path(out_dir, "effect_sizes_binary_groups.csv"))
fwrite(pairwise_perf, file.path(out_dir, "effect_sizes_group_verbal_perf.csv"))
fwrite(class_metrics, file.path(out_dir, "classifier_metrics_asd_vs_td.csv"))
fwrite(class_metrics_multi, file.path(out_dir, "classifier_metrics_asd_vs_td_multifeature.csv"))
fwrite(correlations, file.path(out_dir, "pooled_vs_withingroup_correlations.csv"))

# concise console summary
cat("\n=== Additional Stats Summary ===\n")
cat("Main analysis N:", nrow(main_dt), "\n")
cat("Classifier dataset N:", nrow(class_dt), "(ASD:", class_dt[binary_group == "ASD", .N], ", TD:", class_dt[binary_group == "TD", .N], ")\n")
cat("\n[Effect sizes: binary groups]\n")
print(effects)
cat("\n[Effect sizes: group_verbal pairwise for perf]\n")
print(pairwise_perf)
cat("\n[Classifier metrics]\n")
print(class_metrics)
cat("\n[Classifier metrics - multifeature]\n")
print(class_metrics_multi)
cat("\n[Correlations: pooled and within-group]\n")
print(correlations)

# short markdown report for manuscript drafting
report_path <- file.path(out_dir, "reviewer_stats_response.md")

fmt_num <- function(x, k = 3) ifelse(is.na(x), "NA", format(round(x, k), nsmall = k))

pooled_srs_perf <- correlations[x == "srs2total" & y == "perf" & grouping == "pooled" & group == "Pooled"]
asd_srs_perf <- correlations[x == "srs2total" & y == "perf" & grouping == "binary_group" & group == "ASD"]
td_srs_perf <- correlations[x == "srs2total" & y == "perf" & grouping == "binary_group" & group == "TD"]
asd_low_srs_perf <- correlations[x == "srs2total" & y == "perf" & grouping == "group_verbal" & group == "ASD (low verbal)"]

perf_td_asd <- pairwise_perf[group1 == "TD" & group2 == "ASD"]
perf_td_asdlow <- pairwise_perf[group1 == "TD" & group2 == "ASD (low verbal)"]
slope_td_asd <- effects[outcome == "Slope" & group1 == "ASD" & group2 == "TD"]
kline <- if (nrow(slope_td_asd) == 1) {
  paste0("- `Psychometric slope`: ASD vs TD mean(SD) = ",
         fmt_num(slope_td_asd$mean1), " (", fmt_num(slope_td_asd$sd1), ") vs ",
         fmt_num(slope_td_asd$mean2), " (", fmt_num(slope_td_asd$sd2), "), Cohen's d = ",
         fmt_num(slope_td_asd$cohen_d), ".")
} else {
  "- `Psychometric slope`: added in effect-size table."
}

lines <- c(
  "# Additional Statistics Requested by Reviewer",
  "",
  "This file was generated by `additional_stats_geode.R` without modifying `analysisCode.R`.",
  "",
  "## 1) Groupwise effect size (with SD and Cohen's d)",
  "",
  paste0("- `perf`: TD vs ASD: mean(SD) = ", fmt_num(perf_td_asd$mean1), " (", fmt_num(perf_td_asd$sd1), ") vs ",
         fmt_num(perf_td_asd$mean2), " (", fmt_num(perf_td_asd$sd2), "), Cohen's d = ", fmt_num(perf_td_asd$cohen_d), "."),
  paste0("- `perf`: TD vs ASD (low verbal): mean(SD) = ", fmt_num(perf_td_asdlow$mean1), " (", fmt_num(perf_td_asdlow$sd1), ") vs ",
         fmt_num(perf_td_asdlow$mean2), " (", fmt_num(perf_td_asdlow$sd2), "), Cohen's d = ", fmt_num(perf_td_asdlow$cohen_d), "."),
  kline,
  "",
  "Full effect size tables:",
  "- `effect_sizes_binary_groups.csv`",
  "- `effect_sizes_group_verbal_perf.csv`",
  "",
  "## 2) Classifier accuracy/F1 (ASD vs TD)",
  "",
  paste0("- Model: ", class_metrics$model[1]),
  paste0("- Accuracy = ", fmt_num(class_metrics$accuracy[1]), ", F1 = ", fmt_num(class_metrics$f1[1]),
         ", Macro-F1 = ", fmt_num(class_metrics$macro_f1[1]), "."),
  paste0("- Balanced accuracy = ", fmt_num(class_metrics$balanced_accuracy[1]),
         ", Sensitivity(ASD) = ", fmt_num(class_metrics$recall[1]),
         ", Specificity(TD) = ", fmt_num(class_metrics$specificity[1]), "."),
  paste0("- Majority-class baseline accuracy = ", fmt_num(class_metrics$majority_class_accuracy[1]), "."),
  paste0("- Multifeature model (perf + trainperf + slope/history/timing): Accuracy = ", fmt_num(class_metrics_multi$accuracy[1]),
         ", Macro-F1 = ", fmt_num(class_metrics_multi$macro_f1[1]),
         ", Balanced accuracy = ", fmt_num(class_metrics_multi$balanced_accuracy[1]), "."),
  "",
  "Full classifier table:",
  "- `classifier_metrics_asd_vs_td.csv`",
  "- `classifier_metrics_asd_vs_td_multifeature.csv`",
  "",
  "## 3) Within-group correlations (to avoid pooled-group conflation)",
  "",
  paste0("- Pooled `srs2total` vs `perf`: r = ", fmt_num(pooled_srs_perf$r), ", p = ", fmt_num(pooled_srs_perf$p), ", n = ", pooled_srs_perf$n, "."),
  paste0("- ASD-only `srs2total` vs `perf`: r = ", fmt_num(asd_srs_perf$r), ", p = ", fmt_num(asd_srs_perf$p), ", n = ", asd_srs_perf$n, "."),
  paste0("- ASD (low verbal) `srs2total` vs `perf`: r = ", fmt_num(asd_low_srs_perf$r), ", p = ", fmt_num(asd_low_srs_perf$p), ", n = ", asd_low_srs_perf$n, "."),
  paste0("- TD-only `srs2total` vs `perf`: r = ", fmt_num(td_srs_perf$r), ", p = ", fmt_num(td_srs_perf$p), ", n = ", td_srs_perf$n, "."),
  "",
  "Full correlation table:",
  "- `pooled_vs_withingroup_correlations.csv`"
)

writeLines(lines, report_path)
cat("\nWrote report:", report_path, "\n")
