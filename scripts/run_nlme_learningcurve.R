# Nonlinear mixed-effects learning-curve models (p0, alpha) by group
# Fits separate models for trials 1-10 and 11-200

suppressPackageStartupMessages({
  library(data.table)
  library(nlme)
})

load("/Users/doq/Geode_Data/preprocessed_data.RData")

alldata <- as.data.table(alldata)

# Group definitions
alldata[, binary_group := fifelse(grepl("ASD", group, ignore.case = TRUE), "ASD", "TD")]
alldata[, group_verbal := fcase(
  is.na(vinabcstd), "flag",
  binary_group == "ASD" & vinabcstd >= 75, "ASD",
  binary_group == "TD" & vinabcstd >= 75, "TD",
  binary_group == "ASD" & vinabcstd < 75, "ASD (low verbal)",
  binary_group == "TD" & vinabcstd < 75, "flag_low_td",
  default = "Uncategorized"
)]

# Inclusion criteria
alldata[, srs_normal := fifelse(
  (binary_group == "ASD" & srs2total > 20) | (binary_group == "TD" & srs2total < 120), 1, 0
)]
alldata[, above_crit := fifelse((nTrials == 200 & perf > .616), 1, 0)]

verbal_groups <- c("TD", "ASD", "ASD (low verbal)")

df <- alldata[game_version == "ft" & srs_normal == 1 & above_crit == 1]
df <- df[group_verbal %in% verbal_groups]
df <- df[, .(subID, trial, score, group_verbal)]
df[, group_verbal := factor(group_verbal, levels = verbal_groups)]

run_nlme <- function(dat, label, alpha_start = 0.05) {
  dat <- copy(dat)
  dat <- dat[, .(mean_score = mean(score)), by = .(subID, group_verbal, trial)]
  dat[, trial0 := trial - min(trial) + 1]

  nl_fun <- mean_score ~ 1 - (1 - p0) * exp(-alpha * trial0)

  fit <- nlme(
    nl_fun,
    data = dat,
    fixed = list(p0 ~ group_verbal, alpha ~ group_verbal),
    random = pdDiag(p0 + alpha ~ 1),
    groups = ~ subID,
    start = c(
      p0 = rep(0.7, length(levels(dat$group_verbal))),
      alpha = rep(alpha_start, length(levels(dat$group_verbal)))
    ),
    na.action = na.exclude,
    control = nlmeControl(maxIter = 200, pnlsMaxIter = 50, msMaxIter = 200, msVerbose = FALSE)
  )

  cat("\n=====", label, "(Nonlinear Mixed Effects, pinf=1) =====\n")
  print(summary(fit)$tTable)
  cat("\nRandom effects SD:\n")
  print(VarCorr(fit))
}

# Trials 1-10
run_nlme(df[trial %in% 1:10], "Trials 1-10", alpha_start = 0.1)

# Trials 11-200
run_nlme(df[trial %in% 11:200], "Trials 11-200", alpha_start = 0.01)
