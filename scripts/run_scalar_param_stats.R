#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

setwd('/Users/doq/RCode')
load('preprocessed_data.RData')

# Match core preprocessing used in analysisCode.R
alldata[, binary_group := fifelse(grepl('ASD', group, ignore.case = TRUE), 'ASD', 'TD')]
alldata[, group_verbal := fcase(
  is.na(vinabcstd), 'flag',
  binary_group == 'ASD' & vinabcstd >= 75, 'ASD',
  binary_group == 'TD' & vinabcstd >= 75, 'TD',
  binary_group == 'ASD' & vinabcstd < 75, 'ASD (low verbal)',
  binary_group == 'TD' & vinabcstd < 75, 'flag_low_td',
  default = 'Uncategorized'
)]
alldata[, srs_normal := fifelse((binary_group == 'ASD' & srs2total > 20) | (binary_group == 'TD' & srs2total < 120), 1, 0)]
alldata[, above_crit := fifelse((nTrials == 200 & perf > .616), 1, 0)]

# Same inclusion window as Figure 4 scalar model block
DT <- alldata[
  game_version == 'ft' & srs_normal == 1 & above_crit == 1 &
    group_verbal != 'flag' & group_verbal != 'flag_low_td' &
    trial %in% 11:200,
  .(binary_group, choice, numflashright, numflashleft, winstay, loseswitch)
]
DT[, `:=`(r = numflashright, l = numflashleft)]

negloglik_scalar <- function(params, data) {
  # params: k0, k1, bs, b1, b2
  k0 <- params[1]
  k1 <- params[2]
  bs <- params[3]
  b1 <- params[4]
  b2 <- params[5]

  sd_term <- sqrt(k1 * (data$r^2 + data$l^2) + k0)
  p_left <- pnorm(0, mean = data$r - data$l + bs + b1 * data$winstay + b2 * data$loseswitch, sd = sd_term)
  p_right <- 1 - p_left

  p_obs <- ifelse(data$choice == 1, p_right, p_left)
  p_obs[p_obs < 1e-6] <- 1e-6
  -sum(log(p_obs))
}

fit_scalar_group <- function(data) {
  # L-BFGS-B is deterministic and much faster than SANN for this summary analysis
  fit <- optim(
    par = c(0.5, 0.1, 0, 0, 0),
    fn = negloglik_scalar,
    data = data,
    method = 'L-BFGS-B',
    lower = c(1e-6, 1e-6, -Inf, -Inf, -Inf),
    upper = c(Inf, Inf, Inf, Inf, Inf),
    hessian = TRUE,
    control = list(maxit = 4000)
  )

  se <- rep(NA_real_, 5)
  inv_h <- tryCatch(solve(fit$hessian), error = function(e) NULL)
  if (!is.null(inv_h)) {
    se <- sqrt(diag(inv_h))
  }

  list(par = fit$par, se = se, convergence = fit$convergence, nll = fit$value)
}

asd_data <- as.data.frame(DT[binary_group == 'ASD'])
td_data <- as.data.frame(DT[binary_group == 'TD'])

fit_asd <- fit_scalar_group(asd_data)
fit_td <- fit_scalar_group(td_data)

params <- c('k0', 'k1', 'bs', 'b1', 'b2')
res <- data.table(
  parameter = params,
  ASD = fit_asd$par,
  ASD_se = fit_asd$se,
  TD = fit_td$par,
  TD_se = fit_td$se
)

res[, diff := ASD - TD]
res[, se_diff := sqrt(ASD_se^2 + TD_se^2)]
res[, z := diff / se_diff]
res[, p := 2 * pnorm(abs(z), lower.tail = FALSE)]
res[, ci_low := diff - 1.96 * se_diff]
res[, ci_high := diff + 1.96 * se_diff]
res[, `:=`(ASD_convergence = fit_asd$convergence, TD_convergence = fit_td$convergence)]

out_dir <- '/Users/doq/RCode/results/stats_additional'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_csv <- file.path(out_dir, 'scalar_sdt_group_param_differences.csv')
fwrite(res, out_csv)

cat('\nScalar SDT parameter contrasts (ASD - TD)\n')
print(res)
cat('\nSaved:', out_csv, '\n')

fmt <- function(x, d = 4) format(round(x, d), nsmall = d)
for (row in 1:nrow(res)) {
  r <- res[row]
  cat(sprintf(
    "%s: ASD=%s, TD=%s; Delta=%s, 95%% CI [%s, %s], z=%s, p=%s\n",
    r$parameter,
    fmt(r$ASD), fmt(r$TD),
    fmt(r$diff), fmt(r$ci_low), fmt(r$ci_high),
    fmt(r$z, 2), format(r$p, scientific = TRUE, digits = 4)
  ))
}
