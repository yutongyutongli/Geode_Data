# Three-group bootstrap learning curve analysis with KS tests (ASD_high, ASD_low, TD)
# Hard phase: trial > 10

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(minpack.lm)
})

load("/Users/doq/Geode_Data/preprocessed_data.RData")

alldata <- as.data.table(alldata)

# ft only
alldata <- alldata[game_version == "ft"]

# binary group
alldata[, binary_group := fifelse(grepl("ASD", group, ignore.case = TRUE), "ASD", "TD")]

# SRS typical
alldata[, srs_typical := fifelse(
  (binary_group == "ASD" & srs2total >= 20) | (binary_group == "TD" & srs2total <= 120), 1, 0
)]

# task success: nTrials==200 and mean accuracy trials 11-200 >= 0.616
acc_11_200 <- alldata[trial >= 11 & trial <= 200,
                      .(acc_11_200 = mean(score, na.rm = TRUE)),
                      by = .(subID)]
subsum <- unique(alldata[, .(subID, binary_group, vinabcstd, nTrials, srs_typical)])
subsum <- merge(subsum, acc_11_200, by = "subID", all.x = TRUE)
subsum[, task_success := fifelse(nTrials == 200 & acc_11_200 >= 0.616, 1, 0)]

eligible_ids <- subsum[srs_typical == 1 & task_success == 1, subID]

hard_df <- alldata[subID %in% eligible_ids & trial > 10,
                   .(subID, trial, score, binary_group, vinabcstd)]

# group split by VABS (exclude ASD missing VABS for split)
hard_df[, group3 := fifelse(binary_group == "TD", "TD",
                            fifelse(!is.na(vinabcstd) & vinabcstd >= 75, "ASD_high", "ASD_low"))]

hard_df <- hard_df[!(binary_group == "ASD" & is.na(vinabcstd))]

# hard-phase trial index
hard_df[, trial_hard := trial - min(trial) + 1L]

fit_learning_curve <- function(data) {
  tryCatch({
    model <- nlsLM(
      mean_score ~ A0 + (P - A0) * (1 - exp(L * trial_hard)),
      data = data,
      start = list(A0 = 0.5, P = 0.9, L = -0.05),
      lower = c(A0 = 0, P = 0, L = -10),
      upper = c(A0 = 1, P = 1, L = 0),
      control = nls.lm.control(maxiter = 2000)
    )
    as.list(coef(model))
  }, error = function(e) {
    list(A0 = NA_real_, P = NA_real_, L = NA_real_)
  })
}

filter_fits <- function(df) {
  df %>%
    filter(
      !is.na(A0), !is.na(P), !is.na(L),
      L < -1e-6, L > -10 + 1e-6,
      A0 >= 0, A0 <= 1,
      P  >= 0, P  <= 1
    )
}

bootstrap_curve_fit <- function(data, n_boot = 1000, seed = 42) {
  set.seed(seed)
  results <- vector("list", n_boot)

  for (i in seq_len(n_boot)) {
    boot_data <- data %>%
      group_by(group3) %>%
      group_modify(~ {
        ids <- unique(.x$subID)
        boot_ids <- sample(ids, size = length(ids), replace = TRUE)
        bind_rows(lapply(boot_ids, function(id) filter(.x, subID == id)))
      }) %>%
      ungroup()

    boot_avg <- boot_data %>%
      group_by(group3, trial_hard) %>%
      summarise(mean_score = mean(score), .groups = "drop")

    fits <- boot_avg %>%
      group_by(group3) %>%
      nest() %>%
      mutate(params = map(data, fit_learning_curve)) %>%
      unnest_wider(params) %>%
      select(-data)

    fits$iter <- i
    results[[i]] <- fits
  }

  bind_rows(results)
}

boot_params_raw <- bootstrap_curve_fit(hard_df, n_boot = 1000, seed = 42)
boot_params <- filter_fits(boot_params_raw)

# Summary statistics
summary_tbl <- boot_params %>%
  pivot_longer(cols = c(A0, P, L), names_to = "param", values_to = "value") %>%
  group_by(group3, param) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    lower = quantile(value, 0.025, na.rm = TRUE),
    upper = quantile(value, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_tbl)

# KS tests pairwise
pairs <- list(c("ASD_high", "ASD_low"), c("ASD_high", "TD"), c("ASD_low", "TD"))

ks_out <- lapply(pairs, function(p) {
  list(
    pair = p,
    A0 = ks.test(
      boot_params$A0[boot_params$group3 == p[1]],
      boot_params$A0[boot_params$group3 == p[2]]
    )$p.value,
    P = ks.test(
      boot_params$P[boot_params$group3 == p[1]],
      boot_params$P[boot_params$group3 == p[2]]
    )$p.value,
    L = ks.test(
      boot_params$L[boot_params$group3 == p[1]],
      boot_params$L[boot_params$group3 == p[2]]
    )$p.value
  )
})

for (x in ks_out) {
  cat(paste0("Pair ", x$pair[1], " vs ", x$pair[2], "\n"))
  cat(paste0("  A0 p=", format(x$A0, digits = 4), "\n"))
  cat(paste0("  P  p=", format(x$P, digits = 4), "\n"))
  cat(paste0("  L  p=", format(x$L, digits = 4), "\n"))
}
