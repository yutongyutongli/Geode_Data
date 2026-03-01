# Bootstrapped learning curve fit (ASD vs TD) for hard phase
# Model: y = A0 + (P - A0) * (1 - exp(L * trial))

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(minpack.lm)
})

# Load data
load("/Users/doq/RCode/preprocessed_data.RData")

# Prepare data
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
alldata[, above_crit := fifelse((nTrials == 200 & perf > .616), 1, 0)]

# Filter to analysis set and hard phase
hard_df <- alldata[game_version == "ft" & srs_normal == 1 & above_crit == 1 &
                     group_verbal != "flag" & group_verbal != "flag_low_td" &
                     trial > 10,
                   .(subID, trial, score, group_verbal)]

# Collapse ASD verbal groups
hard_df <- hard_df %>%
  mutate(group = ifelse(grepl("ASD", group_verbal), "ASD", "TD"))

# Hard-phase trial index starting at 1
hard_df <- hard_df %>%
  group_by() %>%
  mutate(trial_hard = trial - min(trial) + 1L) %>%
  ungroup()

# Fit learning curve to mean accuracy per trial
fit_learning_curve <- function(data) {
  # data: columns trial_hard, mean_score
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

# Filter out failed fits or boundary-hitting solutions
filter_fits <- function(df) {
  df %>%
    filter(
      !is.na(A0), !is.na(P), !is.na(L),
      L < -1e-6, L > -10 + 1e-6,
      A0 >= 0, A0 <= 1,
      P  >= 0, P  <= 1
    )
}

# Bootstrap function
bootstrap_curve_fit <- function(data, n_boot = 1000, seed = 42) {
  set.seed(seed)
  results <- vector("list", n_boot)

  for (i in seq_len(n_boot)) {
    boot_data <- data %>%
      group_by(group) %>%
      group_modify(~ {
        ids <- unique(.x$subID)
        boot_ids <- sample(ids, size = length(ids), replace = TRUE)
        bind_rows(lapply(boot_ids, function(id) filter(.x, subID == id)))
      }) %>%
      ungroup()

    boot_avg <- boot_data %>%
      group_by(group, trial_hard) %>%
      summarise(mean_score = mean(score), .groups = "drop")

    fits <- boot_avg %>%
      group_by(group) %>%
      nest() %>%
      mutate(params = map(data, fit_learning_curve)) %>%
      unnest_wider(params) %>%
      select(-data)

    fits$iter <- i
    results[[i]] <- fits
  }

  bind_rows(results)
}

# Run bootstrap
boot_params_raw <- bootstrap_curve_fit(hard_df, n_boot = 1000, seed = 42)
boot_params <- filter_fits(boot_params_raw)

# Summary statistics
boot_summary <- boot_params %>%
  pivot_longer(cols = c(A0, P, L), names_to = "param", values_to = "value") %>%
  group_by(group, param) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    lower = quantile(value, 0.025, na.rm = TRUE),
    upper = quantile(value, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

print(boot_summary)

# KS tests for group differences
ks_results <- list(
  A0 = ks.test(boot_params$A0[boot_params$group == "ASD"],
               boot_params$A0[boot_params$group == "TD"]),
  P  = ks.test(boot_params$P[boot_params$group == "ASD"],
              boot_params$P[boot_params$group == "TD"]),
  L  = ks.test(boot_params$L[boot_params$group == "ASD"],
              boot_params$L[boot_params$group == "TD"])
)

ks_pvals <- sapply(ks_results, function(x) x$p.value)
print(ks_pvals)

# Bootstrap difference distributions (TD - ASD)
boot_wide <- boot_params %>%
  select(iter, group, A0, P, L) %>%
  pivot_wider(names_from = group, values_from = c(A0, P, L))

boot_diff <- boot_wide %>%
  transmute(
    A0_diff = A0_TD - A0_ASD,
    P_diff  = P_TD  - P_ASD,
    L_diff  = L_TD  - L_ASD
  )

boot_diff_summary <- boot_diff %>%
  pivot_longer(cols = everything(), names_to = "param", values_to = "diff") %>%
  group_by(param) %>%
  summarise(
    mean = mean(diff, na.rm = TRUE),
    sd = sd(diff, na.rm = TRUE),
    lower = quantile(diff, 0.025, na.rm = TRUE),
    upper = quantile(diff, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

boot_diff_pvals <- boot_diff %>%
  summarise(
    p_A0 = 2 * min(mean(A0_diff < 0, na.rm = TRUE), mean(A0_diff > 0, na.rm = TRUE)),
    p_P  = 2 * min(mean(P_diff  < 0, na.rm = TRUE), mean(P_diff  > 0, na.rm = TRUE)),
    p_L  = 2 * min(mean(L_diff  < 0, na.rm = TRUE), mean(L_diff  > 0, na.rm = TRUE))
  )

print(boot_diff_summary)
print(boot_diff_pvals)

# Permutation test on participant labels (TD vs ASD)
perm_test <- function(data, n_perm = 10000, seed = 42) {
  set.seed(seed)
  compute_diffs <- function(df) {
    out <- df %>%
      group_by(group, trial_hard) %>%
      summarise(mean_score = mean(score), .groups = "drop") %>%
      group_by(group) %>%
      nest() %>%
      mutate(params = map(data, fit_learning_curve)) %>%
      unnest_wider(params) %>%
      select(-data) %>%
      filter_fits() %>%
      pivot_wider(names_from = group, values_from = c(A0, P, L))

    if (!all(c("A0_TD", "A0_ASD", "P_TD", "P_ASD", "L_TD", "L_ASD") %in% names(out))) {
      return(tibble(A0_diff = NA_real_, P_diff = NA_real_, L_diff = NA_real_))
    }

    out %>% transmute(A0_diff = A0_TD - A0_ASD, P_diff = P_TD - P_ASD, L_diff = L_TD - L_ASD)
  }

  # compute observed differences
  obs <- compute_diffs(data)

  perm_diffs <- vector("list", n_perm)
  ids <- unique(data$subID)

  for (i in seq_len(n_perm)) {
    perm_labels <- data %>%
      distinct(subID, group) %>%
      mutate(group_perm = sample(group))

    perm_data <- data %>%
      left_join(perm_labels %>% select(subID, group_perm), by = "subID") %>%
      mutate(group = group_perm) %>%
      select(-group_perm)

    perm_diffs[[i]] <- compute_diffs(perm_data)
  }

  perm_diffs <- bind_rows(perm_diffs)

  pvals <- perm_diffs %>%
    summarise(
      p_A0 = mean(abs(A0_diff) >= abs(obs$A0_diff), na.rm = TRUE),
      p_P  = mean(abs(P_diff)  >= abs(obs$P_diff),  na.rm = TRUE),
      p_L  = mean(abs(L_diff)  >= abs(obs$L_diff),  na.rm = TRUE)
    )

  list(observed = obs, perm_pvals = pvals)
}

perm_results <- perm_test(hard_df, n_perm = 10000, seed = 42)
print(perm_results$observed)
print(perm_results$perm_pvals)
