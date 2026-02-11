# Load required package
#fit learning glm
library(lme4)
library(dplyr)
df <- df %>%
  mutate(
    group = case_when(
      group_verbal %in% c("ASD", "ASD (low verbal)") ~ "ASD_all",
      group_verbal == "TD" ~ "TD"
    ),
    group = factor(group)  # Make sure it's a factor
  )
# Assuming your data frame is called `df`
# Create phase variable
df <- df %>%
  mutate(
    phase = ifelse(trial <= 10, "easy", "hard"),
    trial_within_phase = ifelse(phase == "easy", trial, trial - 10),
    phase = factor(phase),  # treat as categorical
    group_verbal = factor(group)  # make sure it's a factor
  )

# Fit mixed-effects logistic regression
model <- glmer(
  score ~ phase * trial_within_phase * group + (1 | subID),
  data = df,
  family = binomial
)

# View summary
summary(model)
df$trial_scaled <- scale(df$trial_within_phase, center = TRUE, scale = TRUE)
model_scaled <- glmer(
  score ~ phase * trial_scaled * group + (1 | subID),
  data = df,
  family = binomial
)
summary(model)

hard_df <- df %>%
  filter(trial > 10) %>%
  mutate(
    group = factor(binary_group),
    trial_centered = trial - mean(trial)  # optional centering
  )


hard_df <- df %>%
  filter(trial > 11) %>%
  mutate(
    group = factor(case_when(
      group_verbal %in% c("ASD", "ASD (low verbal)") ~ "ASD_all",
      group_verbal == "TD" ~ "TD"
    )),
    trial_centered = trial - mean(trial)
  )
# Fit a basic logistic regression (no random effects)
model_simple <- glm(score ~ trial_centered * group, data = hard_df, family = binomial)

summary(model_simple)

#####approach 2#######

library(dplyr)
library(purrr)
library(tidyr)
install.packages("dplyr")  # if not already latest
library(dplyr)
# Original data prep

df = alldata[game_version == "ft" & srs_normal == 1 & above_crit == 1 & group_verbal != "flag" & group_verbal != "flag_low_td"]
df = df[,c("subID","binary_group","trial","score")]

hard_df <- df %>%
  filter(trial > 11) 
hard_df$group = hard_df$binary_group

# Function to fit learning curve
fit_learning_curve <- function(data) {
  tryCatch({
    model <- nls(
      mean_score ~ pinf - (pinf - p0) * exp(-alpha * (trial - min(trial))),
      data = data,
      start = list(pinf = 1, p0 = 0.5, alpha = 0.05),
      control = nls.control(maxiter = 10000)
    )
    as.list(coef(model))
  }, error = function(e) {
    return(list(p0 = NA, pinf = NA, alpha = NA))
  })
}
# Bootstrap function
bootstrap_curve_fit <- function(data, n_boot = 1000) {
  results <- list()
  
  for (i in 1:n_boot) {
    # Resample participants with replacement within group
    boot_data <- data %>%
      group_by(group) %>%
      group_modify(~ {
        boot_ids <- sample(unique(.x$subID), size = n_distinct(.x$subID), replace = TRUE)
        bind_rows(lapply(boot_ids, function(id) filter(.x, subID == id)))
      }) %>%
      ungroup()
    
    # Compute mean accuracy per trial
    boot_avg <- boot_data %>%
      group_by(group, trial) %>%
      dplyr::summarize(mean_score = mean(score), .groups = "drop")
    
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

# Run the bootstrap
set.seed(42)
boot_params <- bootstrap_curve_fit(hard_df, n_boot = 1000)
library(ggplot2)


ggplot(boot_params, aes(x = group, y = p0)) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.3, color = "red", fatten = 2) +
  stat_compare_means(method = "t.test", label.y = 0.5) +  # or use "wilcox.test"
  theme_minimal() +
  ylim(0,1)+
  labs(title = "Parameter a (Initial Accuracy)", x = "Group", y = "a")

ggplot(boot_params, aes(x = group, y = pinf)) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.3, color = "red", fatten = 2) +
  stat_compare_means(method = "t.test") +
  theme_minimal() +
  ylim(0,1)+
  labs(title = "Parameter b (Plateau Accuracy)", x = "Group", y = "b")

ggplot(boot_params, aes(x = group, y = alpha)) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.3, color = "red", fatten = 2) +
  stat_compare_means(method = "t.test", label.y = 0.06) +
  theme_minimal() +


  labs(title = "Parameter r (Learning Rate)", x = "Group", y = "r")

library(ggpubr)

# Compare initial accuracy (p0)
compare_p0 <- compare_means(p0 ~ group, data = boot_params, method = "t.test")

# Compare plateau accuracy (pinf)
compare_pinf <- compare_means(pinf ~ group, data = boot_params, method = "t.test")

# Compare learning rate (alpha)
compare_alpha <- compare_means(alpha ~ group, data = boot_params, method = "t.test")

# Show results
compare_p0
compare_pinf
compare_alpha


boot_summary <- boot_params %>%
  pivot_longer(cols = c(p0, pinf, alpha), names_to = "param", values_to = "value") %>%
  group_by(group, param) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    lower = quantile(value, 0.025, na.rm = TRUE),
    upper = quantile(value, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

boot_summary


boot_wide <- boot_params %>%
  pivot_wider(names_from = group, values_from = c(p0, pinf, alpha))

# Compute differences
boot_diff <- boot_wide %>%
  transmute(
    alpha_diff = alpha_TD - alpha_ASD,
    p0_diff = p0_TD - p0_ASD,
    pinf_diff = pinf_TD - pinf_ASD
  )

# Calculate two-tailed p-values
pvals <- boot_diff %>%
  summarise(
    p_alpha = 2 * min(mean(alpha_diff < 0, na.rm = TRUE), mean(alpha_diff > 0, na.rm = TRUE)),
    p_p0    = 2 * min(mean(p0_diff    < 0, na.rm = TRUE), mean(p0_diff    > 0, na.rm = TRUE)),
    p_pinf  = 2 * min(mean(pinf_diff  < 0, na.rm = TRUE), mean(pinf_diff  > 0, na.rm = TRUE))
  )

print(pvals)


boot_long <- boot_params %>%
  pivot_longer(cols = c(p0, pinf, alpha), names_to = "parameter", values_to = "value")

# Plot histograms
ggplot(boot_long, aes(x = value, fill = group)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 100) +
  facet_wrap(~parameter, scales = "free", ncol = 1) +
  labs(title = "Bootstrapped Parameter Distributions by Group",
       x = "Parameter Value", y = "Count", fill = "Group") +
  theme_minimal()


ks_results <- list(
  alpha = ks.test(
    boot_params$alpha[boot_params$group == "ASD"],
    boot_params$alpha[boot_params$group == "TD"]
  ),
  p0 = ks.test(
    boot_params$p0[boot_params$group == "ASD"],
    boot_params$p0[boot_params$group == "TD"]
  ),
  pinf = ks.test(
    boot_params$pinf[boot_params$group == "ASD"],
    boot_params$pinf[boot_params$group == "TD"]
  )
)

# Print p-values
lapply(ks_results, function(test) test$p.value)