##############################################################################.
# Purpose:
# CDL analysis using Bayesian logistic regression in brms
#
#
# By: Eliza
# Created: 9/9/2026
# Last modified: 9/9/2026
#
##############################################################################.




# Packages ----------------------------------------------------------------

library(tidyverse)
library(brms)
library(readxl)
library(ggpubr)
library(tidybayes)
library(bayesplot)


select <- dplyr::select

# Load Data ---------------------------------------------------------------

# 2020
data_2020_raw <- read_excel("plasticity_experiment_data/2020diapause_datasheet_final_data.xlsx", col_types = c("guess","guess","guess","guess","guess", "date","date","date","guess","date", "guess", "guess"))

# 2022
data_2022_raw <- read_excel("plasticity_experiment_data/2022_temp_photoperiod_data.xlsx")


# pop_labels <- c("De" = "Delta\n(39°N)", "Sg" = "St. George\n(37°N)", "Bi" = "Big Bend\n(35°N)", "Ci" = "Cibola\n(33°N)")
# temp_labels = c("38" = "38°", "28" = "28°")


# Clean Data -----------------------------------------------------------------

# turn 2022 data into long, uncounted dataframe, instead of grouped

data_2022_long <- data_2022_raw %>%
  select(-no_pairs, -prop_repro) %>%
  pivot_longer(
    cols = c(no_reproductive, no_diapause),
    names_to = "outcome",
    values_to = 'count'
  ) %>%
  uncount(count) %>%
  mutate(eggs_present = if_else(outcome == "no_reproductive", 1, 0),
         year = 2022) %>%
  select(year, population, daylength, daylength_hr, daylength_block, temperature, eggs_present)

# prep 2020 data
data_2020_raw1 <- data_2020_raw %>%
  mutate(year = 2020,
         daylength_hr = if_else(daylength == 1030, 10.5,
                                if_else(daylength==1125, 11.42,
                                        if_else(daylength==1220, 12.33,
                                                if_else(daylength==1315, 13.25,
                                                        if_else(daylength==1410, 14.17,
                                                                if_else(daylength==1505, 15.08,111111111)))))),
         daylength_block = 2020
  ) %>%
  select(year, population, daylength, daylength_hr, daylength_block, temperature, eggs_present)

# Make sure these match before combining datasets
names(data_2022_long)
names(data_2020_raw1)


# combine 2020 and 2022 datasets

all_data <- rbind(data_2020_raw1, data_2022_long)

# format columns how I want them

str(all_data)
all_data <- all_data %>%
  mutate(temperature = if_else(temperature == "28/13", '28', if_else(
    temperature == "38/23", '38', temperature))) %>%
  filter(!is.na(eggs_present)) %>%
  mutate(
    # center daylength for interpretability and better sampling geometry;
    # save the mean to back-transform later
    daylength_mean = mean(daylength_hr, na.rm = TRUE),
    daylength_centered = daylength_hr - daylength_mean
  )

all_data$year <- factor(all_data$year)
all_data$population <- factor(all_data$population, levels = c("De", "Sg", "Bi", "Ci"))
all_data$daylength <- factor(all_data$daylength)
all_data$daylength_block <- factor(all_data$daylength_block)
all_data$temperature <- factor(all_data$temperature, levels = c("38", "28"))

str(all_data)

all_data %>%
  count(population, temperature, daylength_block) %>%
  print(n = 100)

# summary/grouped data - with proportion reproductive
grouped_summary <- all_data %>%
  group_by(daylength_hr, temperature, population) %>%
  summarise(n = n(),
            eggs_present_sum = sum(eggs_present)) %>%
  mutate(prop = eggs_present_sum / n)


# Sample sizes ---------------------------------------------------------------
all_data %>%
  group_by(population, temperature, daylength) %>%
  summarise(n = n())

all_data %>% distinct(daylength) %>% nrow()

all_data %>% distinct(daylength_hr) %>% range()


# Logistic Regression with BRMS -----------------------------------------------


# Priors
# Weakly informative priors on the logit scale. Adjust sd if you have
# substantive prior knowledge (e.g., from prior published critical
# photoperiod studies) about how strong these effects should be.

bpriors <- c(
  prior(normal(0, 1.5), class = "b"),          # slopes / interactions
  prior(normal(0, 2),   class = "Intercept")#,  # baseline
  # prior(exponential(1), class = "sd")          # batch random-effect SD
)

# Fit model
fit1 <- brm(eggs_present ~ temperature * daylength_centered * population +
              daylength_block,
            data = all_data,
            family = bernoulli(link = "logit"),
            prior = bpriors,
            iter = 5000, warmup = 2000,
            chains = 4, cores = 4,
            seed = 123,
            control = list(adapt_delta = 0.95),
            file = "models/brms_fit1",
            file_refit = 'on_change')

summary(fit1)

plot(fit1)
plot(conditional_effects(fit1, effects = "daylength_centered:population"),
     conditions = 'temperature')

# ggarrange(plotlist = list(plot(conditional_effects(fit1))))

# Model checking -------------------------------------------------------------



## Other checks --------------------------------------------------------------
pp_check(fit1, ndraws = 100)                       # posterior predictive check
pp_check(fit1, type = "error_binned")              # binned residuals vs predictors

# Check whether the linear-in-logit assumption for daylength is reasonable
# by comparing to a smooth term per population:
# fit_smooth <- brm(
#   diapause ~ s(daylength_centered, by = population, k = 5) * temperature + (1 | daylength_block),
#   data = df, family = bernoulli(), chains = 4, cores = 4,
#   iter = 3000, warmup = 1000, control = list(adapt_delta = 0.95), seed = 123
# )
# # Compare fit via approximate leave-one-out CV
# fit        <- add_criterion(fit, "loo")
# fit_smooth <- add_criterion(fit_smooth, "loo")
# loo_compare(fit, fit_smooth)
# If the smooth model doesn't fit meaningfully better, the linear logistic
# model is an adequate description of the daylength response.



# Population differences --------------------------------------------------

# Conditional effects plot: predicted diapause probability vs daylength,
#     separately by population and temperature
conditional_effects(
  fit1,
  effects = "daylength_centered:population",
  conditions = data.frame(temperature = levels(all_data$temperature))
)


# (b) Estimate critical photoperiod (CPP = daylength at 50% diapause) for
#     each population x temperature combination, with full posterior
#     uncertainty. CDL is where logit(p) = 0, i.e.
#     Intercept + pop_effect + temp_effect + ... = -(slope terms) * daylength_c = 0

draws <- as_draws_df(fit1)

# Build a grid over which to compute predicted probabilities, then
# back-solve for CDL numerically (robust to the exact parameterization)
grid <- expand_grid(
  daylength_centered = seq(min(all_data$daylength_centered), max(all_data$daylength_centered), length.out = 200),
  population  = levels(all_data$population),
  temperature = levels(all_data$temperature),
  daylength_block = levels(all_data$daylength_block)
)

preds_by_block <- fitted(
  fit1, newdata = grid,
  summary = FALSE                          # averaging over batch effect
)  # returns a draws x rows matrix of predicted probabilities


# average across blocks, per draw, per daylength_centered x population x
# temperature combination -- this gives the block-averaged prediction
# equivalent to what re_formula = NA gave with the random-effect model
grid_key <- grid %>%
  mutate(row = row_number()) %>%
  group_by(daylength_centered, population, temperature) %>%
  summarize(rows = list(row), .groups = "drop")

n_draws <- nrow(preds_by_block)
preds <- matrix(
  NA_real_,
  nrow = n_draws,
  ncol = nrow(grid_key)
)
for (j in seq_len(nrow(grid_key))) {
  cols <- grid_key$rows[[j]]
  preds[, j] <- rowMeans(preds_by_block[, cols, drop = FALSE])
}
# `preds` is now a draws x (daylength_centered x population x temperature)
# matrix of block-averaged predicted probabilities, matching the shape the
# rest of the CPP-finding / boundary-probability code expects. Rebuild
# `grid` to match preds' column ordering (drop daylength_block, now that
# it's been averaged over):
grid <- grid_key %>% select(-rows)




# CDL Calculation -----------------------------------------------------------

# For each posterior draw and each population x temperature combination,
# find the daylength_c at which predicted probability crosses 0.5
# IMPORTANT: when a population's diapause probability never gets close to
# 0.5 across the tested daylength range, the crossing point is not
# identified by the data -- the slope for that group is weakly estimated,
# and only a biased subset of posterior draws (the unusually steep ones)
# will happen to cross 0.5 within the grid. Averaging with na.rm = TRUE
# over just those draws silently conditions on that subset and produces a
# misleading number. So instead of just recording cpp_c, we also record
# WHETHER each draw crossed at all, and in which direction it failed to
# cross if not -- this lets you see and report the non-identification
# rather than papering over it.

cdl_by_draw <- grid %>%
  mutate(row = row_number()) %>%
  group_by(population, temperature) %>%
  group_modify(function(g, key) {
    n_draws <- nrow(preds)
    results <- map(seq_len(n_draws), function(d) {
      p <- preds[d, g$row]
      idx <- which(diff(sign(p - 0.5)) != 0)[1]
      if (is.na(idx)) {
        # curve never crosses 0.5 in the tested range -- record which side
        status <- if (mean(p) < 0.5) "below_range" else "above_range"
        return(tibble(cdl_c = NA_real_, status = status))
      }
      x0 <- g$daylength_centered[idx]; x1 <- g$daylength_centered[idx + 1]
      p0 <- p[idx]; p1 <- p[idx + 1]
      cdl_val <- x0 + (0.5 - p0) * (x1 - x0) / (p1 - p0)
      tibble(cdl_c = cdl_val, status = "crosses")
    })
    bind_rows(results) %>% mutate(draw = seq_len(n_draws))
  }) %>%
  ungroup() %>%
  mutate(cdl = cdl_c + unique(all_data$daylength_mean))

# Diagnostic: what fraction of posterior draws actually cross 0.5 within
# your tested daylength range, per population x temperature? Low values
# here (e.g. well under ~90-95%) mean the CPP is poorly identified for
# that group and any point estimate/interval below should be treated with
# real caution -- or not reported at all.
crossing_diagnostic <- cdl_by_draw %>%
  count(temperature, population, status) %>%
  group_by(population, temperature) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()
crossing_diagnostic

# Only summarize CPP for draws that actually cross -- and report the
# crossing fraction alongside so the reader knows how much of the
# posterior this interval represents. For groups where "crosses" is a
# small fraction, treat the interval as unreliable/non-identified rather
# than reporting it as a normal credible interval.
cdl_summary <- cdl_by_draw %>%
  group_by(population, temperature) %>%
  summarize(
    frac_crosses = mean(status == "crosses"),
    cdl_mean  = mean(cdl[status == "crosses"], na.rm = TRUE),
    cdl_lower = quantile(cdl[status == "crosses"], 0.025, na.rm = TRUE),
    cdl_upper = quantile(cdl[status == "crosses"], 0.975, na.rm = TRUE),
    .groups = "drop"
  )
cdl_summary <- cdl_summary %>%
  mutate(
    reliable = if_else(frac_crosses > 0.95, 1, 0),
    cdl_reliable = if_else(frac_crosses > 0.95, cdl_mean, 10.3)
  )
cdl_summary
cdl_summary$population <- factor(cdl_summary$population, levels = c("De", "Sg", "Bi", "Ci"))


ggplot(data = cdl_summary, aes(x = temperature, y = cdl_reliable, color = population,
                               group = population)) +
  geom_point(aes(shape = as.factor(reliable), size = as.factor(reliable)),
             position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = cdl_lower, ymax = cdl_upper, linetype = as.factor(reliable)),
                width = 0,
                position = position_dodge(width = 0.5)) +
  geom_line(position = position_dodge(width = 0.5)) +
  labs(x = 'Temperature',
       y = 'Critical daylength (hrs)',
       color = "Site") +
  scale_shape_manual(values = c(6, 16)) +
  scale_linetype_manual(values = c(0, 1)) +
  scale_size_manual(values = c(3.5, 2.2)) +
  scale_color_viridis_d() +
  guides(shape = 'none', linetype = 'none', size = 'none') +
  theme_classic(base_size = 16)

# # 6b. Robust fallback: boundary probabilities instead of forcing a CPP

# # For groups where CPP is poorly identified (low frac_crosses above), a more
# # honest and fully-identified summary is the predicted diapause probability
# # AT your tested daylength extremes, with credible intervals -- this is
# # well-estimated by the data regardless of whether the curve crosses 0.5.
#
boundary_grid <- expand_grid(
  daylength_centered = c(min(all_data$daylength_centered), max(all_data$daylength_centered)),
  population  = levels(all_data$population),
  temperature = levels(all_data$temperature),
  daylength_block = levels(all_data$daylength_block)
) %>%
  mutate(daylength = daylength_centered + unique(all_data$daylength_mean))

boundary_preds <- fitted(fit1, newdata = boundary_grid)
boundary_summary <- bind_cols(boundary_grid, as_tibble(boundary_preds))
boundary_summary



# Direct posterior comparison of CPP between populations (within a
# temperature), e.g. population A vs B at the low temperature:
cdl_wide <- cdl_by_draw %>%
  filter(temperature == levels(all_data$temperature)[1]) %>%
  select(draw, population, cdl) %>%
  pivot_wider(names_from = population, values_from = cdl)

# Replace "PopA","PopB" with your actual population labels
# diff_draws <- cpp_wide$PopA - cpp_wide$PopB
# mean(diff_draws); quantile(diff_draws, c(0.025, 0.975))
# mean(diff_draws > 0)   # posterior probability PopA has higher CPP than PopB
#####



grid2 <- expand_grid(
  daylength_centered = seq(min(all_data$daylength_centered), max(all_data$daylength_centered), length.out = 100),
  population  = levels(all_data$population),
  temperature = levels(all_data$temperature),
  daylength_block = levels(all_data$daylength_block)
) %>%
  add_epred_draws(fit1) %>%
  ungroup() %>%
  group_by(daylength_centered, population, temperature) %>%
  summarise(
    value = median(.epred),
    .lower = quantile(.epred, probs = 0.025),
    .upper = quantile(.epred, probs = 0.975),
    .groups = "drop"
  )
#

grid2 <- grid2 %>%
  mutate(daylength = daylength_centered + unique(all_data$daylength_mean))
grid2$population <- factor(grid2$population, levels = c("De", "Sg", "Bi", "Ci"))
grid2$temperature <- factor(grid2$temperature, levels = c("38", "28"))


ggplot(grid2, aes(x = daylength, y = value, color = population)) +
  geom_ribbon(aes(x = daylength, ymin = .lower, ymax = .upper, fill = population),
              alpha = 0.2, color = NA) +
  geom_line() +
  facet_wrap(~ factor(temperature, levels = c('28', '38'))) +
  # add cdl
  # geom_rect(data = cdl_summary %>%
  #             filter(frac_crosses > 0.95),
  #              aes(xmin = cdl_lower, xmax = cdl_upper, ymin = -Inf, ymax = Inf, fill = population),
  #             inherit.aes = FALSE, alpha = 0.2) +
  # geom_vline(data = cdl_summary %>% filter(frac_crosses > 0.95),
  #            aes(xintercept = cdl_mean, color = population)) +
  # geom_hline(yintercept = 0.5, linetype = 2) +
  geom_point(data = grouped_summary, aes(x = daylength_hr, y = prop)) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  labs(x = 'Daylength (hrs)',
       y = 'Proportion reproductive',
       color = "Site",
       fill = "Site") +
  theme_classic(base_size = 16)


ggplot(grid2, aes(x = daylength, y = value, color = temperature)) +
  geom_ribbon(aes(x = daylength, ymin = .lower, ymax = .upper, fill = temperature),
              alpha = 0.2, color = NA) +
  geom_line() +
  facet_wrap(~ population, ncol = 1) +
  # add cdl
  geom_rect(data = cdl_summary %>%
              filter(frac_crosses > 0.95),
               aes(xmin = cdl_lower, xmax = cdl_upper, ymin = -Inf, ymax = Inf, fill = temperature),
              inherit.aes = FALSE, alpha = 0.2) +
  geom_vline(data = cdl_summary %>% filter(frac_crosses > 0.95),
             aes(xintercept = cdl_mean, color = temperature)) +
  # geom_hline(yintercept = 0.5, linetype = 2) +
  geom_point(data = grouped_summary, aes(x = daylength_hr, y = prop)) +
  labs(x = 'Daylength (hrs)',
       y = 'Proportion reproductive',
       color = "Temperature",
       fill = "Temperature") +
  scale_color_manual(values = c("#E69F00", "#0072B2"), aesthetics = c('color', 'fill')) +
  theme_classic(base_size = 16)




# ggplot() +
# geom_point(data = all_data %>% filter(eggs_present == 1),
#            aes(y = 1.2, x = daylength_hr, color = population),
#            position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0.005,
#                                            dodge.width = 1200)) +
#   # geom_boxplot(data = all_data %>% filter(eggs_present == 1), aes(y = 1, x = daylength_hr,
#   #                                   color = population, fill = population), alpha = 0.2,
#   #              position = position_dodge(width = 0.8)) +
#   # geom_boxplot(data = all_data %>% filter(eggs_present == 0), aes(y = 0, x = daylength_hr,
#   #                                                                 color = population, fill = population), alpha = 0.2,
#   #              position = position_dodge(width = 0.8)) +
#   facet_wrap(~ temperature) +
#   lims(y = c(0, 1.5)) +
#   theme_classic()


# Hypothesis tests on interaction terms directly -----------------------------

# Quick look at whether the daylength:population and 3-way terms are
# credibly different from zero (check which coefficient names your fit uses):
fixef(fit1)

# e.g. hypothesis() lets you test specific combinations, such as whether
# the daylength slope differs between two named populations:
hypothesis(fit1, "daylength_centered:population = daylength_centered:populationCi")


emmeans::emmeans(fit1, pairwise ~ daylength_block|temperature|population, type = 'response')


