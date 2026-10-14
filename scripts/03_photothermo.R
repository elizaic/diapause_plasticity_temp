###############################################################################.
# Purpose:
# photo-thermal graph analysis of diapause timing of Diorhabda
#
#
# By: Eliza
# Created: 9/30/2026
# Last modified: 9/30/2026
#
##############################################################################.



# Packages ----------------------------------------------------------------

library(tidyverse)
library(degday)

# Useful function

# create a function that will reorder a factor
reorder_populations <- function(df, col = population,
                                new_names = c("De" = "De", "Sg" = "Sg", "Bi" = "Bi", "Ci" = "Ci")) {
  df %>%
    mutate({{col}} := factor({{col}}, levels = names(new_names), labels = new_names))
}

# Climate data ---------------------------------------------------
daymet_clim <- read_csv("data_derived/daymet_clim.csv")


# Calculate degree days
tbase = 11.1
tmax = 36.7
daymet_clim <- daymet_clim %>%
  mutate(week = week(date), .after = yday) %>%
  mutate(GDD_modave = pmax(0, (ifelse(tmax_c > tmax, tmax, tmax_c) +
                                 ifelse(tmin_c < tbase, tbase, tmin_c)) / 2 - tbase),
         GDD_ss = dd_calc(tmin_c, tmax_c,
                          thresh_low = tbase, thresh_up = tmax,
                          method = 'sng_sine')
         ) %>%
  group_by(site_name, year) %>%
  mutate(
    cum_GDD_ss = cumsum(GDD_ss)
  ) %>%
  ungroup()


# minimum daylength for all populations
daymet_clim %>%
  group_by(site_name) %>%
  summarise(
    min = min(dayl_hr),
    ave_gdd = mean(max(cum_GDD_ss))
  ) %>%
  mutate(
    min_tested = 10.3, .after = min
  )

# CDL data ----------------------------------------------------------
cdl_summary <- read.csv("data_derived/new_cdl_summary.csv")

plast_summary <- read.csv("data_derived/plast_summary.csv")

plast_by_draw <- read.csv("data_derived/plast_by_draw.csv")

bean_cdl_estimates <- data.frame(
  population = c("De", 'Sg', 'Bi', "Ci"),
  cdl_2019 = c(14.44, 13.46, 11.94, 12.29)
) %>%
  reorder_populations()


# Relate CDL to field temperatures -----------------------------------------
#what is the difference between high and low temps?
daymet_clim %>% group_by(site_name) %>% summarize(
  diff = mean(tmax_c - tmin_c)
) #for each site, ranges from 15.6 to 17.0

daymet_clim %>% summarize(
  diff = mean(tmax_c - tmin_c)
) #across all sites, all years = 16

#what is the average high temperature at each site each week?
weekly_temperature <- daymet_clim %>%
  group_by(site_name, week, year) %>%
  summarise(
    week.avehigh = mean(tmax_c),
    week_cumGDD = mean(cum_GDD_ss),
    week_dayl = mean(dayl_hr)
  )

ggplot(data = weekly_temperature, aes(x = week, y = week.avehigh, color = site_name)) +
  geom_line(aes(group = interaction(as.factor(year), site_name)))

ggplot(data = weekly_temperature, aes(x = week, y = week_cumGDD, color = site_name)) +
  geom_line(aes(group = interaction(as.factor(year), site_name)))

# Temperature over the year for each site
ggplot(data = daymet_clim %>% reorder_populations(col = site_name), aes(x = yday, y = tmean_c, color = site_name)) +
  geom_line(aes(group = interaction(as.factor(year), site_name)), alpha = 0.2) +
  facet_wrap(~site_name, ncol = 4)



# Calculate weekly CDL by draw

temp_by_draw <- full_join(plast_by_draw, weekly_temperature, by = c("population" = "site_name"))

weekly_cdl_by_draw <- temp_by_draw %>%
  mutate(
    weekly_CDL = week.avehigh * cdl_plast + cdl_intercept,
    temp_out_of_range = if_else(week.avehigh > 38, "high",
                                ifelse(week.avehigh < 28, 'low', 'in_range'))
  )

week_cdl_summary <- weekly_cdl_by_draw %>%
  group_by(population, week, year) %>%
  summarise(
    week.avehigh = unique(week.avehigh),
    week_cumGDD = unique(week_cumGDD),
    week_dayl = unique(week_dayl),
    temp_out_of_range = unique(temp_out_of_range),
    median_qi(weekly_CDL)
  )
week_cdl_summary <- week_cdl_summary %>%
  reorder_populations()

# week_cdl_summary %>%
#   full_join(., weekly_cumulative_GDD, by = c('population' = 'site_name', 'week', 'year'))

week_cdl_summary_acr_yr <- weekly_cdl_by_draw %>%
  group_by(population, week) %>%
  summarise(
    mean_week.avehigh = mean(week.avehigh),
    mean_week_cumGDD = mean(week_cumGDD),
    mean_week_dayl = mean(week_dayl),
    median_qi(weekly_CDL)
  )
week_cdl_summary_acr_yr <- week_cdl_summary_acr_yr %>%
  mutate(
    temp_out_of_range = if_else(mean_week.avehigh > 38, "high",
                                ifelse(mean_week.avehigh < 28, 'low', 'in_range'))
  ) %>%
  reorder_populations()

# Fitness consequences of plasticity ------------------------------------------


# use daily daylength data
week_cdl_summary %>%
  select(population, week, year, y, ymin, ymax) %>%
  left_join(., daymet_clim, by = c('week', 'year', 'population' = 'site_name'))

# Find cross-over point, first week where CDL is less than daylength
plast_gained <- week_cdl_summary %>%
  select(population, week, year, y, ymin, ymax) %>%
  left_join(., daymet_clim, daymet_clim,
            by = c('week', 'year', 'population' = 'site_name')) %>%
  filter(week > 26) %>%
  filter(dayl_hr < y) %>%
  group_by(population, year) %>%
  slice_min(yday, n = 1, with_ties = FALSE) %>%
  select(population, yday, year, cum_GDD_ss) %>%
  rename(diap_init_yday_plastic = yday,
         diap_init_GDD_plastic = cum_GDD_ss)

# when does CDL happen each year for constant CDL
const_gained <- week_cdl_summary %>%
  select(population, week, year, y, ymin, ymax) %>%
  left_join(., daymet_clim, daymet_clim,
            by = c('week', 'year', 'population' = 'site_name')) %>%
  left_join(., bean_cdl_estimates) %>%
  filter(week > 26) %>%
  filter(dayl_hr < cdl_2019) %>%
  group_by(population, year) %>%
  slice_min(yday, n = 1, with_ties = FALSE) %>%
  select(population, yday, year, cum_GDD_ss) %>%
  rename(diap_init_yday_constant = yday,
         diap_init_GDD_constant = cum_GDD_ss)

first_frost <- daymet_clim %>%
  filter(week > 26) %>%
  filter(tmin_c <= 0) %>%
  group_by(site_name, year) %>%
  slice_min(yday, n = 1, with_ties = FALSE) %>%
  select(site_name, year, yday) %>%
  rename(first_frost = yday) %>%
  ungroup() %>%
  complete(site_name, year, fill = list(first_frost = 365))

gained_yearly <- full_join(plast_gained, const_gained, by = c('population', 'year')) %>%
  full_join(., first_frost, by = c('year', 'population' = 'site_name')) %>%
  # find difference between plastic and constant diapause initiation days and summarise
  mutate(
    diff_yday = diap_init_yday_plastic - diap_init_yday_constant,
    diff_GDD = diap_init_GDD_plastic - diap_init_GDD_constant,
    frost_plast = first_frost- diap_init_yday_plastic,
    frost_constant = first_frost - diap_init_yday_constant
  ) %>%
  reorder_populations()


# summary
gained_summary <- gained_yearly %>%
  group_by(population) %>%
  summarise(
    across(c(diff_yday, diff_GDD, frost_plast, frost_constant),
           .fns =
             list(mean = mean, min = min, max = max),
           .names = "{.col}_{.fn}"),
  )

gained_yearly %>%
  ggplot(aes(x = year, y = diff_GDD, color = population)) +
  geom_line() +
  labs(y = "Degree Days Gained") +
  theme_classic()


gained_summary %>%
  pivot_longer(cols = starts_with('frost') & ends_with('mean'),
               names_to = 'metric',
               values_to = 'days_before_frost') %>%
  ggplot() +
  geom_point(aes(x = population, y = days_before_frost, color = metric)) +
  theme_classic()


gained_yearly %>%
  ggplot() +
  geom_jitter(aes(x = population, y = diff_GDD, color = population),
              width = 0.2, alpha = 0.5, shape = 16) +
  geom_point(data = gained_summary, aes(x = population, y = diff_GDD_mean,
                                        color = population),
             size = 4, shape = 16) +
  geom_hline(yintercept = 0, color = 'grey50', linetype = 2) +
  theme_classic()


gained_yearly %>%
  ggplot() +
  geom_jitter(aes(x = population, y = diap_init_yday_plastic, color = population),
              width = 0.2, alpha = 0.5, shape = 16) +
  # geom_point(data = gained_summary, aes(x = population, y = diff_GDD_mean,
  #                                       color = population),
  # size = 4, shape = 16) +
  geom_jitter(aes(x = population, y = first_frost, color = population),
              width = 0.2, color = 'cornflowerblue') +
  scale_y_date(date_breaks = '1 month', date_labels = '%b') +
  theme_classic()




# Photothermographs -----------------------------------------------------------

const_cross_point <- week_cdl_summary %>%
  select(population, week, year, y, ymin, ymax) %>%
  left_join(., daymet_clim, daymet_clim,
            by = c('week', 'year', 'population' = 'site_name')) %>%
  left_join(., bean_cdl_estimates) %>%
  filter(week > 26) %>%
  filter(dayl_hr < cdl_2019) %>%
  group_by(population, year) %>%
  slice_min(yday, n = 1, with_ties = FALSE) %>%
  select(population, yday, year, cum_GDD_ss, dayl_hr) %>%
  group_by(population) %>%
  summarise(
    dayl_hr = unique(dayl_hr),
    mean_gdd_constant = mean (cum_GDD_ss)) %>%
  reorder_populations()

plast_cross_point <- week_cdl_summary %>%
  select(population, week, year, y, ymin, ymax) %>%
  left_join(., daymet_clim, daymet_clim,
            by = c('week', 'year', 'population' = 'site_name')) %>%
  filter(week > 26) %>%
  filter(dayl_hr < y) %>%
  group_by(population, year) %>%
  slice_min(yday, n = 1, with_ties = FALSE) %>%
  select(population, yday, year, cum_GDD_ss, dayl_hr) %>%
  group_by(population) %>%
  summarise(
    dayl_hr = mean(dayl_hr),
    mean_gdd_plast = mean (cum_GDD_ss)) %>%
  reorder_populations()


# Photothermograph for each year
ggplot(data = week_cdl_summary, aes(x = week_cumGDD, y = y, group = year)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = 'steelblue3', alpha = 0.2) +
  geom_line(aes(group = year, color = temp_out_of_range), alpha = 0.5) +
  geom_line(aes(y = week_dayl, color = 'zdaylength')) +
  facet_wrap( ~ population, nrow = 1) +
  theme_classic(base_size = 16)


# Photothermograph for all years
week_cdl_summary_acr_yr %>%
  ggplot(aes(x = mean_week_cumGDD, y = y)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = 'steelblue3', alpha = 0.2) +
  geom_line(aes(group = population, color = 'Plastic CDL')) +
  # add daylength line
  geom_line(aes(y = mean_week_dayl, color = 'Daylength')) +
  # cross point for constant cdl
  geom_point(data = const_cross_point,
             aes(x = mean_gdd_constant, y = dayl_hr, shape = 'without plasticity')) +
  # geom_hline(data = bean_cdl_estimates, aes(yintercept = cdl_2019),
  #            linetype = 2 ) +
  # cross point for plastic cdl
  geom_point(data = plast_cross_point,
             aes(x = mean_gdd_plast, y = dayl_hr, shape = 'with plasticity')) +
  geom_vline(aes(xintercept = 627), linetype = 3) +
  facet_wrap( ~ population, nrow = 1, scales = "free_x") +
  scale_shape_manual(values = c('without plasticity' = 1, 'with plasticity' = 16)) +
  scale_color_manual(values = c('Plastic CDL' = 'steelblue3', 'Daylength' = 'goldenrod1')) +
  labs(y = "Daylength (hr)",
       x = "Cumulative GDD",
       shape = "Diapause initiation",
       color = "") +
  theme_classic(base_size = 16)




