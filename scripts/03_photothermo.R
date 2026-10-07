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


# Load data -------------------------------------------------------


## Climate data ---------------------------------------------------
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
                          method = 'sng_sine'))


# minimum daylength for all populations
daymet_clim %>%
  group_by(site_name) %>%
  summarise(
    min = min(dayl_hr)
  ) %>%
  mutate(
    min_tested = 10.3
  )

# CDL data ----------------------------------------------------------
cdl_summary <- read.csv("data_derived/new_cdl_summary.csv")


# Relate CDL to field temperatures -----------------------------------------
#what is the difference between high and low temps?
daymet_clim %>% group_by(site_name) %>% summarize(
  diff = mean(tmax_c - tmin_c)
) #for each site, ranges from 15.6 to 17.0

daymet_clim %>% summarize(
  diff = mean(tmax_c - tmin_c)
) #across all sites, all years = 16

#what is the average high temperature at each site each week?
plast_CDL <- daymet_clim %>%
  group_by(site_name, week, year) %>%
  summarise(
    week.avehigh = mean(tmax_c)
  ) %>%
  mutate(
    pop.index = case_when (
      site_name == "De" ~ 1,
      site_name == "Sg" ~ 2,
      site_name == "Bi" ~ 3,
      site_name == "Ci" ~ 4,
    ), .after = "site_name")

ggplot(data = plast_CDL, aes(x = week, y = week.avehigh, color = site_name)) +
  geom_line(aes(group = interaction(as.factor(year), site_name)))

# Temperature over the year for each site
ggplot(data = daymet_clim, aes(x = yday, y = tmean_c, color = site_name)) +
  geom_line(aes(group = interaction(as.factor(year), site_name)), alpha = 0.2) +
  facet_wrap(~site_name, ncol = 4)


# # CDL estimates
# str(cdl_summary)
# # CDL_estimates$population <- factor(CDL_estimates$population, levels = c("De", "Sg", "Bi", "Ci"), labels = c("Delta", "StGeorge", "BigBend", "Cibola") )
# # CDL_estimates$temperature <- factor(CDL_estimates$temperature, levels = c("38", "28"))
#
# # CDL_estimates <- CDL_estimates %>% mutate(temperature.cont = as.double(as.character(temperature)))
# # str(CDL_estimates)
#
# slopes_by_group <- cdl_summary %>%
#   group_by(population) %>%
#   group_modify(~ broom::tidy(lm(cdl_reliable ~ as.double(temperature), data = .x))) %>%
#   filter(term == "as.double(temperature)")
#
# lm(cdl_reliable ~ as.double(temperature), data = cdl_summary) %>%
#   emmeans::emtrends(~population)
#
# plast_slope_delta <- lm(daylength_hr ~ temperature.cont, data = CDL_estimates %>% filter(population == 'Delta'))
# plast_slope_StGeorge <- lm(daylength_hr ~ temperature.cont, data = CDL_estimates %>% filter(population == 'StGeorge'))
# plast_slope_BigBend <- lm(daylength_hr ~ temperature.cont, data = CDL_estimates %>% filter(population == 'BigBend'))
# plast_slope_Cibola <- lm(daylength_hr ~ temperature.cont, data = CDL_estimates %>% filter(population == 'Cibola'))
#
# plast.estimates <- rbind(
#   coef(plast_slope_delta) %>% as.array(),
#   coef(plast_slope_StGeorge) %>% as.array(),
#   coef(plast_slope_BigBend) %>% as.array(),
#   coef(plast_slope_Cibola) %>% as.array()
# ) %>% as.data.frame() %>% mutate(population = c("Delta", "StGeorge", "BigBend", "Cibola"), .before = "(Intercept)") %>%
#   mutate(pop.index = case_when (
#     population == "Delta" ~ 1,
#     population == "StGeorge" ~ 2,
#     population == "BigBend" ~ 3,
#     population == "Cibola" ~ 4,
#   ), .before = "(Intercept)")
#
# ggplot(data = CDL_estimates, aes(x = temperature.cont, y = daylength_hr, color = population, linetype = factor(year))) +
#   geom_line() +
#   lims(x = c(0,38), y = c(0,22))

