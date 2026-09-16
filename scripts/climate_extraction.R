##############################################################################.
#
# Purpose: Download climate data from daymet to use in subsequent analyses
#
# By: Eliza
# Date: 9/16/2026
# Last modified: 9/16/2026
#
##############################################################################.

library(daymetr)
library(tidyverse)

site_coords <- data.frame(
  site = c("De", "Sg", "Bi", "Ci"),
  latitude = c(39.10004, 37.08671, 35.06076, 33.30337),
  longitude = c(-112.641, -113.562, -114.65, -114.672)
)

# Write to a temporary CSV file (required input format for batch mode)
write.csv(site_coords, "data_derived/site_coords.csv", row.names = FALSE)

# Download data for all points simultaneously
daymet_clim_raw <- download_daymet_batch(
  file_location = "data_derived/site_coords.csv",
  start = 2010,
  end = 2025,
  internal = TRUE # Returns a nested list of data frames
)

comb_df <- map_df(daymet_clim_raw, function(site_list) {
  df <- site_list$data
  df <- df %>%
    mutate(site_name = site_list$site,
           longitude = site_list$longitude,
           latitude = site_list$latitude,
           .before = 'year')
  # df$site_name <- site_list$site # Add site identifier column
  return(df)
})

# 3. Clean dates and columns for the entire dataset
daymet_clim <- comb_df %>%
  mutate(date = ymd(paste0(year, "-01-01")) - 1 + days(yday), .before = 'yday') %>%
  rename(dayl_sec = dayl..s.,
         prcp_mm = prcp..mm.day.,
         srad_wm2 = srad..W.m.2.,
         swe_kgm2 = swe..kg.m.2.,
         tmax_c = tmax..deg.c.,
         tmin_c = tmin..deg.c.,
         vp_pa = vp..Pa.) %>%
  mutate(tmean_c = (tmax_c + tmin_c)/2, .before = 'tmax_c') %>%
  mutate(dayl_hr = dayl_sec/60/60, .after = dayl_sec)

write.csv(daymet_clim, "data_derived/daymet_clim.csv", row.names = FALSE)
