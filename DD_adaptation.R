library(lubridate)
library(geosphere)
library(degday)
library(RColorBrewer)
library(lme4)
library(lmerTest)


month.letters = c('J','F','M','A','M','J','J','A','S','O','N','D')


# Import USPEST data ----
# create list of file names for downloaded USPEST csv files
list_of_files <- list.files(path = "Deg_days_USPESTmodel",
                            recursive = FALSE,
                            pattern = "\\.csv$",
                            full.names = TRUE)

# import all the files, add file name as a column
uspest <- read_csv(list_of_files, id = "file_name", skip = 1, show_col_types = F, name_repair = "universal")

# add data columns and reformat
uspest <- uspest %>%
  separate(col = file_name, sep = c(42, -12, -8), into = c(NA, 'station', 'year', NA)) %>% # separate out file name to get year and station
  mutate(population = case_when(station == "CBRA3_CIBOLA_AZ_" ~ "Cibola", # add population column, from station name
                          station == "CMP08_Delta_UT_" ~ "Delta",
                          station == "DLTU1_USRCRN_SITE_AT_DELTA_M_UT_" ~ "Delta",
                          station == "KIFP_Bullhead_Cty_LaughlinB_AZ_" ~ "BigBend",
                          station == "KSGU_St_George_Muni_Apt_UT_" ~ "StGeorge"
            ), .before = "year") %>%
  mutate(latitude = case_when(population == "Delta" ~ 39.100040, # add column for latitude, for daylength calculation
                              population == "StGeorge" ~ 37.086714,
                              population == "BigBend" ~ 35.060760,
                              population == "Cibola" ~ 33.303369)
         ) %>%
  unite(col = 'date', year:day, sep = '-', remove = F) %>% # create date column from year, month, day columns
  mutate(daylength = daylength(lat = latitude, doy = date)) # calculate daylength with geosphere package fn
uspest$date <- ymd(uspest$date)
uspest$year <- as.factor(uspest$year)
uspest$week <- week(uspest$date)
uspest$population <- factor(uspest$population, levels = c("Delta", 'StGeorge', 'BigBend', "Cibola"))
str(uspest)

# Relate CDL to field temperatures ----
#what is the difference between high and low temps?
uspest %>% group_by(population) %>% summarize(
  diff = mean(maximum.temperature.C - minimum.temperature.C)
) #for each site, ranges from 12.5 to 19.2

uspest %>% summarize(
  diff = mean(maximum.temperature.C - minimum.temperature.C)
) #across all sites, all years = 15.2

#what is the average high temperature at each site each week?
plast_CDL <- uspest %>% group_by(population, week, year) %>% summarise(
  week.avehigh = mean(maximum.temperature.C)
) %>% mutate(pop.index = case_when (
  population == "Delta" ~ 1,
  population == "StGeorge" ~ 2,
  population == "BigBend" ~ 3,
  population == "Cibola" ~ 4,
), .after = "population")

ggplot(data = plast_CDL, aes(x = week, y = week.avehigh, color = population, linetype = year)) +
  geom_line()

#extract slope and intercept from each population
###!!!!!! CDL_estimates file uploaded in 2022... file----
CDL_estimates <- CDL_estimates %>% mutate(temperature.cont = as.double(as.character(temperature)))
str(CDL_estimates)

plast_slope_delta <- lm(daylength_hr ~ temperature.cont, data = CDL_estimates %>% filter(population == 'Delta'))
plast_slope_StGeorge <- lm(daylength_hr ~ temperature.cont, data = CDL_estimates %>% filter(population == 'StGeorge'))
plast_slope_BigBend <- lm(daylength_hr ~ temperature.cont, data = CDL_estimates %>% filter(population == 'BigBend'))
plast_slope_Cibola <- lm(daylength_hr ~ temperature.cont, data = CDL_estimates %>% filter(population == 'Cibola'))

plast.estimates <- rbind(
  coef(plast_slope_delta) %>% as.array(),
  coef(plast_slope_StGeorge) %>% as.array(),
  coef(plast_slope_BigBend) %>% as.array(),
  coef(plast_slope_Cibola) %>% as.array()
) %>% as.data.frame() %>% mutate(population = c("Delta", "StGeorge", "BigBend", "Cibola"), .before = "(Intercept)") %>%
  mutate(pop.index = case_when (
    population == "Delta" ~ 1,
    population == "StGeorge" ~ 2,
    population == "BigBend" ~ 3,
    population == "Cibola" ~ 4,
  ), .before = "(Intercept)")

ggplot(data = CDL_estimates, aes(x = temperature.cont, y = daylength_hr, color = population, linetype = factor(year))) +
  geom_line() +
  lims(x = c(0,38), y = c(0,22))

#write a function for a line that represents how CDL changes with temperature

plast_CDL <- plast_CDL %>% mutate(week.CDL = 0)

for (i in 1:nrow(plast_CDL)) {
  plast_CDL[i,6] <- as.numeric(plast_CDL[i,5]) * plast.estimates[as.numeric(plast_CDL[i,2]),4] + plast.estimates[as.numeric(plast_CDL[i,2]),3]
}
plast_CDL$week.CDL <- as.numeric(plast_CDL$week.CDL)
str(plast_CDL)

#Plot the weekly CDLs
ggplot(data = plast_CDL, aes(x = week, y = week.CDL, color = population, linetype = year)) +
  geom_line()

#relate weekly CDL to cumulative degree days
uspest <- uspest %>% mutate(dayofyear = yday(date), .after = "day")
ave.cum.dd <- uspest %>% dplyr::select(population, date, year, week, dayofyear, cumulative.degree.days.C)

#   group_by(population, dayofyear) %>% summarise(
#   ave.cum.degdays = mean(cumulative.degree.days.C)
# )
# ave.cum.dd <- ave.cum.dd %>% mutate(date = as.Date(dayofyear-1, origin = "2020-01-01")) %>%
#   mutate(week = week(date))
week.cum.dd <- ave.cum.dd %>% group_by(population, week, year) %>% summarise(
  week.cumDD = max(cumulative.degree.days.C)
)

plast_CDL <- left_join(x = plast_CDL, y = week.cum.dd)

#make table of constant CDLs in 2019 from Bean et al. manuscript in prep Voltinism Range Expansion 2023
bean_cdl_estimates <- data.frame(
  population = c("Delta", 'StGeorge', 'BigBend', "Cibola"),
  cdl.2019 = c(14.44, 13.46, 11.94, 12.29)
  )
bean_cdl_estimates$population <- factor(bean_cdl_estimates$population, levels = c("Delta", 'StGeorge', 'BigBend', "Cibola"))

#Make data frame of points for the first day of the month
first_ofthe_month_all <- uspest %>% filter(day == 1) %>% filter(year == 2022)

#Make annotation data frame
G1.annotation <- data.frame(population = "Delta", lab = "G1 adult emergence")
G1.annotation$population <- factor(G1.annotation$population, levels = c("Delta", 'StGeorge', 'BigBend', "Cibola"))

#Make population labels
pop.names <- c("Delta" = "Delta (39°N)",
               "StGeorge" = "St. George (37°N)",
               "BigBend" = "Big Bend (35°N)",
               "Cibola" = "Cibola (33°N)")


# Photothermographs ----

## Sites separately - facet plot, with constant CDL

ggplot() +
  geom_line(data = uspest, aes(x = cumulative.degree.days.C, y = daylength, color = year)) +
  # labs(title = "Delta") +
  geom_hline(data = CDL_estimates, aes(yintercept = daylength_hr, linetype = temperature), size = 1) +
  facet_wrap(~ population, nrow = 4) +
  theme_bw(base_size = 15)

## Sites separately - facet plot, with plastic CDL
ggplot() +
  #overwintering emergence vertical line
  geom_vline(xintercept = 627, color = "grey 50", size = 0.9) +
  geom_text(data = G1.annotation, aes(x = 655, y = 6.75, group = population),
            label = "G1 adult emergence", hjust = "left") +

  #Deg days & daylength blue lines
  geom_line(data = uspest, aes(x = cumulative.degree.days.C, y = daylength, color = year)) +
  scale_color_manual(values = colorRampPalette(brewer.pal(8, "Blues"))(11) )+

  #plastic CDL lines
  geom_line(data = plast_CDL, aes(x = week.cumDD, y = week.CDL, alpha = year), size = .7) +
  geom_text(data = plast_CDL %>% filter(week == 53 & population == "Delta" & year == 2020), aes(x = week.cumDD, y = week.CDL, group = population),
            label = "Plastic CDL", hjust = "left", nudge_y = 0.7, nudge_x = 100) +

  #constant CDL lines
  geom_hline(data = bean_cdl_estimates, aes(yintercept = cdl.2019), linetype = 2, size = 1) +
  geom_text(data = bean_cdl_estimates %>% filter(population == "Delta"), aes(x = 4500, y = cdl.2019, group = population),
            label = "Non-plastic CDL", hjust = "center", nudge_y = -0.5) +

  #Month labels
  geom_point(data = first_ofthe_month_all, aes(x = cumulative.degree.days.C, y = daylength)) +
  geom_text(data = first_ofthe_month_all, aes(x = cumulative.degree.days.C, y = daylength),
            label = rep(month.letters,4), nudge_y = -0.5, nudge_x = 60) +

  #DD gained labels
  # geom_text(data = DD.gained.results, aes(x = 820, y = 16.75, group = population),
  #           label = "Degree days gained with plasticity:", hjust = "left", size = 4.5) +
  # geom_text(data = DD.gained.results, aes(x = 3700, y = 16.75, group = population),
  #           label = DD.gained.results$DD.gained, hjust = "left", size = 4.7) +

  #formatting
  facet_wrap(~ population, nrow = 4, strip.position = "right", labeller = as_labeller(pop.names)) +
  labs(x = "Cumulative Degree Days (C)", y = "Daylength", color = "Year") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(color = 'black', fill="white", size = 1))

#export 718 x 710 (or 714 x 680)

# DD gained with plasticity ----
### for plastic CDL, find row in data frame below where daylength is less than weekly CDL in that row, find degree days of that row
### for constant CDL, find row where daylength is less than the constant CDL, find degree days of that row

DD.gained.data <-left_join(uspest %>% dplyr::select(population, dayofyear, week, year, daylength) %>% distinct(), plast_CDL) %>%
  left_join(ave.cum.dd)

#figuring out when CDL happens each year for plastic CDL

DD.gained.data.afterJuly <- DD.gained.data %>% filter(week > 26)

DD.gained.results_fromfn <- DD.gained.data.afterJuly[DD.gained.data.afterJuly$daylength < DD.gained.data.afterJuly$week.CDL, ] %>%
  group_by(population, year) %>%
  filter(rank(dayofyear, ties.method = 'first') == 1)
DD.gained.results_fromfn$population <- factor(DD.gained.results_fromfn$population, levels = c("Delta", 'StGeorge', 'BigBend', "Cibola") )#, labels = c("Delta", 'St. George', 'Big Bend', "Cibola"))

#when does CDL happen each year for constant CDL
DD.gained.data.afterJuly_withconstant <- left_join(DD.gained.data.afterJuly, bean_cdl_estimates)

DD.gained.results_fromfn_constant <- DD.gained.data.afterJuly_withconstant[DD.gained.data.afterJuly_withconstant$daylength < DD.gained.data.afterJuly_withconstant$cdl.2019, ] %>%
  group_by(population, year) %>%
  filter(rank(dayofyear, ties.method = 'first') == 1) %>%
  filter(population != "BigBend" | year != 2015 ) # remove Big Bend 2015, since we don't have plastic CDL data for that year
DD.gained.results_fromfn_constant$population <- factor(DD.gained.results_fromfn_constant$population, levels = c("Delta", 'StGeorge', 'BigBend', "Cibola") )#, labels = c("Delta", 'St. George', 'Big Bend', "Cibola"))


#Old results, from an average of 10 years
# DD.gained.results <- data.frame(
#   population = c("Delta", 'StGeorge', 'BigBend', "Cibola"),
#   DD.gained = c(-137.1,197.3,315.9,538.5),
#   plast.CDL = c(14.6422, 13.05374,11.571149,11.299),
#   plast.dayofyear = c(197, 243, 295, 302),
#   calendar.gained = c(-9, 11, 23, 40)
# )

#exploratory plot of difference between diapause timing with and without plasticity
ggplot() +
  geom_point(data = DD.gained.results_fromfn, aes(x = population, y = dayofyear),
             position = position_jitter(width = 0.2), color = 'blue') +
  geom_point(data = DD.gained.results_fromfn_constant, aes(x = population, y = dayofyear),
             position = position_jitter(width = 0.2), color = 'black')

#figure out difference between constant and plastic CDL
DDGainedFullResults <- data.frame(population = DD.gained.results_fromfn$population,
           year = DD.gained.results_fromfn$year,
           plastic.CDL = DD.gained.results_fromfn$week.CDL,
           plastic.DD = DD.gained.results_fromfn$cumulative.degree.days.C,
           plastic.day = DD.gained.results_fromfn$dayofyear,
           constant.CDL = DD.gained.results_fromfn_constant$week.CDL,
           constant.DD = DD.gained.results_fromfn_constant$cumulative.degree.days.C,
           constant.day = DD.gained.results_fromfn_constant$dayofyear,
           days.gained = DD.gained.results_fromfn$dayofyear - DD.gained.results_fromfn_constant$dayofyear,
           dd.gained = DD.gained.results_fromfn$cumulative.degree.days.C - DD.gained.results_fromfn_constant$cumulative.degree.days.C)
DDGainedFullResults$population <- factor(DDGainedFullResults$population, levels = c("Delta", 'StGeorge', 'BigBend', "Cibola"), labels = c("Delta", 'St. George', 'Big Bend', "Cibola"))

#run a model to get model estimates and 95% CI for degree days gained
mod.ddgained <- lmer(dd.gained ~ population + (1|year), data = DDGainedFullResults)
summary(mod.ddgained)
emmeans(mod.ddgained, pairwise ~ population)
ddgained.emout <- as.data.frame(emmeans(mod.ddgained, pairwise ~ population)$emmeans)

#model for calendar days gained
mod.daysgained <- lmer(days.gained ~ population + (1|year), data = DDGainedFullResults)
summary(mod.daysgained)
emmeans(mod.daysgained, pairwise ~ population)
daysgained.emout <- as.data.frame(emmeans(mod.daysgained, pairwise ~ population)$emmeans)

ggplot() +
  geom_point(data = DDGainedFullResults, aes(x = population, y = dd.gained, color = population),
             size = 3, alpha = 0.5, position = position_jitter(width = 0.2)) +
  geom_point(data = ddgained.emout, aes(x = population, y = emmean, color = population),
                  size = 5) +
  geom_errorbar(data = ddgained.emout, aes(x = population, ymin = lower.CL, ymax = upper.CL, color = population),
                width = 0, size = 1.2) +
  #labels
  geom_text(data = ddgained.emout, aes(x = population, y = 0), position = position_nudge(y = -30),
            label = c("-11 days" ,"10 days", '19 days', '35 days'),color = 'black') +
  geom_text(data = ddgained.emout, aes(x = population, y = 0), position = position_nudge(y = 55),
            label = c("-171 deg-days" ,"194 deg-days", '280 deg-days', '485 deg-days'),color = 'black') +
  geom_hline(yintercept = 0, color = 'grey 75') +
  #format
  labs(x = "Population", y = "Degree Days Gained", color = "Population") +
  scale_color_viridis_d() +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(), legend.position = 'none')
#export 500 x 400



################OLDER VERSIONS#####################
# Photothermographs with USPEST Data ----
Delta_DD_2018 <- read.csv("Deg_days_USPESTmodel//Diorhabda_carinulata_CMP08_Delta_UT_2018_1_1.csv", skip = 1)
Delta_DD_2018$month <- as.numeric(Delta_DD_2018$month)
Delta_DD_2018$day <- as.numeric(Delta_DD_2018$day)
Delta_DD_2018 <- Delta_DD_2018 %>% mutate(year = rep(2018, nrow(Delta_DD_2018)), .before = 'month') %>%
  unite(col = 'date', year:day, sep = '-', remove = F)
Delta_DD_2018$date <- ymd(Delta_DD_2018$date)
Delta_DD_2018 <- Delta_DD_2018 %>% mutate(daylength = daylength(lat = 39.100040, doy = date))

str(Delta_DD_2018)

first_ofthe_month <- Delta_DD_2018 %>% filter(day == 1)

ggplot() +
  geom_line(data = Delta_DD_2018, aes(x = cumulative.degree.days.C, y = daylength)) +
  geom_point(data = first_ofthe_month, aes(x = cumulative.degree.days.C, y = daylength)) +
  geom_text(data = first_ofthe_month, aes(x = cumulative.degree.days.C, y = daylength),
            label = month.letters, nudge_y = 0.2, nudge_x = 25)

ggplot(Delta_DD_2018, aes(x = maximum.temperature.C, y = daylength)) +
  geom_line()



# Photothermographs with NOAA temp Data ----

Delta_2018_temp <- Delta %>% filter(YEAR == 2018)
Delta_2018_temp$today_DD <- dd_sng_sine(daily_min = Delta_2018_temp$TMIN, daily_max = Delta_2018_temp$TMAX, thresh_low = 11.11, thresh_up = 36.7, cumulative = F)
Delta_2018_temp$cum_DD <- dd_sng_sine(daily_min = Delta_2018_temp$TMIN, daily_max = Delta_2018_temp$TMAX, thresh_low = 11.11, thresh_up = 36.7, cumulative = T)

daily <- dd_sng_sine(daily_min = Delta_2018_temp$TMIN, daily_max = Delta_2018_temp$TMAX, thresh_low = 11.11, thresh_up = 36.7, cumulative = F)

# Degree day and daylength calculations ----
## Delta
Delta.latitude = 39.100040
Delta_photothermo_data <- DD_allyears(Delta, "Delta", Delta.latitude)
Delta_photothermo_data$YEAR <- as.factor(Delta_photothermo_data$YEAR)

Delta_photothermo_data %>% group_by(YEAR) %>%
  summarise(
    total_DD <- sum(daily_DD)
  )

##St. George
StGeorge.latitude = 37.1
StGeorge_photothermo_data <- DD_allyears(StGeorge, "StGeorge", StGeorge.latitude)
StGeorge_photothermo_data$YEAR <- as.factor(StGeorge_photothermo_data$YEAR)

StGeorge_photothermo_data %>% group_by(YEAR) %>%
  summarise(
    total_DD <- sum(daily_DD)
  )


# Photothermographs ----
## Delta
first_ofthe_month_noaa <- Delta_photothermo_data %>% filter(DAY == 1)

ggplot() +
  geom_line(data = Delta_photothermo_data, aes(x = cumulative_DD, y = daylength, color = YEAR)) +
  geom_point(data = first_ofthe_month_noaa, aes(x = cumulative_DD, y = daylength)) +
  # geom_text(data = first_ofthe_month_noaa, aes(x = cumulative_DD, y = daylength),
  #           label = month.letters, nudge_y = 0.2, nudge_x = 25) +
  theme_bw()
#2019,2015, 2016

## StGeorge
StG_first_ofthe_month <- StGeorge_photothermo_data %>% filter(DAY == 1)

ggplot() +
  geom_line(data = StGeorge_photothermo_data, aes(x = cumulative_DD, y = daylength, color = YEAR)) +
  geom_point(data = StG_first_ofthe_month, aes(x = cumulative_DD, y = daylength)) +
  # geom_text(data = first_ofthe_month, aes(x = cumulative_DD, y = daylength),
  #           label = month.letters, nudge_y = 0.2, nudge_x = 25) +
  theme_bw()

# Functions----
## Degree-days functions ----
Diorhabda_DD <- function(data) {
  daily_DD <- dd_calc(daily_min = data$TMIN, daily_max = data$TMAX, thresh_low = 11.11, thresh_up = 36.7,
                   method = "sng_sine", cumulative = F, quiet = T, interpolate_na = T)
  cumulative_DD <- dd_calc(daily_min = data$TMIN, daily_max = data$TMAX, thresh_low = 11.11, thresh_up = 36.7,
                        method = "sng_sine", cumulative = T, quiet = T, interpolate_na = T)
  out <- cbind(data, daily_DD, cumulative_DD)
  return(out)
}

DD_allyears <- function(data, site, latitude){
  listyear <- split(data, data$YEAR)
  out <- lapply(listyear, Diorhabda_DD)
  df <- bind_rows(out)
  out2 <- cbind(site = rep(site, nrow(df)), df) %>%
    mutate(daylength = daylength(lat = latitude, doy = dayofyear))
  return(out2)
}


