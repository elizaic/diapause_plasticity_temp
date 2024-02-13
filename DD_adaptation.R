library(lubridate)
library(geosphere)
library(degday)
library(RColorBrewer)
library(ggnewscale)
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
         ) %>% mutate(pop.index = case_when (
    population == "Delta" ~ 1,
    population == "StGeorge" ~ 2,
    population == "BigBend" ~ 3,
    population == "Cibola" ~ 4,
  ), .after = "population") %>%
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

#Remaining degree days after CDL
#total DD each year
totalDD <- uspest %>% group_by(population, year) %>% summarize(
  total.DD = max(cumulative.degree.days.C)) %>%
  filter(population != "BigBend" | year != 2015) # remove Big Bend 2015, since we don't have plastic CDL data for that year
totalDD_plastic <- left_join(DD.gained.results_fromfn, totalDD)
totalDD_constant <- left_join(DD.gained.results_fromfn_constant, totalDD)

# totalDD$population <- factor(totalDD$population, levels = c("Delta", 'StGeorge', 'BigBend', "Cibola") )#, labels = c("Delta", 'St. George', 'Big Bend', "Cibola"))


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
             position = position_jitter(width = 0.2, height = 0), color = 'blue') +
  geom_point(data = DD.gained.results_fromfn_constant, aes(x = population, y = dayofyear),
             position = position_jitter(width = 0.2, height = 0), color = 'black')

#figure out difference between constant and plastic CDL
DDGainedFullResults <- data.frame(population = DD.gained.results_fromfn$population,
           year = DD.gained.results_fromfn$year,
           plastic.CDL = DD.gained.results_fromfn$week.CDL,
           plastic.DD = DD.gained.results_fromfn$cumulative.degree.days.C,
           plastic.day = DD.gained.results_fromfn$dayofyear,
           plastic.DDremaining = totalDD_plastic$total.DD - totalDD_plastic$cumulative.degree.days.C,
           constant.CDL = DD.gained.results_fromfn_constant$week.CDL,
           constant.DD = DD.gained.results_fromfn_constant$cumulative.degree.days.C,
           constant.day = DD.gained.results_fromfn_constant$dayofyear,
           constant.DDremaining = totalDD_constant$total.DD - totalDD_constant$cumulative.degree.days.C,
           dd.leftover.reduced = (totalDD_constant$total.DD - totalDD_constant$cumulative.degree.days.C) - (totalDD_plastic$total.DD - totalDD_plastic$cumulative.degree.days.C),
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

#model for days remaining plastic
mod.daysremaining <- lmer(plastic.DDremaining ~ population + (1|year), data = DDGainedFullResults)
summary(mod.daysremaining)
emmeans(mod.daysremaining, pairwise ~ population)
# daysgained.emout <- as.data.frame(emmeans(mod.daysremaining, pairwise ~ population)$emmeans)

#model for days remaining - constant
mod.daysremaining.constant <- lmer(constant.DDremaining ~ population + (1|year), data = DDGainedFullResults)
summary(mod.daysremaining.constant)
emmeans(mod.daysremaining.constant, pairwise ~ population)
# daysgained.emout <- as.data.frame(emmeans(mod.daysremaining.constant, pairwise ~ population)$emmeans)

#model for plastic day
mod.plasticday <- lmer(plastic.day ~ population + (1|year), data = DDGainedFullResults)
summary(mod.plasticday)
emmeans(mod.plasticday, pairwise ~ population)
plasticday.emout <- as.data.frame(emmeans(mod.plasticday, pairwise ~ population)$emmeans)
# plasticday.emout$population <- factor(plasticday.emout$population, levels = c("Delta", "StGeorge", "BigBend", "Cibola"), labels = c("Delta", 'St. George', 'Big Bend', "Cibola"))

#model for leftover dd reduced
mod.leftoverdd <- lmer(dd.leftover.reduced ~ population + (1|year), data = DDGainedFullResults)
summary(mod.leftoverdd)
emmeans(mod.leftoverdd, pairwise ~ population)


daysGainedPlot <- ggplot() +
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
  labs(x = "Population", y = "Degree Days Gained\nwith Plasticity", color = "Population") +
  scale_color_viridis_d() +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(), legend.position = 'none')
daysGainedPlot
#export 500 x 400

ggarrange(daysGainedPlot, frostVSCDL, labels = "AUTO")
#export 1010 x 475

#summaries of days gained for photothermographs
average.ddgained <- data.frame(population = c("Delta", 'StGeorge', 'BigBend', "Cibola"),
                               dd.gained = c(c("-171" ,"194", '280', '485')))
average.ddgained$population <- factor(average.ddgained$population, levels = c("Delta", 'StGeorge', 'BigBend', "Cibola"))


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
  scale_color_manual(values = colorRampPalette(brewer.pal(8, "Blues"))(11), guide = 'none' )+

  #plastic CDL lines
  geom_line(data = plast_CDL, aes(x = week.cumDD, y = week.CDL, alpha = year), size = .7) +
  # new_scale_colour() +
  # scale_color_manual(values = colorRampPalette(brewer.pal(8, "Greys"))(11) )+
  geom_text(data = plast_CDL %>% filter(week == 53 & population == "Delta" & year == 2020),
            aes(x = week.cumDD, y = week.CDL, group = population),
            label = "Plastic CDL", hjust = "left", nudge_y = 0.7, nudge_x = 100) +

  #constant CDL lines
  geom_hline(data = bean_cdl_estimates, aes(yintercept = cdl.2019), linetype = 2, size = 1) +
  geom_text(data = bean_cdl_estimates %>% filter(population == "Delta"),
            aes(x = 4500, y = cdl.2019, group = population),
            label = "Non-plastic CDL", hjust = "center", nudge_y = -0.5) +

  #Month labels
  geom_point(data = first_ofthe_month_all, aes(x = cumulative.degree.days.C, y = daylength)) +
  geom_text(data = first_ofthe_month_all, aes(x = cumulative.degree.days.C, y = daylength),
            label = rep(month.letters,4), nudge_y = -0.5, nudge_x = 60) +

  #DD gained labels
  geom_text(data = average.ddgained, aes(x = 820, y = 19, group = population),
            label = "Degree days gained with plasticity:", hjust = "left", size = 4.5) +
  geom_text(data = average.ddgained, aes(x = 3610, y = 19, group = population),
            label = average.ddgained$dd.gained, hjust = "left", size = 4.7) +

  #formatting
  facet_wrap(~ population, nrow = 4, strip.position = "right", labeller = as_labeller(pop.names)) +
  labs(x = "Cumulative Degree Days (C)", y = "Daylength", alpha = "Year") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(color = 'black', fill="white", size = 1))


#export 718 x 710 (or 714 x 680)


# If high (cibola) plasticity is adaptive at all populations----

##write a function for a line that represents how CDL changes with temperature

Cibola_plast_CDL <- plast_CDL %>% mutate(week.CDL = 0)

for (i in 1:nrow(Cibola_plast_CDL)) {
  Cibola_plast_CDL[i,6] <- as.numeric(Cibola_plast_CDL[i,5]) * plast.estimates[4,4] + plast.estimates[as.numeric(plast_CDL[i,2]),3]
}
Cibola_plast_CDL$week.CDL <- as.numeric(Cibola_plast_CDL$week.CDL)
str(Cibola_plast_CDL)

#Plot the weekly CDLs
ggplot(data = Cibola_plast_CDL, aes(x = week, y = week.CDL, color = population, linetype = year)) +
  geom_line()

#relate weekly CDL to cumulative degree days (join week.cum.dd dataframe to plastic CDLs)
Cibola_plast_CDL <- left_join(x = Cibola_plast_CDL, y = week.cum.dd)


#Days gained when all pops have plasticity of Cibola

Cibola_DD.gained.data <-left_join(uspest %>% dplyr::select(population, dayofyear, week, year, daylength) %>% distinct(), Cibola_plast_CDL) %>%
  left_join(ave.cum.dd)

#figuring out when CDL happens each year for plastic CDL

Cibola_DD.gained.data.afterJuly <- Cibola_DD.gained.data %>% filter(week > 26)

Cibola_DD.gained.results_fromfn <- Cibola_DD.gained.data.afterJuly[Cibola_DD.gained.data.afterJuly$daylength < Cibola_DD.gained.data.afterJuly$week.CDL, ] %>%
  group_by(population, year) %>%
  filter(rank(dayofyear, ties.method = 'first') == 1)
Cibola_DD.gained.results_fromfn$population <- factor(Cibola_DD.gained.results_fromfn$population, levels = c("Delta", 'StGeorge', 'BigBend', "Cibola") )#, labels = c("Delta", 'St. George', 'Big Bend', "Cibola"))

#when does CDL happen each year for constant CDL
Cibola_DD.gained.data.afterJuly_withconstant <- left_join(Cibola_DD.gained.data.afterJuly, bean_cdl_estimates)

Cibola_DD.gained.results_fromfn_constant <- Cibola_DD.gained.data.afterJuly_withconstant[Cibola_DD.gained.data.afterJuly_withconstant$daylength < Cibola_DD.gained.data.afterJuly_withconstant$cdl.2019, ] %>%
  group_by(population, year) %>%
  filter(rank(dayofyear, ties.method = 'first') == 1) #%>%
  # filter(population != "BigBend" | year != 2015 ) # remove Big Bend 2015, since we don't have plastic CDL data for that year
Cibola_DD.gained.results_fromfn_constant$population <- factor(Cibola_DD.gained.results_fromfn_constant$population, levels = c("Delta", 'StGeorge', 'BigBend', "Cibola") )#, labels = c("Delta", 'St. George', 'Big Bend', "Cibola"))

#Remaining degree days after CDL
#total DD each year
Cibola_totalDD <- uspest %>% group_by(population, year) %>% summarize(
  total.DD = max(cumulative.degree.days.C)) #%>%
  # filter(population != "BigBend" | year != 2015) # remove Big Bend 2015, since we don't have plastic CDL data for that year
Cibola_totalDD_plastic <- left_join(Cibola_DD.gained.results_fromfn, Cibola_totalDD)
Cibola_totalDD_constant <- left_join(Cibola_DD.gained.results_fromfn_constant, Cibola_totalDD)

#exploratory plot of difference between diapause timing with and without plasticity
ggplot() +
  geom_point(data = Cibola_DD.gained.results_fromfn, aes(x = population, y = dayofyear),
             position = position_jitter(width = 0.2, height = 0), color = 'blue') +
  geom_point(data = Cibola_DD.gained.results_fromfn_constant, aes(x = population, y = dayofyear),
             position = position_jitter(width = 0.2, height = 0), color = 'black')

#figure out difference between constant and plastic CDL
Cibola_DDGainedFullResults <- data.frame(population = Cibola_DD.gained.results_fromfn$population,
                                  year = Cibola_DD.gained.results_fromfn$year,
                                  # plastic.CDL = DD.gained.results_fromfn$week.CDL,
                                  # plastic.DD = DD.gained.results_fromfn$cumulative.degree.days.C,
                                  plastic.day = Cibola_DD.gained.results_fromfn$dayofyear,
                                  # plastic.DDremaining = totalDD_plastic$total.DD - totalDD_plastic$cumulative.degree.days.C,
                                  # constant.CDL = DD.gained.results_fromfn_constant$week.CDL,
                                  # constant.DD = DD.gained.results_fromfn_constant$cumulative.degree.days.C,
                                  # constant.day = DD.gained.results_fromfn_constant$dayofyear,
                                  # constant.DDremaining = totalDD_constant$total.DD - totalDD_constant$cumulative.degree.days.C,
                                  # dd.leftover.reduced = (totalDD_constant$total.DD - totalDD_constant$cumulative.degree.days.C) - (totalDD_plastic$total.DD - totalDD_plastic$cumulative.degree.days.C),
                                  # days.gained = DD.gained.results_fromfn$dayofyear - DD.gained.results_fromfn_constant$dayofyear,
                                  dd.gained = Cibola_DD.gained.results_fromfn$cumulative.degree.days.C - Cibola_DD.gained.results_fromfn_constant$cumulative.degree.days.C)
Cibola_DDGainedFullResults$population <- factor(Cibola_DDGainedFullResults$population, levels = c("Delta", 'StGeorge', 'BigBend', "Cibola"), labels = c("Delta", 'St. George', 'Big Bend', "Cibola"))


ggplot() +
  geom_point(data = Cibola_DDGainedFullResults, aes(x = population, y = dd.gained, color = population),
             size = 3, alpha = 0.5, position = position_jitter(width = 0.2)) +
  # geom_point(data = ddgained.emout, aes(x = population, y = emmean, color = population),
  #            size = 5) +
  # geom_errorbar(data = ddgained.emout, aes(x = population, ymin = lower.CL, ymax = upper.CL, color = population),
  #               width = 0, size = 1.2) +
  #labels
  # geom_text(data = ddgained.emout, aes(x = population, y = 0), position = position_nudge(y = -30),
  #           label = c("-11 days" ,"10 days", '19 days', '35 days'),color = 'black') +
  # geom_text(data = ddgained.emout, aes(x = population, y = 0), position = position_nudge(y = 55),
  #           label = c("-171 deg-days" ,"194 deg-days", '280 deg-days', '485 deg-days'),color = 'black') +
  geom_hline(yintercept = 0, color = 'grey 75') +
  #format
  labs(x = "Population", y = "Degree Days Gained\nwith Plasticity", color = "Population") +
  scale_color_viridis_d() +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(), legend.position = 'none')

#export 500 x 400

#frost vs. diapause figure
ggplot() +
  # First Frost
  geom_pointrange(data = emout2, mapping = aes(x = site, y = emmean, ymin = lower.CL, ymax = upper.CL),
                  size = 1, shape = 15, color = 'cornflowerblue') +
  geom_jitter(data = firstFrost2022, mapping = aes(x = site, y = firstfrostday),
              width = 0.2, height = 0, alpha = 0.4, shape = 15, size = 2, color = 'cornflowerblue') +
  annotate(geom = 'text', x = 1, y = 305, label = "First Frost", size = 5) +

  # plastic CDL
  # geom_pointrange(data = plasticday.emout, mapping = aes(x = population, y = emmean, ymin = lower.CL, ymax = upper.CL, color = population),
                  # size = 1, shape = 16) +
  geom_point(data = Cibola_DDGainedFullResults, aes(x = population, y = plastic.day, color = population),
             position = position_jitter(width = 0.2, height = 0), size = 2, shape = 16, alpha = 0.4) +
  # annotate(geom = 'text', x = 1.55, y = 200, label = "Diapause\ninitiation\n(CDL)", size = 5) +

  # formatting
  scale_y_continuous(name = "Day of Year", breaks = month.breaks, labels = month.name) +
  scale_x_discrete(name = "Population") +
  scale_color_viridis_d(guide = 'none') +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank())


#Correlations of frost and temp and daylength -----

#Correlations and variation within sites for whole year
uspest %>% group_by(population) %>% summarize(
  mean.max = mean(maximum.temperature.C),
  sd.max = sd(maximum.temperature.C),
  mean.min = mean(minimum.temperature.C),
  sd.min = sd(minimum.temperature.C),
  CV.max = sd.max/mean.max * 100,
  CV.min = sd.min/mean.min * 100
)

#CV for only fall
uspest %>% filter(month > 6) %>% group_by(population) %>% summarize(
  mean.max = mean(maximum.temperature.C),
  sd.max = sd(maximum.temperature.C),
  mean.min = mean(minimum.temperature.C),
  sd.min = sd(minimum.temperature.C),
  CV.max = sd.max/mean.max * 100,
  CV.min = sd.min/mean.min * 100
)


# First frost with uspest data

uspest.aferJuly <- uspest %>% filter(week > 26)

frost.uspest <- uspest.aferJuly[uspest.aferJuly$minimum.temperature.C < 0, ] %>%
  group_by(population, year) %>%
  filter(rank(dayofyear, ties.method = 'first') == 1)
frost.uspest$population <- factor(frost.uspest$population, levels = c("Delta", 'StGeorge', 'BigBend', "Cibola") )#, labels = c("Delta", 'St. George', 'Big Bend', "Cibola"))

#Getting upstream (14 days) dates and temps those days
frost.uspest.upstreamdates <- data.frame(population = frost.uspest$population,
                                         pop.index = frost.uspest$pop.index,
                                         year = frost.uspest$year,
                                         dayofyear = frost.uspest$dayofyear)
frost.uspest.upstreamdates <- frost.uspest.upstreamdates %>% mutate(upstream.start = frost.uspest$dayofyear-15) %>%
  mutate(upstream.end = frost.uspest$dayofyear-1)

upstream.data <- data.frame()
for (i in 1:nrow(frost.uspest)) {
  dat <- uspest %>% filter(pop.index == as.numeric(frost.uspest.upstreamdates[i,2]) &
                      year == frost.uspest.upstreamdates[i,3] &
                      between(dayofyear, frost.uspest.upstreamdates[i,5], frost.uspest.upstreamdates[i,6]))
  dat2 <- data.frame(dat)
  upstream.data <- rbind(upstream.data, dat2)
}

upstream.summary <- upstream.data %>% group_by(population, year) %>%
  summarise(
    ave.high = mean(maximum.temperature.C),
    ave.low = mean(minimum.temperature.C)
  )
frost.uspest <- left_join(frost.uspest, upstream.summary)

#Autocorrelation of temperature for each population
lag = 31

auto.cors <- rbind(
acf(uspest %>% filter(population == 'Delta') %>% select(maximum.temperature.C), lag.max = lag, pl = F)$acf,
acf(uspest %>% filter(population == 'StGeorge') %>% select(maximum.temperature.C), lag.max = lag, pl = F)$acf,
acf(uspest %>% filter(population == 'BigBend') %>% select(maximum.temperature.C), lag.max = lag, pl = F)$acf,
acf(uspest %>% filter(population == 'Cibola') %>% select(maximum.temperature.C), lag.max = lag, pl = F)$acf,

acf(uspest %>% filter(population == 'Delta') %>% select(minimum.temperature.C), lag.max = lag, pl = F)$acf,
acf(uspest %>% filter(population == 'StGeorge') %>% select(minimum.temperature.C), lag.max = lag, pl = F)$acf,
acf(uspest %>% filter(population == 'BigBend') %>% select(minimum.temperature.C), lag.max = lag, pl = F)$acf,
acf(uspest %>% filter(population == 'Cibola') %>% select(minimum.temperature.C), lag.max = lag, pl = F)$acf
) %>% as.data.frame() %>%
  mutate(population = rep(c('Delta', 'StGeorge', 'BigBend', 'Cibola'),2), .before = 'V1') %>%
  mutate(min.max.temp = c(rep(c('max'),4), rep(c('min'),4)), .before = 'V1')

library(matrixStats)
week1auto = rowMeans(auto.cors[ , c(7: 13)])
week1auto.sd = rowSds(as.matrix(auto.cors[ , c(7: 13)]))
week2auto = rowMeans(auto.cors[ , c(14: 20)])
week2auto.sd = rowSds(as.matrix(auto.cors[ , c(14: 20)]))
week3auto = rowMeans(auto.cors[ , c(21: 27)])
week3auto.sd = rowSds(as.matrix(auto.cors[ , c(21: 27)]))
week4auto = rowMeans(auto.cors[ , c(28: 34)])
week4auto.sd = rowSds(as.matrix(auto.cors[ , c(28: 34)]))

auto.cors <- cbind(auto.cors, week1auto, week1auto.sd, week2auto, week2auto.sd, week3auto, week3auto.sd, week4auto, week4auto.sd)

ggplot(data = auto.cors, aes(x = population, y = week1auto)) +
  geom_point()

lm(dayofyear ~ ave.low * daylength, data = frost.uspest %>% filter(population == 'Cibola')) %>%
  # emtrends(~population, var = 'ave.low')


lm(dayofyear ~ ave.low * daylength, data = frost.uspest %>% filter(population == 'Cibola')) %>%
  # emtrends(~population, var = 'ave.low')
  summary()
  # Anova(type = 3)
lm(dayofyear ~ ave.low * daylength, data = frost.uspest %>% filter(population == 'Delta')) %>%
  summary()
lm(dayofyear ~ ave.low * daylength, data = frost.uspest %>% filter(population == 'StGeorge')) %>%
  summary()



# frost.uspest.upstreamdates = NA
# for (i in 1:14) {
#   frost.uspest.upstreamdates <- cbind(frost.uspest.upstreamdates, frost.uspest$date - i)
# }
# frost.uspest.upstreamdates <- cbind(frost.uspest, frost.uspest.upstreamdates)
#
# frost.uspest.upstreamdates$V2 <- as_date(frost.uspest.upstreamdates$V2, origin = lubridate::origin)
# frost.uspest.upstreamdates$V3 <- as_date(frost.uspest.upstreamdates$V3, origin = lubridate::origin)
# frost.uspest.upstreamdates$V4 <- as_date(frost.uspest.upstreamdates$V4, origin = lubridate::origin)
# frost.uspest.upstreamdates$V5 <- as_date(frost.uspest.upstreamdates$V5, origin = lubridate::origin)
# frost.uspest.upstreamdates$V6 <- as_date(frost.uspest.upstreamdates$V6, origin = lubridate::origin)
# frost.uspest.upstreamdates$V7 <- as_date(frost.uspest.upstreamdates$V7, origin = lubridate::origin)
# frost.uspest.upstreamdates$V8 <- as_date(frost.uspest.upstreamdates$V8, origin = lubridate::origin)
# frost.uspest.upstreamdates$V9 <- as_date(frost.uspest.upstreamdates$V9, origin = lubridate::origin)
# frost.uspest.upstreamdates$V10 <- as_date(frost.uspest.upstreamdates$V10, origin = lubridate::origin)
# frost.uspest.upstreamdates$V11 <- as_date(frost.uspest.upstreamdates$V11, origin = lubridate::origin)
# frost.uspest.upstreamdates$V12 <- as_date(frost.uspest.upstreamdates$V12, origin = lubridate::origin)
# frost.uspest.upstreamdates$V13 <- as_date(frost.uspest.upstreamdates$V13, origin = lubridate::origin)
# frost.uspest.upstreamdates$V14 <- as_date(frost.uspest.upstreamdates$V14, origin = lubridate::origin)
# frost.uspest.upstreamdates$V15 <- as_date(frost.uspest.upstreamdates$V15, origin = lubridate::origin)
#
#
# uspest %>% filter(pop.index == as.numeric(frost.uspest.upstreamdates[1,3]) & between(date, as_date(frost.uspest.upstreamdates[1,20]), frost.uspest.upstreamdates[1,33]))
#
#                   #)date == as_date(frost.uspest.upstreamdates[1,20]))
#
#
#
#
#
#  kmnmn<- uspest %>% group_by(population, year) %>% filter
#   summarise(
#
#   )
#
# before.frost.uspest <- uspest.aferJuly[frost.uspest$dayofyear, ] %>%
#   group_by(population, year)
#
#
# before.frost.uspest = data.frame(NA)
#
# for (j in 1:nrow(frost.uspest)) {
#   for (i in (as.numeric(frost.uspest[j,7] - 14)) : (as.numeric(frost.uspest[j,7]))) {
#
#     rbind(before.frost.uspest, uspest.aferJuly[i,])
#     }
# }
#
#
#
#
# frost.uspest <- frost.uspest %>% mutate(upstream.date = frost.uspest$date - 14)

lm(dayofyear ~ maximum.temperature.C + daylength + population + maximum.temperature.C:population + daylength:population, data = frost.uspest %>% filter(population != 'BigBend')) %>% Anova(type = 3)





cor.test(~ daylength + ave_max, data = uspest, subset = site == "Delta (39°N)")
cor.test(~ daylength + ave_max, data = daylength.temp, subset = site == "St. George (37°N)")
cor.test(~ daylength + ave_max, data = daylength.temp, subset = site == "Big Bend (35°N)")
cor.test(~ daylength + ave_max, data = daylength.temp, subset = site == "Cibola (33°N)")

lm(maximum.temperature.C ~ daylength * population, data = uspest) %>%
  summary()
emmeans()





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


