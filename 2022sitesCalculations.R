library(tidyverse)
library(car)
library(emmeans)
library(ggpubr)
source("functions.R")

#Constants/other stuff for plotting
month.breaks = c(1,32,60,91,121,152,182,213,244,274,305,335)+15
month.letters = c('J','F','M','A','M','J','J','A','S','O','N','D')

#Load and Prep Data ----
#OLDER Datasets, incomplete cibola
# dataAllSites <- read.csv("3048279.csv")
dataAllSites <- read.csv("3378358_downloaded_6_29_23.csv")

#Preliminary Clean up
dataAllSites <- drop_na(dataAllSites, TMIN)
dataAllSites$NAME <- as.factor(dataAllSites$NAME)
dataAllSites <- dataAllSites %>% mutate(NAME_group = if_else(dataAllSites$LATITUDE > 39, "Delta",
                                       if_else(dataAllSites$LATITUDE < 39 & dataAllSites$LATITUDE > 37, "StGeorge",
                                               if_else(dataAllSites$LATITUDE < 37 & dataAllSites$LATITUDE > 35, "BigBend",
                                                       if_else(dataAllSites$LATITUDE <35, "Cibola", "XXXXXXXXX")))))
dataAllSites$NAME_group <- as.factor(dataAllSites$NAME_group)
dataAllSites$NAME_group <- factor(dataAllSites$NAME_group, levels = c("Delta", "StGeorge", "BigBend", "Cibola"))
unique(dataAllSites$NAME)
unique(dataAllSites$NAME_group)
dataAllSites <- prep1(dataAllSites) #add day of year, year, and month columns to dataset

#Split up dataset into sites
listAllSites <- split(dataAllSites, dataAllSites$NAME_group)
# names(listAllSites) <- c("Big Bend", "Delta", "Blythe (Cibola)", "St George")
Delta <- listAllSites[[1]]
StGeorge <- listAllSites[[2]]
BigBend <- listAllSites[[3]]
Cibola <- listAllSites[[4]]

#First Frost ----
#Get first frost of each site, put in one data frame
firstFrost2022 <- rbind(frostNEW(Delta, "Delta"),
      frostNEW(StGeorge, "StGeorge"),
      frostNEW(BigBend, "BigBend"),
      frostNEW(Cibola, "Cibola")
)
firstFrost2022$site <- factor(firstFrost2022$site, levels = c("Delta", "StGeorge", "BigBend", "Cibola"), labels = c("Delta", 'St. George', 'Big Bend', "Cibola"))

summaryFirstFrost <- firstFrost2022 %>%
  group_by(site) %>%
  summarize(meanDay = mean(firstfrostday, na.rm = TRUE), meanTemp = mean(firstfrosttemp, na.rm = TRUE))

## first frost day
mod <- lm(firstfrostday ~ site + year, data = firstFrost2022)
Anova(mod, type = 3)
plot(mod)
emout <- emmeans(mod, pairwise ~ site)
emout2 <- as.data.frame(emout$emmean)
emout2$site <- factor(emout2$site, levels = c("Delta", "StGeorge", "BigBend", "Cibola"), labels = c("Delta", 'St. George', 'Big Bend', "Cibola"))

frostVSCDL <- ggplot() +
  # First Frost
  geom_pointrange(data = emout2, mapping = aes(x = site, y = emmean, ymin = lower.CL, ymax = upper.CL, color = site), size = 1) +
  geom_jitter(data = firstFrost2022, mapping = aes(x = site, y = firstfrostday, color = site),
              width = 0.2, height = 0, alpha = 0.4, shape = 16, size = 2) +
  annotate(geom = 'text', x = 1, y = 305, label = "First Frost", size = 5) +

  # plastic CDL
  geom_pointrange(data = plasticday.emout, mapping = aes(x = population, y = emmean, ymin = lower.CL, ymax = upper.CL, color = population),
                  size = 1, shape = 15) +
  geom_point(data = DDGainedFullResults, aes(x = population, y = plastic.day, color = population),
             position = position_jitter(width = 0.2, height = 0), size = 2, shape = 15, alpha = 0.4) +
  annotate(geom = 'text', x = 1.55, y = 200, label = "Diapause\ninitiation\n(CDL)", size = 5) +

  # formatting
  scale_y_continuous(name = "Day of Year", breaks = month.breaks, labels = month.name) +
  scale_x_discrete(name = "Population") +
  scale_color_viridis_d(guide = 'none') +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank())
frostVSCDL
#export 454 x 496

##OLDER analyses ----

#warm & cool CDL with First frost OLD ----
CDL_estimates <- read.csv("CDL_estimates.csv")
CDL_estimates$population <- factor(CDL_estimates$population, levels = c("De", "Sg", "Bi", "Ci"), labels = c("Delta", "StGeorge", "BigBend", "Cibola") )
CDL_estimates$temperature <- factor(CDL_estimates$temperature, levels = c("38", "28"))

ggplot() +
  geom_pointrange(data = emout2, mapping = aes(x = site, y = emmean, ymin = lower.CL, ymax = upper.CL), size = 1) +
  geom_jitter(data = firstFrost2022, mapping = aes(x = site, y = firstfrostday),
              width = .15, height = 0, alpha = .4) +
  geom_jitter(data = CDL_estimates, aes(x = population, y = dayofyear, color = temperature),
             shape = 17, size = 3, width = 0.1) +
  scale_y_continuous(name = "First frost day +/- CI", breaks = month.breaks, labels = month.name) +
  scale_x_discrete(name = "Site") +
  labs(color = "Temp.\nof CDL") +
  theme_bw(base_size = 15)

#Difference between CDL and first frost - OLD ----
frost.CDL.diffs <- read.csv("CDL_resultstable_diff_frostCDL.csv")
frost.CDL.diffs$Population <- factor(frost.CDL.diffs$Population, levels = c("Delta", "StGeorge", "BigBend", "Cibola"))
frost.CDL.diffs$Temperature <- factor(frost.CDL.diffs$Temperature, levels = c("warm", "cool"))

ggplot() +
  geom_jitter(data = frost.CDL.diffs, aes(x = Population, y = Diff.frost.dayofCDL.estim,
                                         color = Temperature), width = 0.1, size = 3) +
  scale_y_continuous(name = "Days between diapause initiation (CDL date)\nand average first frost date (degree of mismatch)") +
  theme_bw(base_size = 15)

#Delta Temp Plots ----
# plot with max and min temps
ggplot(data = Delta %>% filter(YEAR == c(2018,2020,2021))) +
  geom_line(aes(x = dayofyear, y = TMIN, group = YEAR, color = factor(YEAR)), size = 1) +
  geom_line(aes(x = dayofyear, y = TMAX, group = YEAR, color = factor(YEAR)), size = 1) +
  # facet_grid( ~ site) +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank()) +
  # scale_color_discrete(name = "Year") +
  scale_color_brewer(name = "Year", palette = "Set2") +
  scale_x_continuous(name = "Month", breaks = month.breaks, labels = month.letters) +
  scale_y_continuous(name = "Daily Low and High Temperature (C)")

ggplot(data = Delta %>% filter(YEAR == 2021)) +
  geom_line(aes(x = dayofyear, y = TMAX), color = "indianred1", size = 1) +
  geom_line(aes(x = dayofyear, y = TMIN), color = "deepskyblue", size = 1) +
  # facet_grid( ~ site) +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank()) +
  # scale_color_discrete(name = "Year") +
  scale_x_continuous(name = "Month", breaks = month.breaks, labels = month.letters) +
  scale_y_continuous(name = "Daily High and Low Temperature (C)")

#Days above 38 degrees, by year
Delta %>% filter(TMAX > 38) %>% group_by(YEAR) %>% summarize(
  count = n()
)

#Ave temps ----
#Ave min and max temp calculations

aveTemp <- rbind(
  Delta.summary <- avetempsfn(Delta, "Delta"),
  StGeorge.summary <- avetempsfn(StGeorge, "StGeorge"),
  BigBend.summary <- avetempsfn(BigBend, "BigBend"),
  Cibola.summary <- avetempsfn(Cibola, "Cibola")
)

str(aveTemp)
aveTemp$site <- factor(aveTemp$site, levels = c("Delta", "StGeorge", "BigBend", "Cibola"),
                       labels = c("Delta" = "Delta (39°N)","StGeorge" = "St. George (37°N)", "BigBend" = "Big Bend (35°N)", "Cibola" = "Cibola (33°N)"))
aveTemp_longer <- aveTemp %>% pivot_longer(cols = ave_min:sd_max,
                               names_to = c(".value", "min_or_max"),
                               names_sep = "_")

#plot with ave max and min temps
ggplot(data = aveTemp) +
  geom_ribbon(aes(x = dayofyear, ymin = ave_min - sd_min, ymax = ave_min + sd_min), fill = 'deepskyblue') +
  geom_ribbon(aes(x = dayofyear, ymin = ave_max - sd_max, ymax = ave_max + sd_max), fill = 'indianred1') +
  # geom_ribbon(aes(x = dayofyear, ymin = 23, ymax = 38), color = "red", fill = "red", alpha = 0.2) +
  # geom_ribbon(aes(x = dayofyear, ymin = 13, ymax = 28), color = "blue", fill = "blue", alpha = 0.2) +
  facet_grid( ~ site) +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        strip.background = element_blank()) +
  scale_x_continuous(name = "Month", breaks = month.breaks, labels = month.letters) +
  scale_y_continuous(name = "Average High and Low Temperature (C), 2012-2021")

# dataAllSites_longer <- dataAllSites %>% pivot_longer(cols = c("TMAX", "TMIN"), names_to = "min_or_max", values_to = "TEMP")
#
# #Using raw data
# ggplot(data = dataAllSites_longer, aes(x = dayofyear, y = TEMP, group = NAME_group, color = min_or_max)) +
#   geom_line(alpha = 0.2) +
#   facet_grid( ~ NAME_group) +
#   theme_bw()

# treatline <- data.frame(Treatment = c("Warm", "Cool"), y = c(38, 28))
treatline.annotation <- data.frame(site = rep("Delta (39°N)",4),
                                   lab = c("Warm treat. high", "Cool treat. high", "Warm treat. low", "Cool treat.\n   low"),
                                   y = c(38, 28, 23, 12))
treatline.annotation$site <- factor(treatline.annotation$site,
                                    levels = c("Delta (39°N)", "St. George (37°N)", "Big Bend (35°N)", "Cibola (33°N)"))

#using longer form data for a better graph
site_temps <- ggplot(data = aveTemp_longer, aes(x = dayofyear, y = ave, color = min_or_max)) +

  #horizontal 0 line
  geom_hline(yintercept = c(0), color = 'gray50', size = 0.8) +

  #temp lines
  geom_ribbon(aes(ymax = ave + sd, ymin = ave-sd, fill = min_or_max), color = NA, alpha = 0.4) +
  # scale_fill_manual(name = "", labels = c("Daily Max", "Daily Min"), values = c("#CC79A7", "#56B4E9")) +
  scale_fill_manual(name = "", labels = c("Daily Max", "Daily Min"), values = c("grey70", "grey30")) +


  geom_line() +
  scale_color_manual(name = "", labels = c("Daily Max", "Daily Min"), values = c("grey70", "grey30")) +

  facet_grid(~ site) +

  #treatment reference lines
  geom_hline(yintercept = c(38), color = '#E69F00', size = 1.2) +
  geom_hline(yintercept = c(23), color = '#E69F00', size = 1) +

  geom_hline(yintercept = c(28), color = '#0072B2', size = 1.2) +
  geom_hline(yintercept = c(13), color = '#0072B2', size = 1) +

  geom_text(data = treatline.annotation, aes(x = 1, y = y, group = site, label = lab),
            hjust = "left", size = 3.4, inherit.aes = F, nudge_y = 1.35) +
  #formatting
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank()) +
  scale_x_continuous(name = "Month", breaks = month.breaks, labels = month.letters) +
  scale_y_continuous(name = "Temperature (°C)")
site_temps

#####
#Put map and temp plots together
ggarrange(full_map, site_temps, ncol = 2, labels = "AUTO", widths = c(1,3), legend = 'bottom')
##Export: 956 x 556 px (It looks horrible until it's at this size, so don't edit without adjusting size!!)

#Daylengths ----
##Daylength caclulations
oneyr = seq(1,365)

lat.delta = 39.376
lat.stgeorge = 37.116
lat.bigbend = 35.124
lat.cibola = 33.392

daylength.all <- rbind(
Delta.daylengths <- daylengthfn(lat.delta, oneyr, "Delta"),
StGeorge.daylengths <- daylengthfn(lat.stgeorge, oneyr, "StGeorge"),
BigBend.daylengths <- daylengthfn(lat.bigbend, oneyr, "BigBend"),
Cibola.daylengths <- daylengthfn(lat.cibola, oneyr, "Cibola")
)
daylength.all$site <- factor(daylength.all$site, levels = c("Delta", "StGeorge", "BigBend", "Cibola"),
                             labels = c("Delta" = "Delta (39°N)","StGeorge" = "St. George (37°N)", "BigBend" = "Big Bend (35°N)", "Cibola" = "Cibola (33°N)"))


# Transformation formula for 2 axes
ylim.prim = range(aveTemp_longer$ave)

ylim.prim = range(StGeorge.summary$ave_min)
ylim.sec = range(StGeorge.daylengths$daylength)
b.scale <- diff(ylim.prim)/diff(ylim.sec)
a.scale <- ylim.prim[1] - b.scale*ylim.sec[1]

#Plot with temp and daylength, two axes
ggplot(data = aveTemp) +
  geom_ribbon(aes(x = dayofyear, ymin = ave_min - sd_min, ymax = ave_min + sd_min), fill = 'skyblue') +
  geom_ribbon(aes(x = dayofyear, ymin = ave_max - sd_max, ymax = ave_max + sd_max), fill = 'pink') +
  geom_line(data = daylength.all, aes(x = day.new, y = a.scale + daylength*b.scale), color = 'gold', size = 3) +
  facet_grid( ~ site) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        strip.background = element_blank()) +
  scale_x_continuous(name = "Month", breaks = month.breaks, labels = month.letters) +
  scale_y_continuous(name = "Average High and Low Temperature (C), 2012-2021",
                     sec.axis = sec_axis(~ (. - a.scale)/b.scale, name = 'Daylength (hours)'))

# Correlations ----
##Correlations between daylength and average temperatures

daylength.temp <- right_join(x = aveTemp, y = daylength.all, by = c("dayofyear" = "day.new", "site" = "site"))
daylength.temp.delta <- daylength.temp %>% filter(site == "Delta (39°N)")
daylength.temp.stgeorge <- daylength.temp %>% filter(site == "St. George (37°N)")
daylength.temp.bigbend <- daylength.temp %>% filter(site == "Big Bend (35°N)")
daylength.temp.cibola <- daylength.temp %>% filter(site == "Cibola (33°N)")


#correlation tests
#max temp Pearson Method
cor.test(~ daylength + ave_max, data = daylength.temp, subset = site == "Delta (39°N)")
cor.test(~ daylength + ave_max, data = daylength.temp, subset = site == "St. George (37°N)")
cor.test(~ daylength + ave_max, data = daylength.temp, subset = site == "Big Bend (35°N)")
cor.test(~ daylength + ave_max, data = daylength.temp, subset = site == "Cibola (33°N)")

#min temp Pearson Method
cor.test(~ daylength + ave_min, data = daylength.temp, subset = site == "Delta (39°N)")
cor.test(~ daylength + ave_min, data = daylength.temp, subset = site == "St. George (37°N)")
cor.test(~ daylength + ave_min, data = daylength.temp, subset = site == "Big Bend (35°N)")
cor.test(~ daylength + ave_min, data = daylength.temp, subset = site == "Cibola (33°N)")

#max temp Kendall Method
cor.test(~ daylength + ave_max, data = daylength.temp, subset = site == "Delta (39°N)", method = "kendall")
cor.test(~ daylength + ave_max, data = daylength.temp, subset = site == "St. George (37°N)", method = "kendall")
cor.test(~ daylength + ave_max, data = daylength.temp, subset = site == "Big Bend (35°N)", method = "kendall")
cor.test(~ daylength + ave_max, data = daylength.temp, subset = site == "Cibola (33°N)", method = "kendall")

#min temp Kendall Method
cor.test(~ daylength + ave_min, data = daylength.temp, subset = site == "Delta (39°N)", method = "kendall")
cor.test(~ daylength + ave_min, data = daylength.temp, subset = site == "St. George (37°N)", method = "kendall")
cor.test(~ daylength + ave_min, data = daylength.temp, subset = site == "Big Bend (35°N)", method = "kendall")
cor.test(~ daylength + ave_min, data = daylength.temp, subset = site == "Cibola (33°N)", method = "kendall")

# bootstrap correlations MAX TEMP
ITER = 10000
method.choice <- "kendall" ## or "pearson"
boot.cor <- data.frame(delta = rep(0, ITER), stgeorge = rep(0, ITER), bigbend = rep(0, ITER), cibola = rep(0, ITER))
for (i in 1:ITER) {
  samp.d <- slice_sample(daylength.temp.delta, n = 365, replace = TRUE)
  boot.cor[i,1] <- cor(samp.d$ave_max, samp.d$daylength, method = method.choice)
  samp.s <- slice_sample(daylength.temp.stgeorge, n = 365, replace = TRUE)
  boot.cor[i,2] <- cor(samp.s$ave_max, samp.s$daylength, method = method.choice)
  samp.b <- slice_sample(daylength.temp.bigbend, n = 365, replace = TRUE)
  boot.cor[i,3] <- cor(samp.b$ave_max, samp.b$daylength, method = method.choice)
  samp.c <- slice_sample(daylength.temp.cibola, n = 365, replace = TRUE)
  boot.cor[i,4] <- cor(samp.c$ave_max, samp.c$daylength, method = method.choice)
}

boot.cor.long <- boot.cor %>% pivot_longer(cols = 1:4, names_to = "site", values_to = "correlation")
boot.cor.long %>% group_by(site) %>%
  summarise(
    mean = mean(correlation),
    sd = sd(correlation),
    CI.low = quantile(correlation, c(0.025)),
    CI.high = quantile(correlation, c(0.975))
  )
max.cor.plot <- ggplot(boot.cor.long, aes(x = correlation, color = site)) +
  geom_density() +
  scale_x_continuous(name = "Correlation between average daily maximum temperature and daylength") +
  theme_bw()
max.cor.plot

# bootstrap correlations MIN TEMP
boot.cor.min <- data.frame(delta = rep(0, ITER), stgeorge = rep(0, ITER), bigbend = rep(0, ITER), cibola = rep(0, ITER))
for (i in 1:ITER) {
  samp.d <- slice_sample(daylength.temp.delta, n = 365, replace = TRUE)
  boot.cor.min[i,1] <- cor(samp.d$ave_min, samp.d$daylength, method = method.choice)
  samp.s <- slice_sample(daylength.temp.stgeorge, n = 365, replace = TRUE)
  boot.cor.min[i,2] <- cor(samp.s$ave_min, samp.s$daylength, method = method.choice)
  samp.b <- slice_sample(daylength.temp.bigbend, n = 365, replace = TRUE)
  boot.cor.min[i,3] <- cor(samp.b$ave_min, samp.b$daylength, method = method.choice)
  samp.c <- slice_sample(daylength.temp.cibola, n = 365, replace = TRUE)
  boot.cor.min[i,4] <- cor(samp.c$ave_min, samp.c$daylength, method = method.choice)
}

boot.cor.min.long <- boot.cor.min %>% pivot_longer(cols = 1:4, names_to = "site", values_to = "correlation")
boot.cor.min.long %>% group_by(site) %>%
  summarise(
    mean = mean(correlation),
    sd = sd(correlation),
    CI.low = quantile(correlation, c(0.025)),
    CI.high = quantile(correlation, c(0.975))
  )
min.cor.plot <- ggplot(boot.cor.min.long, aes(x = correlation, color = site)) +
  geom_density() +
  scale_x_continuous(name = "Correlation between average daily minimum temperature and daylength") +
  theme_bw()
min.cor.plot

ggarrange(max.cor.plot, min.cor.plot, nrow = 2, ncol = 1, common.legend = T, legend = "bottom", labels = "AUTO")

# ggplot(data = Delta.summary, aes(x = dayofyear, y = ave_min)) +
#   geom_ribbon(aes(x = dayofyear, ymin = ave_min - sd_min, ymax = ave_min + sd_min), fill = 'skyblue') +
#   geom_line(data = Delta.daylengths, aes(x = day.new, y = a.scale + daylength*b.scale), color = 'gold', size = 4) +
#   scale_x_continuous(name = "Month", breaks = month.breaks, labels = month.letters) +
#   scale_y_continuous(name = "Average Low Temperature (°C) \n2010-2019", sec.axis = sec_axis(~ (. - a.scale)/b.scale, name = 'Daylength (hours)')) +
#   theme_classic() +
#   theme(plot.margin=unit(c(0.5,0.5,0.2,0.5), "cm"),
#         text = element_text(size = 20),
#         axis.title.y.left = element_text(color = 'skyblue', vjust = 2),
#         axis.title.y.right = element_text(color = 'gold', vjust = 3))
#
# temp.daylengthfn <- function(tempdata, daylengthdata){
#   # ylim.prim = range(tempdata$ave_min)
#   ylim.prim = c(-20,40)
#   # ylim.sec = range(daylengthdata$daylength)
#   ylim.sec = c(6,18)
#   b.scale <- diff(ylim.prim)/diff(ylim.sec)
#   a.scale <- ylim.prim[1] - b.scale*ylim.sec[1]
#
#   ggplot(data = tempdata, aes(x = dayofyear, y = ave_min)) +
#     geom_ribbon(aes(x = dayofyear, ymin = ave_min - sd_min, ymax = ave_min + sd_min), fill = 'skyblue') +
#     geom_line(data = daylengthdata, aes(x = day.new, y = a.scale + daylength*b.scale), color = 'gold', size = 4) +
#     scale_x_continuous(name = "Month", breaks = month.breaks, labels = month.letters) +
#     scale_y_continuous(limits = c(-20,40), name = "Average Low Temperature (°C) \n2010-2019",
#                        sec.axis = sec_axis(~ (. - a.scale)/b.scale, name = 'Daylength (hours)')) +
#     theme_classic() +
#     theme(plot.margin=unit(c(0.5,0.5,0.2,0.5), "cm"),
#           text = element_text(size = 20),
#           axis.title.y.left = element_text(color = 'skyblue', vjust = 2),
#           axis.title.y.right = element_text(color = 'gold', vjust = 3))
# }
#
# ggarrange(
# temp.daylengthfn(Delta.summary, Delta.daylengths),
# temp.daylengthfn(StGeorge.summary, StGeorge.daylengths),
# temp.daylengthfn(BigBend.summary, BigBend.daylengths),
# temp.daylengthfn(Cibola.summary, Cibola.daylengths),
# nrow = 1)

###########
###########
###########
#New Average Temperature plot ----

uspest.longertemp <- uspest %>% pivot_longer(cols = c(minimum.temperature.C, maximum.temperature.C),
                                             names_to = "min.max", values_to = "temperature") %>%
  group_by(population, week, min.max) %>%
  summarize(
    ave.temp = mean(temperature),
    max.temp = max(temperature),
    min.temp = min(temperature)
  )

ggplot() +
  geom_line(data = uspest.longertemp, aes(x = week, y = ave.temp, group = min.max, color = min.max)) +
  # scale_color_manual(values = colorRampPalette(brewer.pal(8, "Blues"))(11) ) +
  geom_line(data = uspest.longertemp %>% filter(min.max == "maximum.temperature.C"), aes(x = week, y = max.temp, group = min.max, color = min.max)) +
  geom_line(data = uspest.longertemp %>% filter(min.max == "minimum.temperature.C"), aes(x = week, y = min.temp, group = min.max, color = min.max)) +

  facet_wrap(~ population, ncol = 4)
















