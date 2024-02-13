#Load temperature data ----

dataAllSites <- read.csv("3378358_downloaded_6_29_23.csv")

month.breaks = c(1,32,60,91,121,152,182,213,244,274,305,335)+15
month.letters = c('J','F','M','A','M','J','J','A','S','O','N','D')

#Preliminary Clean up ----

prep1 <- function(data){
  if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
  if (!require('lubridate')) install.packages('lubridate'); library('lubridate')
  data$DATE <- ymd(data$DATE)
  data$dayofyear <- yday(data$DATE)
  data$YEAR <- year(data$DATE)
  data$MONTH <- month(data$DATE)
  data$DAY <- day(data$DATE)
  return(data)
}

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


#Calculate average temperatures ----
avetempsfn <- function(data, site){
  tempout <- data %>% group_by(dayofyear) %>%
    summarise(ave_min = mean(TMIN, na.rm = T),
              sd_min = sd(TMIN, na.rm = T),
              ave_max = mean(TMAX, na.rm = T),
              sd_max = sd(TMAX, na.rm = T))
  tempout2 <- cbind(site = rep(site, (nrow(tempout))), tempout)
  return(tempout2)
}

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

#Treatment values from experiment ----
treatline.annotation <- data.frame(site = rep("Delta (39°N)",4),
                                   lab = c("Warm treat. high", "Cool treat. high", "Warm treat. low", "Cool treat.\n   low"),
                                   y = c(38, 28, 23, 12))
treatline.annotation$site <- factor(treatline.annotation$site,
                                    levels = c("Delta (39°N)", "St. George (37°N)", "Big Bend (35°N)", "Cibola (33°N)"))

#Fig - Temperature at each site ----
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
#Export: 875 x 590
